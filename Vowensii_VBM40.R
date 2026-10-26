### Packages Pathogen survival ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(broom)
library(ggpubr)
library(survival)
library(coxme)
### Theme set ####
theme_set(theme_bw())
## Make the degree symbol render in the figures
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

### Import and set everything ####
## Working directory
setwd("pathogen results")
## Folder for the outputs
out <- "FINAL_FINAL_BALANCED/v4_outputs/"
dir.create(out, showWarnings = FALSE)

## Import mortality datasheet
mortality <- read.csv("data_mortality.csv", header = TRUE, fileEncoding = "UTF-8-BOM",
                      colClasses = c(Fragment = "character"))
head(mortality)

## Transform dates (dd/mm/yy)
mortality$Start <- as.Date(mortality$Start, format = "%d/%m/%y")
mortality$End <- as.Date(mortality$End, format = "%d/%m/%y")
mortality$Mortality <- as.Date(mortality$Mortality, format = "%d/%m/%y")

## Order variable levels (by default is alphabetic order)
mortality$Treatment <- ifelse(mortality$Treatment == "Control", "Control", "VBM40")
mortality$Treatment <- ordered(mortality$Treatment, levels = c("Control", "VBM40"))
levels(mortality$Treatment)

mortality$Temperature <- ordered(mortality$Temperature, levels = c("25", "33"), labels = c("25C", "33C"))
levels(mortality$Temperature)

mortality$Colony <- factor(mortality$Colony)
mortality$Tank <- factor(paste0("Tank", sprintf("%02d", mortality$Tank)))

## Dead = mortality date on or before the end of the experiment; otherwise alive (right-censored at the end)
## time = days from Start to death (dead) or to End (alive)
mortality <- mortality %>%
  mutate(event = as.integer(!is.na(Mortality) & Mortality <= End),
         time = ifelse(event == 1, as.numeric(Mortality - Start), as.numeric(End - Start)),
         Group = factor(paste(Treatment, Temperature),
                        levels = c("Control 25C", "Control 33C", "VBM40 25C", "VBM40 33C")))
end_day <- max(as.numeric(mortality$End - mortality$Start))
end_day

## Fragments per tank
table(mortality$Tank, mortality$Group)
write.csv(mortality, file = paste0(out, "Data_time_to_event.csv"), row.names = FALSE)

### Mortality summary ####
## Proportion of dead fragments with exact (Clopper-Pearson) 95% CI
mortality_summary <- mortality %>%
  group_by(Temperature, Treatment, Group) %>%
  summarise(n = n(), dead = sum(event), .groups = "drop") %>%
  mutate(prop_dead = dead / n,
         ci_lower = map2_dbl(dead, n, ~ binom.test(.x, .y)$conf.int[1]),
         ci_upper = map2_dbl(dead, n, ~ binom.test(.x, .y)$conf.int[2]))
mortality_summary
write.csv(mortality_summary, file = paste0(out, "Mortality_summary.csv"), row.names = FALSE)

## Proportion of dead fragments per tank
mortality_tank <- mortality %>%
  group_by(Tank, Treatment, Temperature, Group) %>%
  summarise(n = n(), dead = sum(event), .groups = "drop") %>%
  mutate(prop_dead = dead / n)
mortality_tank
write.csv(mortality_tank, file = paste0(out, "Mortality_per_tank.csv"), row.names = FALSE)

### Kaplan-Meier and log-rank test ####
## Kaplan-Meier estimates per group (median survival with 95% CI)
km_fit <- survfit(Surv(time, event) ~ Group, data = mortality)
km_fit
## First VBM40 death (day) per temperature
mortality %>%
  filter(Treatment == "VBM40", event == 1) %>%
  group_by(Temperature) %>%
  summarise(first_death = min(time))
## Survival at the end of the experiment
summary(km_fit, times = end_day)

## Log-rank test between Control and VBM40, within each temperature, stratified by donor colony
logrank_results <- map_dfr(levels(mortality$Temperature), function(tp) {
  sub_df <- mortality %>% filter(Temperature == tp)
  fit <- survdiff(Surv(time, event) ~ Treatment + strata(Colony), data = sub_df)
  data.frame(Temperature = tp,
             chisq = fit$chisq,
             df = 1,
             p.value = pchisq(fit$chisq, df = 1, lower.tail = FALSE))
})
logrank_results
write.csv(logrank_results, file = paste0(out, "Logrank_by_temperature.csv"), row.names = FALSE)

### Fisher's exact test (final mortality) ####
## Three pairwise comparisons: OR = second group vs first (reference)
fisher_comparisons <- list(
  "VBM40 vs Control, 25C" = list(data = mortality %>% filter(Temperature == "25C"),
                                 var = "Treatment", ref = "Control", alt = "VBM40"),
  "VBM40 vs Control, 33C" = list(data = mortality %>% filter(Temperature == "33C"),
                                 var = "Treatment", ref = "Control", alt = "VBM40"),
  "33C vs 25C, VBM40" = list(data = mortality %>% filter(Treatment == "VBM40"),
                             var = "Temperature", ref = "25C", alt = "33C"))

fisher_results <- map_dfr(names(fisher_comparisons), function(nm) {
  cp <- fisher_comparisons[[nm]]
  grp <- factor(as.character(cp$data[[cp$var]]), levels = c(cp$alt, cp$ref))
  tab <- table(Group = grp, Dead = factor(cp$data$event, levels = c(1, 0)))
  print(tab)
  fisher <- fisher.test(tab)
  data.frame(comparison = nm,
             dead_ref = tab[2, 1], n_ref = sum(tab[2, ]),
             dead_alt = tab[1, 1], n_alt = sum(tab[1, ]),
             OR = unname(fisher$estimate),
             OR_lower95 = fisher$conf.int[1],
             OR_upper95 = fisher$conf.int[2],
             p.value = fisher$p.value)
})
fisher_results
write.csv(fisher_results, file = paste0(out, "Fisher_exact_test.csv"), row.names = FALSE)

### Mixed-effects Cox model (temperature among VBM40) ####
## Temperature and donor colony as fixed effects, tank as random effect on the intercept
vbm40 <- mortality %>% filter(Treatment == "VBM40")
vbm40$Tank <- droplevels(vbm40$Tank)
vbm40$Temperature <- relevel(factor(vbm40$Temperature, ordered = FALSE), ref = "25C")

cox_model <- coxme(Surv(time, event) ~ Temperature + Colony + (1 | Tank), data = vbm40)
cox_model

cox_coef <- fixef(cox_model)
cox_se <- sqrt(diag(vcov(cox_model)))
cox_results <- data.frame(term = names(cox_coef),
                          HR = exp(cox_coef),
                          lower95 = exp(cox_coef - 1.96 * cox_se),
                          upper95 = exp(cox_coef + 1.96 * cox_se),
                          z = cox_coef / cox_se,
                          p.value = 2 * pnorm(-abs(cox_coef / cox_se)))
cox_results
## Standard deviation of the tank random effect
sqrt(VarCorr(cox_model)$Tank)
write.csv(cox_results, file = paste0(out, "Cox_mixed_effects_temperature.csv"), row.names = FALSE)

## Proportional hazards assumption (Schoenfeld residuals)
## cox.zph does not work with coxme, so it is run on the same model without the tank term
cox_fixed <- coxph(Surv(time, event) ~ Temperature + Colony, data = vbm40)
schoenfeld <- cox.zph(cox_fixed)
schoenfeld
write.csv(as.data.frame(schoenfeld$table), file = paste0(out, "Schoenfeld_PH_test.csv"))

### Figure: Kaplan-Meier (A) and mortality (B) ####
group_cols <- c("Control 25C" = "#3AAFE0", "Control 33C" = "#2B2B8C",
                "VBM40 25C" = "#F28E1C", "VBM40 33C" = "#E0262B")
group_labels <- c("Control 25C" = "Control 25 °C", "Control 33C" = "Control 33 °C",
                  "VBM40 25C" = "VBM40 25 °C", "VBM40 33C" = "VBM40 33 °C")

#PLOT_KAPLAN_MEIER
## Survival curve per group, from day 0 to the end of the experiment
km_df <- map_dfr(levels(mortality$Group), function(g) {
  fit <- survfit(Surv(time, event) ~ 1, data = mortality %>% filter(Group == g))
  data.frame(Group = g,
             time = c(0, fit$time, end_day),
             surv = c(1, fit$surv, tail(fit$surv, 1)))
})
km_df$Group <- factor(km_df$Group, levels = levels(mortality$Group))

km_plot <- ggplot(km_df, aes(x = time, y = surv, colour = Group)) +
  geom_step(linewidth = 1.1) +
  geom_vline(xintercept = end_day, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = group_cols, labels = group_labels) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(x = "Days", y = "Survival probability", colour = "Treatments") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 16),
                     axis.text.y = element_text(size = 14)) +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.background = element_blank(), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        plot.margin = margin(t = 12, r = 6, b = 6, l = 22),
        text = element_text(size = 14))
km_plot

#PLOT_MORTALITY
## Height of the dead/total label: above both the CI and the highest tank dot
label_height <- mortality_tank %>%
  group_by(Group) %>%
  summarise(max_tank = 100 * max(prop_dead), .groups = "drop")
mortality_summary <- mortality_summary %>%
  left_join(label_height, by = "Group") %>%
  mutate(label_y = pmax(100 * ci_upper, max_tank) + 4)

## Bars = proportion of dead fragments (exact 95% CI), dots = tanks
mortality_plot <- ggplot(mortality_summary, aes(x = Group, y = 100 * prop_dead, fill = Group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper), width = 0.2) +
  geom_text(aes(y = label_y, label = paste0(dead, "/", n)), size = 3.5, vjust = 0) +
  annotate("point", x = 0.6, y = 112, shape = 21, fill = "white", colour = "black", size = 2.5, stroke = 0.8) +
  annotate("text", x = 0.8, y = 112, label = "Tank", hjust = 0, size = 4) +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(labels = c("Control\n25 °C", "Control\n33 °C", "VBM40\n25 °C", "VBM40\n33 °C")) +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Mortality %") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 16),
                     axis.text.y = element_text(size = 14),
                     axis.text.x = element_text(size = 11)) +
  theme(legend.position = "none", text = element_text(size = 14),
        plot.margin = margin(t = 12, r = 6, b = 6, l = 22))
mortality_plot

#MERGE_PLOTS
Fig_KM_mortality <- ggarrange(km_plot, mortality_plot, ncol = 2, labels = c("A", "B"),
                              widths = c(1.6, 1), font.label = list(size = 18, face = "bold"))
Fig_KM_mortality

ggsave(paste0(out, "Fig_KM_mortality.png"), Fig_KM_mortality, width = 2800 / 300, height = 1400 / 300, dpi = 300)
## Vector PDF with editable text and objects
ggsave(paste0(out, "Fig_KM_mortality.pdf"), Fig_KM_mortality, width = 2800 / 300, height = 1400 / 300, device = cairo_pdf)

### Session info ####
sessionInfo()
