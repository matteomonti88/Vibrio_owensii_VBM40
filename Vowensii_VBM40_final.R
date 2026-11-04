### Packages Pathogen survival ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggpubr)
library(survival)
library(coxme)
library(lme4)
library(openxlsx)
### Theme set ####
theme_set(theme_bw())

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

### Import and set everything ####
## Working directory = folder with the datasheet
setwd("analysis")
## Folder for the outputs
out <- "outputs/"
dir.create(out, showWarnings = FALSE)

## Import mortality datasheet (now with the water bath column)
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
mortality$Waterbath <- factor(paste0("Bath", mortality$Waterbath))

## Dead = mortality date on or before the end of the experiment; otherwise alive (right-censored at the end)
## time = days from Start to death (dead) or to End (alive)
mortality <- mortality %>%
  mutate(event = as.integer(!is.na(Mortality) & Mortality <= End),
         time = ifelse(event == 1, as.numeric(Mortality - Start), as.numeric(End - Start)),
         Group = factor(paste(Treatment, Temperature),
                        levels = c("Control 25C", "Control 33C", "VBM40 25C", "VBM40 33C")))
end_day <- max(as.numeric(mortality$End - mortality$Start))
end_day

## Design: which tanks are in which water bath, temperature and treatment
design <- mortality %>%
  distinct(Waterbath, Tank, Temperature, Treatment) %>%
  arrange(Waterbath, Tank)
design
write.csv(design, file = paste0(out, "Design_waterbath_tank.csv"), row.names = FALSE)
write.csv(mortality, file = paste0(out, "Data_time_to_event.csv"), row.names = FALSE)

## Function to turn a model call into a single line of text
call_text <- function(x) paste(trimws(deparse(x, width.cutoff = 500)), collapse = " ")

### Mortality summary ####
## Proportion of dead fragments with exact 95% CI
mortality_summary <- mortality %>%
  group_by(Temperature, Treatment, Group) %>%
  summarise(n = n(), dead = sum(event), .groups = "drop") %>%
  mutate(prop_dead = dead / n,
         ci_lower = map2_dbl(dead, n, ~ binom.test(.x, .y)$conf.int[1]),
         ci_upper = map2_dbl(dead, n, ~ binom.test(.x, .y)$conf.int[2]))
mortality_summary
write.csv(mortality_summary, file = paste0(out, "Mortality_summary.csv"), row.names = FALSE)

## Proportion of dead fragments per tank (with its water bath)
mortality_tank <- mortality %>%
  group_by(Waterbath, Tank, Treatment, Temperature, Group) %>%
  summarise(n = n(), dead = sum(event), .groups = "drop") %>%
  mutate(prop_dead = dead / n)
mortality_tank
write.csv(mortality_tank, file = paste0(out, "Mortality_per_tank.csv"), row.names = FALSE)

## Mortality per aquarium by donor colony (dead/total)
mortality_tank_colony <- mortality %>%
  group_by(Tank, Colony) %>%
  summarise(dead_total = paste0(sum(event), "/", n()), .groups = "drop") %>%
  pivot_wider(names_from = Colony, values_from = dead_total, names_prefix = "Colony")
mortality_tank_colony <- mortality_tank %>%
  left_join(mortality_tank_colony, by = "Tank") %>%
  mutate(dead_total = paste0(dead, "/", n)) %>%
  arrange(Temperature, Treatment, Tank)
mortality_tank_colony
write.csv(mortality_tank_colony, file = paste0(out, "Mortality_per_aquarium_by_colony.csv"), row.names = FALSE)

### Kaplan-Meier and log-rank test ####
## Kaplan-Meier estimates per group (median survival with 95% CI)
km_fit <- survfit(Surv(time, event) ~ Group, data = mortality)
km_fit
## Survival at the end of the experiment
km_end <- summary(km_fit, times = end_day)
km_end

## Log-rank test between Control and VBM40, within each temperature, stratified by donor colony
logrank_fits <- map(levels(mortality$Temperature), function(tp) {
  sub_df <- mortality %>% filter(Temperature == tp)
  survdiff(Surv(time, event) ~ Treatment + strata(Colony), data = sub_df)
})
names(logrank_fits) <- levels(mortality$Temperature)
logrank_results <- map_dfr(names(logrank_fits), function(tp) {
  fit <- logrank_fits[[tp]]
  data.frame(Temperature = tp, chisq = fit$chisq, df = 1,
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
  fisher <- fisher.test(tab)
  data.frame(Comparison = nm,
             Reference_group = cp$ref, Reference_dead = tab[2, 1], Reference_alive = tab[2, 2],
             Test_group = cp$alt, Test_dead = tab[1, 1], Test_alive = tab[1, 2],
             Odds_ratio = unname(fisher$estimate),
             OR_lower95 = fisher$conf.int[1],
             OR_upper95 = fisher$conf.int[2],
             p.value = fisher$p.value)
})
fisher_results
write.csv(fisher_results, file = paste0(out, "Fisher_exact_test.csv"), row.names = FALSE)

### Mixed-effects Cox models (temperature among VBM40) ####
vbm40 <- mortality %>% filter(Treatment == "VBM40")
vbm40$Tank <- droplevels(vbm40$Tank)
vbm40$Waterbath <- droplevels(vbm40$Waterbath)
vbm40$Temperature <- relevel(factor(vbm40$Temperature, ordered = FALSE), ref = "25C")

## Model 1: tank as random intercept (as in the main analysis)
cox_tank <- coxme(Surv(time, event) ~ Temperature + Colony + (1 | Tank), data = vbm40)
cox_tank
## Model 2: tank nested within water bath (random intercepts for both)
cox_bath <- coxme(Surv(time, event) ~ Temperature + Colony + (1 | Waterbath/Tank), data = vbm40)
cox_bath

## Fixed effects of both models
cox_table <- function(model, label) {
  b <- fixef(model)
  se <- sqrt(diag(vcov(model)))
  data.frame(Model = label,
             Term = c("Temperature 33C vs 25C", "Colony 2 vs Colony 1", "Colony 3 vs Colony 1"),
             Coefficient = b, SE = se,
             Hazard_ratio = exp(b),
             HR_lower95 = exp(b - 1.96 * se),
             HR_upper95 = exp(b + 1.96 * se),
             z = b / se,
             p.value = 2 * pnorm(-abs(b / se)),
             row.names = NULL)
}
cox_results <- bind_rows(cox_table(cox_tank, "Tank random intercept"),
                         cox_table(cox_bath, "Water bath / tank random intercepts"))
cox_results
write.csv(cox_results, file = paste0(out, "Cox_mixed_effects_temperature.csv"), row.names = FALSE)

## Random effects and model fit
cox_info <- function(model, label) {
  vc <- VarCorr(model)
  data.frame(Model = label,
             Item = c("Fragments", "Deaths (events)",
                      paste0("Random-intercept variance: ", names(vc)),
                      paste0("Random-intercept SD: ", names(vc)),
                      "Log-likelihood, null model", "Log-likelihood, fitted model (integrated)"),
             Value = c(model$n[2], model$n[1],
                       sapply(vc, as.numeric), sqrt(sapply(vc, as.numeric)),
                       model$loglik[1], model$loglik[2]),
             row.names = NULL)
}
cox_random <- bind_rows(cox_info(cox_tank, "Tank random intercept"),
                        cox_info(cox_bath, "Water bath / tank random intercepts"))
cox_random
write.csv(cox_random, file = paste0(out, "Cox_random_effects_and_fit.csv"), row.names = FALSE)

## Does adding the water bath improve the model? (likelihood-ratio test)
cox_comparison <- anova(cox_tank, cox_bath)
cox_comparison
write.csv(as.data.frame(cox_comparison), file = paste0(out, "Cox_tank_vs_waterbath_LRT.csv"))

## Proportional hazards assumption (Schoenfeld residuals)
cox_fixed <- coxph(Surv(time, event) ~ Temperature + Colony, data = vbm40)
schoenfeld <- cox.zph(cox_fixed)
schoenfeld
write.csv(as.data.frame(schoenfeld$table), file = paste0(out, "Schoenfeld_PH_test.csv"))

### Aquarium-level binomial mixed model (temperature among VBM40) ####
## One count per aquarium and donor colony: number of dead and surviving fragments
aquarium_counts <- mortality %>%
  filter(Treatment == "VBM40") %>%
  group_by(Waterbath, Tank, Temperature, Colony) %>%
  summarise(dead = sum(event), alive = n() - sum(event), .groups = "drop")
aquarium_counts$Temperature <- relevel(factor(aquarium_counts$Temperature, ordered = FALSE), ref = "25C")
aquarium_counts
write.csv(aquarium_counts, file = paste0(out, "Aquarium_counts_VBM40.csv"), row.names = FALSE)

## Model 1: aquarium nested within water bath as random intercepts
glmm_bath <- glmer(cbind(dead, alive) ~ Temperature + Colony + (1 | Waterbath/Tank),
                   data = aquarium_counts, family = binomial)
## Model 2: aquarium as random intercept only
glmm_tank <- glmer(cbind(dead, alive) ~ Temperature + Colony + (1 | Tank),
                   data = aquarium_counts, family = binomial)
summary(glmm_bath)

glmm_table <- function(model, label) {
  cf <- summary(model)$coefficients
  data.frame(Model = label,
             Term = c("Temperature 33C vs 25C", "Colony 2 vs Colony 1", "Colony 3 vs Colony 1"),
             Coefficient = cf[-1, 1], SE = cf[-1, 2],
             Odds_ratio = exp(cf[-1, 1]),
             OR_lower95 = exp(cf[-1, 1] - 1.96 * cf[-1, 2]),
             OR_upper95 = exp(cf[-1, 1] + 1.96 * cf[-1, 2]),
             z = cf[-1, 3], p.value = cf[-1, 4], row.names = NULL)
}
glmm_results <- bind_rows(glmm_table(glmm_bath, "Aquarium nested in water bath random intercepts"),
                          glmm_table(glmm_tank, "Aquarium random intercept"))
glmm_results
write.csv(glmm_results, file = paste0(out, "GLMM_binomial_temperature.csv"), row.names = FALSE)

## Random effects 
glmm_random <- as.data.frame(VarCorr(glmm_bath))[, c("grp", "vcov", "sdcor")]
colnames(glmm_random) <- c("Random_effect", "Variance", "SD")
glmm_random
write.csv(glmm_random, file = paste0(out, "GLMM_binomial_random_effects.csv"), row.names = FALSE)
## Does adding the water bath improve the model? (likelihood-ratio test)
glmm_comparison <- anova(glmm_tank, glmm_bath)
glmm_comparison
write.csv(as.data.frame(glmm_comparison), file = paste0(out, "GLMM_tank_vs_waterbath_LRT.csv"))

### Tank-level exact permutation tests ####
## The tank is the unit that received the treatment: one value per tank (proportion of dead fragments)
## Statistic = difference in mean tank mortality; all possible assignments of the labels to the tanks are enumerated
perm_one <- function(tank_df, group, alt) {
  lab <- tank_df[[group]] == alt
  obs <- mean(tank_df$prop_dead[lab]) - mean(tank_df$prop_dead[!lab])
  null <- apply(combn(nrow(tank_df), sum(lab)), 2,
                function(ix) mean(tank_df$prop_dead[ix]) - mean(tank_df$prop_dead[-ix]))
  list(obs = obs, null = null)
}
p_perm <- function(x) mean(abs(x$null) >= abs(x$obs) - 1e-12)

perm25 <- perm_one(mortality_tank %>% filter(Temperature == "25C"), "Treatment", "VBM40")
perm33 <- perm_one(mortality_tank %>% filter(Temperature == "33C"), "Treatment", "VBM40")
permT <- perm_one(mortality_tank %>% filter(Treatment == "VBM40"), "Temperature", "33C")
## VBM40 effect across both temperatures (labels permuted within temperature)
null_both <- rowSums(expand.grid(perm25$null, perm33$null))
obs_both <- perm25$obs + perm33$obs

permutation_results <- data.frame(
  Comparison = c("VBM40 vs Control, 25C", "VBM40 vs Control, 33C",
                 "33C vs 25C, VBM40", "VBM40 vs Control, both temperatures (stratified)"),
  Observed_difference = c(perm25$obs, perm33$obs, permT$obs, obs_both),
  Allocations = c(length(perm25$null), length(perm33$null), length(permT$null), length(null_both)),
  p.value = c(p_perm(perm25), p_perm(perm33), p_perm(permT),
              mean(abs(null_both) >= abs(obs_both) - 1e-12)))
permutation_results
write.csv(permutation_results, file = paste0(out, "Tank_permutation_tests.csv"), row.names = FALSE)


### Figure: Kaplan-Meier (A) and mortality (B) ####
group_cols <- c("Control 25C" = "#3AAFE0", "Control 33C" = "#2B2B8C",
                "VBM40 25C" = "#F28E1C", "VBM40 33C" = "#E0262B")
group_labels <- c("Control 25C" = "Control 25 °C", "Control 33C" = "Control 33 °C",
                  "VBM40 25C" = "VBM40 25 °C", "VBM40 33C" = "VBM40 33 °C")

#PLOT_KAPLAN_MEIER
## Survival curve per group, from day 0 to the end of the experiment
km_df <- map_dfr(levels(mortality$Group), function(g) {
  fit <- survfit(Surv(time, event) ~ 1, data = mortality %>% filter(Group == g))
  data.frame(Group = g, time = c(0, fit$time, end_day), surv = c(1, fit$surv, tail(fit$surv, 1)))
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
mortality_label <- mortality_summary %>%
  left_join(label_height, by = "Group") %>%
  mutate(label_y = pmax(100 * ci_upper, max_tank) + 4)

## Bars = proportion of dead fragments (exact 95% CI), dots = tanks (shape = water bath)
mortality_plot <- ggplot(mortality_label, aes(x = Group, y = 100 * prop_dead, fill = Group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = 100 * ci_lower, ymax = 100 * ci_upper), width = 0.2) +
  geom_point(data = mortality_tank, aes(x = Group, y = 100 * prop_dead, shape = Waterbath),
             fill = "white", colour = "black", size = 2.8, stroke = 0.8,
             position = position_jitter(width = 0.12, height = 0, seed = 123),
             inherit.aes = FALSE) +
  geom_text(aes(y = label_y, label = paste0(dead, "/", n)), size = 3.5, vjust = 0) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Water bath") +
  scale_x_discrete(labels = c("Control\n25 °C", "Control\n33 °C", "VBM40\n25 °C", "VBM40\n33 °C")) +
  scale_y_continuous(limits = c(0, 128), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Mortality %") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size = 16),
                     axis.text.y = element_text(size = 14),
                     axis.text.x = element_text(size = 11)) +
  theme(legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
        legend.direction = "vertical", legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0, "cm"),
        legend.title = element_text(size = 10), legend.text = element_text(size = 10),
        text = element_text(size = 14),
        plot.margin = margin(t = 12, r = 6, b = 6, l = 22))
mortality_plot

#MERGE_PLOTS
Fig_KM_mortality <- ggarrange(km_plot, mortality_plot, ncol = 2, labels = c("A", "B"),
                              widths = c(1.6, 1), font.label = list(size = 18, face = "bold"))
Fig_KM_mortality

ggsave(paste0(out, "Fig_KM_mortality.pdf"), Fig_KM_mortality, width = 2800 / 300, height = 1400 / 300, device = cairo_pdf)
