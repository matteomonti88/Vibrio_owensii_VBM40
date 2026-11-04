### Packages Pathogen survival (with water bath) ####
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
## Make the degree symbol render in the figures
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

### Import and set everything ####
## Working directory = folder with the datasheet
setwd("pathogen results/Waterbath_analysis_v2")
## Folder for the outputs
out <- "outputs/"
dir.create(out, showWarnings = FALSE)

## Import mortality datasheet (now with the water bath column)
mortality <- read.csv("data_mortality_waterbath.csv", header = TRUE, fileEncoding = "UTF-8-BOM",
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
## Proportion of dead fragments with exact (Clopper-Pearson) 95% CI
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

## Random effects (variance and SD) and model fit
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
## cox.zph does not work with coxme, so it is run on the same model without the random effects
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

## Random effects (variance and SD)
glmm_random <- as.data.frame(VarCorr(glmm_bath))[, c("grp", "vcov", "sdcor")]
colnames(glmm_random) <- c("Random_effect", "Variance", "SD")
glmm_random
write.csv(glmm_random, file = paste0(out, "GLMM_binomial_random_effects.csv"), row.names = FALSE)
## Does adding the water bath improve the model? (likelihood-ratio test)
glmm_comparison <- anova(glmm_tank, glmm_bath)
glmm_comparison
write.csv(as.data.frame(glmm_comparison), file = paste0(out, "GLMM_tank_vs_waterbath_LRT.csv"))

## Strictest version: one value per aquarium, water bath as a fixed block, quasi-binomial
aquarium_totals <- aquarium_counts %>%
  group_by(Waterbath, Tank, Temperature) %>%
  summarise(dead = sum(dead), alive = sum(alive), .groups = "drop")
glm_block <- glm(cbind(dead, alive) ~ Waterbath + Temperature, data = aquarium_totals, family = quasibinomial)
summary(glm_block)
cf_block <- summary(glm_block)$coefficients["Temperature33C", ]
t_crit <- qt(0.975, glm_block$df.residual)
block_results <- data.frame(Model = "Aquarium totals, water bath as fixed block (quasi-binomial)",
                            Term = "Temperature 33C vs 25C",
                            Odds_ratio = exp(cf_block[1]),
                            OR_lower95 = exp(cf_block[1] - t_crit * cf_block[2]),
                            OR_upper95 = exp(cf_block[1] + t_crit * cf_block[2]),
                            t = cf_block[3], df = glm_block$df.residual, p.value = cf_block[4],
                            row.names = NULL)
block_results
write.csv(block_results, file = paste0(out, "GLM_block_temperature.csv"), row.names = FALSE)

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

### Water bath-restricted comparisons ####
## Every water bath holds tanks from both temperatures, so comparisons can be made inside each bath.
## Within each bath the difference in tank mortality is calculated, and the sign of each difference
## is flipped in all possible ways (exact paired permutation / sign-flip test)
signflip <- function(d) {
  obs <- mean(d)
  null <- rowMeans(sweep(as.matrix(expand.grid(rep(list(c(-1, 1)), length(d)))), 2, d, "*"))
  data.frame(Pairs = length(d), Observed_mean_difference = obs,
             Allocations = length(null), p.value = mean(abs(null) >= abs(obs) - 1e-12))
}

## Temperature effect among VBM40 tanks: 33C minus 25C within each bath
bath_temp <- mortality_tank %>%
  filter(Treatment == "VBM40") %>%
  select(Waterbath, Temperature, prop_dead) %>%
  pivot_wider(names_from = Temperature, values_from = prop_dead) %>%
  mutate(difference = `33C` - `25C`)
bath_temp
## VBM40 effect: VBM40 minus Control within each bath and temperature that contains both
bath_treat <- mortality_tank %>%
  group_by(Waterbath, Temperature, Treatment) %>%
  summarise(prop_dead = mean(prop_dead), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = prop_dead) %>%
  filter(!is.na(Control), !is.na(VBM40)) %>%
  mutate(difference = VBM40 - Control)
bath_treat

bath_tests <- bind_rows(
  cbind(Comparison = "33C vs 25C, VBM40 (difference within each water bath)", signflip(bath_temp$difference)),
  cbind(Comparison = "VBM40 vs Control (difference within each water bath and temperature)", signflip(bath_treat$difference)))
bath_tests
## Paired t-test on aquarium mortality across the four water baths (33C vs 25C)
paired_t <- t.test(bath_temp$`33C`, bath_temp$`25C`, paired = TRUE)
paired_t
paired_t_results <- data.frame(Comparison = "33C vs 25C, VBM40 (paired t-test across water baths)",
                               Mean_difference = unname(paired_t$estimate), t = unname(paired_t$statistic),
                               df = unname(paired_t$parameter), CI_lower95 = paired_t$conf.int[1],
                               CI_upper95 = paired_t$conf.int[2], p.value = paired_t$p.value)
write.csv(paired_t_results, file = paste0(out, "Waterbath_paired_t_test.csv"), row.names = FALSE)
write.csv(bath_temp, file = paste0(out, "Waterbath_pairs_temperature.csv"), row.names = FALSE)
write.csv(bath_treat, file = paste0(out, "Waterbath_pairs_treatment.csv"), row.names = FALSE)
write.csv(bath_tests, file = paste0(out, "Waterbath_signflip_tests.csv"), row.names = FALSE)

### Tables for the Supplementary materials ####
km_table <- as.data.frame(summary(km_fit)$table)
km_table$Group <- gsub("Group=", "", rownames(km_table))
km_end_df <- data.frame(Group = gsub("Group=", "", as.character(km_end$strata)),
                        surv_end = km_end$surv, surv_end_lower95 = km_end$lower, surv_end_upper95 = km_end$upper)
Table_km <- km_table %>%
  select(Group, Fragments = records, Deaths = events, Median_survival_days = median,
         Median_lower95 = `0.95LCL`, Median_upper95 = `0.95UCL`) %>%
  left_join(km_end_df, by = "Group") %>%
  left_join(mortality %>% filter(event == 1) %>% group_by(Group = as.character(Group)) %>%
              summarise(First_death_day = min(time), .groups = "drop"), by = "Group")
Table_km
write.csv(Table_km, file = paste0(out, "Table_Kaplan_Meier.csv"), row.names = FALSE)

Table_logrank <- map_dfr(names(logrank_fits), function(tp) {
  fit <- logrank_fits[[tp]]
  obs <- if (is.matrix(fit$obs)) rowSums(fit$obs) else fit$obs
  exp <- if (is.matrix(fit$exp)) rowSums(fit$exp) else fit$exp
  data.frame(Temperature = tp, Group = gsub("Treatment=", "", names(fit$n)),
             Fragments = as.integer(fit$n), Observed_deaths = obs, Expected_deaths = exp,
             Chisq = fit$chisq, df = 1, p.value = pchisq(fit$chisq, df = 1, lower.tail = FALSE))
})
Table_logrank
write.csv(Table_logrank, file = paste0(out, "Table_log_rank.csv"), row.names = FALSE)

Table_models <- as.data.frame(tibble::tribble(
  ~Analysis, ~Question, ~`R function (package)`, ~`Model formula / call`, ~`Data used`, ~Notes,
  "Kaplan-Meier survival curves", "Survival over time in each treatment x temperature group",
  "survfit (survival)", call_text(km_fit$call), paste0("All fragments (n = ", nrow(mortality), ")"),
  "Fragments alive at the end are right-censored",
  "Log-rank test, 25 C", "Control vs VBM40 survival at 25 C", "survdiff (survival)",
  call_text(logrank_fits[["25C"]]$call), "Fragments at 25 C (sub_df)", "Stratified by donor colony",
  "Log-rank test, 33 C", "Control vs VBM40 survival at 33 C", "survdiff (survival)",
  call_text(logrank_fits[["33C"]]$call), "Fragments at 33 C (sub_df)", "Stratified by donor colony",
  "Fisher's exact test (3 comparisons)",
  "Final mortality: (i) VBM40 vs Control at 25 C, (ii) VBM40 vs Control at 33 C, (iii) 33 C vs 25 C in VBM40",
  "fisher.test (stats)", "fisher.test(table(Group, Dead))", "2 x 2 table of the fragments in each comparison",
  "Odds ratio of the test group relative to the reference group",
  "Mixed-effects Cox model, tank random effect", "Effect of temperature on the hazard of death in VBM40-exposed fragments",
  "coxme (coxme)", call_text(cox_tank$call), paste0("VBM40-exposed fragments (n = ", nrow(vbm40), ")"),
  "Fixed effects: temperature (reference 25C) and donor colony (reference colony 1). Random intercept per tank",
  "Mixed-effects Cox model, water bath and tank random effects", "As above, adding the water bath",
  "coxme (coxme)", call_text(cox_bath$call), paste0("VBM40-exposed fragments (n = ", nrow(vbm40), ")"),
  "Tanks nested within water baths; compared with the tank-only model by a likelihood-ratio test (anova)",
  "Binomial mixed-effects model, aquarium level", "Effect of temperature on final mortality of VBM40-exposed fragments, accounting for aquarium and water bath",
  "glmer (lme4)", call_text(glmm_bath@call), paste0("Dead and alive fragments per aquarium and donor colony (", nrow(aquarium_counts), " counts, 8 aquaria)"),
  "Fixed effects: temperature (reference 25C) and donor colony (reference colony 1). Aquarium nested within water bath as random intercepts; compared with the aquarium-only model by a likelihood-ratio test",
  "Quasi-binomial block model, aquarium totals", "Temperature effect using one value per aquarium, water bath as a fixed block",
  "glm (stats)", call_text(glm_block$call), "Dead and alive fragments per aquarium (8 aquaria)",
  "Only 3 residual degrees of freedom; most conservative estimate",
  "Proportional hazards check", "Schoenfeld residual test for the Cox model", "cox.zph (survival)",
  paste0(call_text(cox_fixed$call), " followed by cox.zph()"), paste0("VBM40-exposed fragments (n = ", nrow(vbm40), ")"),
  "cox.zph does not support coxme, so the model without random effects is used",
  "Tank-level exact permutation tests", "Difference in mean tank mortality between groups of tanks",
  "base R (combn)", "difference in mean(prop_dead per tank); all assignments of group labels to tanks (combn)",
  "One value per tank (12 tanks)", "VBM40 vs Control within each temperature (15 allocations each), across temperatures stratified (225), temperature among VBM40 tanks (70)",
  "Water bath-restricted sign-flip tests", "Differences calculated within each water bath, signs flipped in all possible ways",
  "base R (expand.grid)", "mean(difference within water bath); all 2^k sign flips",
  "Tank-level values within water baths", "Temperature among VBM40 (4 baths); VBM40 vs Control (4 bath x temperature pairs)"))
Table_models
write.csv(Table_models, file = paste0(out, "Table_model_formulas.csv"), row.names = FALSE)

supp_tables <- list("Model_formulas" = Table_models,
                    "Kaplan_Meier" = Table_km,
                    "Log_rank" = Table_logrank,
                    "Fisher" = fisher_results,
                    "Cox_fixed_effects" = cox_results,
                    "Cox_random_and_fit" = cox_random,
                    "Cox_tank_vs_bath_LRT" = tibble::rownames_to_column(as.data.frame(cox_comparison), "Model"),
                    "PH_test" = tibble::rownames_to_column(as.data.frame(schoenfeld$table), "Term"),
                    "GLMM_binomial" = glmm_results,
                    "GLMM_random_effects" = glmm_random,
                    "GLMM_block_quasibinomial" = block_results,
                    "Waterbath_paired_t" = paired_t_results,
                    "Tank_permutation" = permutation_results,
                    "Waterbath_signflip" = bath_tests,
                    "Waterbath_pairs_temp" = bath_temp,
                    "Waterbath_pairs_treat" = bath_treat,
                    "Per_aquarium_by_colony" = mortality_tank_colony,
                    "Per_tank" = mortality_tank,
                    "Mortality_summary_CI" = mortality_summary)
wb <- createWorkbook()
walk2(supp_tables, names(supp_tables), function(df, nm) {
  addWorksheet(wb, nm)
  writeData(wb, nm, df, headerStyle = createStyle(textDecoration = "bold"))
  setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
})
saveWorkbook(wb, paste0(out, "Supplementary_tables_waterbath_v2.xlsx"), overwrite = TRUE)

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

ggsave(paste0(out, "Fig_KM_mortality_waterbath.png"), Fig_KM_mortality, width = 2800 / 300, height = 1400 / 300, dpi = 300)
## Vector PDF with editable text and objects
ggsave(paste0(out, "Fig_KM_mortality_waterbath.pdf"), Fig_KM_mortality, width = 2800 / 300, height = 1400 / 300, device = cairo_pdf)

### Session info ####
sessionInfo()
