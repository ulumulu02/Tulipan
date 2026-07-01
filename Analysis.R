# Analysis.R
# Bayesian Weibull survival analysis — main study (roses_analysis2_final.csv).
# Run EDA.R + Prior.R first (surv_rose and prior values must be available).
# Libraries: brms, dplyr, tidyr, ggplot2, purrr.

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
# cmdstanr backend setup for faster and more efficient computing
options(brms.backend = "cmdstanr")
# ============================================================
# 1. Data preparation (mirrors EDA.R logic for pilot data)
# ============================================================
rose_final <- read.csv2("roses_analysis2_final.csv")

surv_final <- rose_final |>
  group_by(id) |>
  slice_max(order_by = time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    event    = ifelse(can.be.presented == FALSE, 1L, 0L),
    time     = pmax(as.numeric(time), 0.5),
    compound = factor(compound),
    species  = factor(species),
    garden   = factor(garden)
  ) |>
  filter(!(compound %in% c(2, 3, 10))) |>
  droplevels()
rm(rose_final)
cat("n flowers:", nrow(surv_final), "\n")
cat("Compounds:", paste(levels(surv_final$compound), collapse = ", "), "\n")
cat("Censoring rate:", round(mean(surv_final$event == 0) * 100, 1), "%\n")

# ============================================================
# 2. Priors calibrated on pilot data (from Prior.R)
# ============================================================
priors_main <- c(
  prior(gamma(15.466, 3.544), class = shape),
  prior(normal(2.918, 0.362), class = Intercept),
  prior(normal(0, 0.24),      class = b)
)

# Sensitivity analysis: weakly informative priors (SAP Section 6)
priors_sensitivity <- c(
  prior(gamma(2, 0.5), class = shape),
  prior(normal(0, 1),  class = Intercept),
  prior(normal(0, 1),  class = b)
)

# ============================================================
# 3. Model formula
#    cens(1 - event): event=1 (died) → cens=0 (observed)
#                     event=0 (alive) → cens=1 (right-censored)
#    brms Weibull: log link on mean mu; scale λ = mu / Γ(1+1/k)
# ============================================================
weibull_formula <- bf(
  time | cens(1 - event) ~ compound + species + garden,
  family = weibull()
)

# ============================================================
# 4. Main model fit
# ============================================================
fit_main <- brm(
  formula = weibull_formula,
  data    = surv_final,
  prior   = priors_main,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 42,
  file    = "models/fit_main"
)

# ============================================================
# 5. Sensitivity model fit
# ============================================================
fit_sens <- brm(
  formula = weibull_formula,
  data    = surv_final,
  prior   = priors_sensitivity,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 42,
  file    = "fit_sensitivity"
)

# ============================================================
# 6. Convergence diagnostics
# ============================================================
cat("\n=== MAIN MODEL SUMMARY ===\n")
print(summary(fit_main))

cat("\nMax R-hat:", round(max(rhat(fit_main), na.rm = TRUE), 4), "\n")
cat("Min Bulk-ESS:", round(min(neff_ratio(fit_main), na.rm = TRUE) * 4000, 0), "\n")

mcmc_plot(fit_main, type = "trace",
  variable = c("b_Intercept", "shape", "b_compound6"))

# ============================================================
# 7. Posterior median survival per compound
#    λ_c = exp(β₀ + α_c) / Γ(1+1/k)  [scale param from brms mu]
#    median_c = λ_c · (ln 2)^{1/k}
#    (reference: species = 1, garden RE = 0)
# ============================================================
draws <- as_draws_df(fit_main)
compound_lvls <- levels(surv_final$compound)

get_lambda <- function(c_lvl, draws) {
  k   <- draws$shape
  b0  <- draws$b_Intercept
  a_c <- if (c_lvl == "1") 0 else draws[[paste0("b_compound", c_lvl)]]
  exp(b0 + a_c) / gamma(1 + 1/k)
}

median_table <- map_dfr(compound_lvls, function(c_lvl) {
  lambda_c <- get_lambda(c_lvl, draws)
  median_c <- lambda_c * log(2)^(1/draws$shape)
  tibble(
    compound    = as.integer(c_lvl),
    post_median = round(median(median_c), 2),
    q025        = round(quantile(median_c, 0.025), 2),
    q975        = round(quantile(median_c, 0.975), 2)
  )
}) |>
  arrange(desc(post_median))

cat("\n=== COMPOUND RANKING — MEDIAN SURVIVAL (DAYS) ===\n")
print(median_table, n = Inf)

# ============================================================
# 8. Survival probabilities at t = 10 and t = 22 days
# ============================================================
surv_probs <- map_dfr(compound_lvls, function(c_lvl) {
  k        <- draws$shape
  lambda_c <- get_lambda(c_lvl, draws)
  tibble(
    compound  = as.integer(c_lvl),
    P_surv_10 = round(mean(exp(-(10 / lambda_c)^k)), 3),
    P_surv_22 = round(mean(exp(-(22 / lambda_c)^k)), 3)
  )
}) |>
  arrange(desc(P_surv_22))

cat("\n=== SURVIVAL PROBABILITIES AT t=10 AND t=22 ===\n")
print(surv_probs, n = Inf)

# ============================================================
# 9. Posterior survival curves
# ============================================================
compound_labels <- c(
  "1" = "1 Distilled Water",        "4"  = "4 Concentrate of Caducues",
  "5" = "5 Distillate of Discovery", "6"  = "6 Essence of Epiphaneia",
  "7" = "7 Four in December",        "8"  = "8 Granules of Geheref",
  "9" = "9 Kar-Hamel Mooh",          "11" = "11 Noospherol",
  "12" = "12 Oil of John's Son",     "13" = "13 Power of Perlimpinpin",
  "14" = "14 Spirit of Scienza",     "15" = "15 Zest of Zen"
)

t_seq <- seq(0, 30, by = 0.5)

surv_curves <- map_dfr(compound_lvls, function(c_lvl) {
  k        <- draws$shape
  lambda_c <- get_lambda(c_lvl, draws)

  n_t      <- length(t_seq)
  surv_mat <- matrix(0, nrow = length(k), ncol = n_t)
  for (j in seq_len(n_t)) {
    surv_mat[, j] <- exp(-(t_seq[j] / lambda_c)^k)
  }

  tibble(
    compound = as.integer(c_lvl),
    label    = compound_labels[c_lvl],
    t        = t_seq,
    q025     = apply(surv_mat, 2, quantile, 0.025),
    q25      = apply(surv_mat, 2, quantile, 0.25),
    q50      = apply(surv_mat, 2, quantile, 0.50),
    q75      = apply(surv_mat, 2, quantile, 0.75),
    q975     = apply(surv_mat, 2, quantile, 0.975)
  )
})

p_surv <- ggplot(surv_curves, aes(x = t)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#C868C8", alpha = 0.15) +
  geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = "#C868C8", alpha = 0.35) +
  geom_line(aes(y = q50), colour = "#7B2D87", linewidth = 0.8) +
  geom_vline(xintercept = c(10, 22), linetype = "dashed",
             colour = "tomato", linewidth = 0.5) +
  facet_wrap(~ label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 30, by = 10)) +
  labs(
    title    = "Posterior predictive survival curves by compound",
    subtitle = "Median with 50% and 95% credible bands | dashed: days 10 and 22",
    x        = "Time (days)", y = "P(survive)"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 7.5),
        plot.title = element_text(face = "bold"))

p_surv

ggsave(p_surv,
  filename = "TeX/plots/posterior_surv_curves.pdf",
  device = cairo_pdf, width = 10, height = 12
)

# ============================================================
# 10. Compound ranking plot
# ============================================================
rank_df <- median_table |>
  mutate(
    compound = factor(compound, levels = rev(median_table$compound)),
    is_control = compound == "1"
  )

ctrl_median <- median_table$post_median[median_table$compound == 1]

p_rank <- ggplot(rank_df, aes(x = post_median, y = compound, colour = is_control)) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.35) +
  geom_point(size = 3) +
  geom_vline(xintercept = ctrl_median, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  scale_colour_manual(values = c("FALSE" = "#7B2D87", "TRUE" = "grey40"),
                      guide = "none") +
  labs(
    title    = "Ranking związków według mediany czasu przeżycia",
    subtitle = paste0("Mediana a posteriori z 95% CrI  |  linia przerywana = kontrola (",
                      round(ctrl_median, 1), " dni)"),
    x = "Mediana przeżycia (dni)", y = "Związek"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_rank

ggsave(p_rank,
  filename = "TeX/plots/compound_ranking.pdf",
  device = cairo_pdf, width = 7, height = 6
)

# ============================================================
# 11. Sensitivity analysis comparison (run after both models fit)
# ============================================================
draws_sens <- as_draws_df(fit_sens)

compare_sensitivity <- map_dfr(compound_lvls, function(c_lvl) {
  k_m   <- draws$shape
  k_s   <- draws_sens$shape
  lam_m <- get_lambda(c_lvl, draws)
  lam_s <- exp(
    draws_sens$b_Intercept +
      (if (c_lvl == "1") 0 else draws_sens[[paste0("b_compound", c_lvl)]])
  ) / gamma(1 + 1/k_s)

  tibble(
    compound      = as.integer(c_lvl),
    median_main   = round(median(lam_m * log(2)^(1/k_m)), 2),
    median_sens   = round(median(lam_s * log(2)^(1/k_s)), 2)
  )
}) |>
  arrange(desc(median_main))

cat("\n=== SENSITIVITY CHECK — RANKING COMPARISON ===\n")
print(compare_sensitivity, n = Inf)
