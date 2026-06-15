# Prior.R
# Requires surv_rose from EDA.R in the current R session (run EDA.R first).
# Libraries: dplyr, ggplot2, survival, tidyverse, patchwork.

# ============================================================
# 1. Per-cell Weibull fits (intercept-only)
#    One model per compound x species x garden cell.
#    Used for: beta_0, alpha_c, gamma_s, sigma_g priors.
# ============================================================
pilot_estimates <- surv_rose |>
  mutate(time = if_else(time == 0, 0.5, as.numeric(time))) |>
  group_by(compound, species, garden) |>
  group_modify(~ {
    fit <- survreg(Surv(time, can.be.presented) ~ 1, data = .x, dist = "weibull")
    tibble(
      k_hat      = 1 / fit$scale,
      mu_hat     = exp(coef(fit)[["(Intercept)"]]),
      log_mu_hat = coef(fit)[["(Intercept)"]]
    )
  }) |>
  ungroup() |>
  mutate(
    compound = factor(compound),
    species  = factor(species),
    garden   = factor(garden)
  )

# ============================================================
# 2. Shape k prior — global model + delta method
#    Single shared-k Weibull model (consistent with model assumption).
#    SE inflated by 2x for conservatism. Moment-matched to Gamma.
# ============================================================
m_adj <- survreg(
  Surv(time, can.be.presented) ~ factor(compound),
  data = filter(surv_rose, time > 0),
  dist = "weibull"
)

se_scale <- summary(m_adj)$table["Log(scale)", "Std. Error"]
k_hat    <- 1 / m_adj$scale
se_k     <- se_scale / m_adj$scale^2
k_var    <- (2 * se_k)^2
k_alpha  <- k_hat^2 / k_var
k_beta   <- k_hat / k_var

cat("k̂ =", round(k_hat, 3), "  SE(k) =", round(se_k, 3), "\n")
cat("k ~ Gamma(", round(k_alpha, 3), ",", round(k_beta, 3), ")\n")
cat(
  "E[k] =", round(k_alpha / k_beta, 3),
  "  SD[k] =", round(sqrt(k_alpha) / k_beta, 3), "\n"
)

# ============================================================
# 3. Intercept beta_0 prior — Normal on log(lambda) scale
#    Mean and SD from per-cell log(mu) estimates; SD inflated 1.5x.
# ============================================================
logmu_mean         <- mean(pilot_estimates$log_mu_hat)
logmu_sd           <- sd(pilot_estimates$log_mu_hat)
logmu_sd_inflated  <- logmu_sd * 1.5

# ============================================================
# 4. Fixed effects prior: compound alpha_c and species gamma_s
#    SD of compound x species cell means (averaged over gardens)
#    so garden noise does not inflate the prior.
#    Same prior applied to both alpha_c and gamma_s.
# ============================================================
cell_means <- pilot_estimates |>
  group_by(compound, species) |>
  summarise(mean_log_mu = mean(log_mu_hat), .groups = "drop")

b_prior_sd <- round(sd(cell_means$mean_log_mu), 2)

# ============================================================
# 5. Garden SD prior sigma_g — data-informed from pilot
#    SD of log(mu) differences between gardens per compound x species pair.
#    Inflated by 1.5x for conservatism.
# ============================================================
garden_diffs <- pilot_estimates |>
  select(compound, species, garden, log_mu_hat) |>
  pivot_wider(
    names_from  = garden,
    values_from = log_mu_hat,
    names_prefix = "g"
  ) |>
  mutate(diff = g2 - g1)

garden_sd_observed <- sd(garden_diffs$diff)
garden_sd_prior    <- round(garden_sd_observed * 1.5, 3)

# ============================================================
# 6. Prior predictive check (numerical summary)
# ============================================================
set.seed(42)
n_sim  <- 10000
k_sim  <- rgamma(n_sim, shape = k_alpha, rate = k_beta)
mu_sim <- exp(rnorm(n_sim, mean = logmu_mean, sd = logmu_sd_inflated))
t_sim  <- rweibull(n_sim, shape = k_sim, scale = mu_sim)

cat("=== Prior Predictive Check ===\n")
cat("Median simulated survival:", round(median(t_sim), 1), "days\n")
cat("% surviving past day 30:  ", round(mean(t_sim > 30) * 100, 1), "%\n")
cat("% dying within day 5:     ", round(mean(t_sim < 5) * 100, 1), "%\n")
cat(
  "90% interval: [",
  round(quantile(t_sim, 0.05), 1), ",",
  round(quantile(t_sim, 0.95), 1), "] days\n"
)

# ============================================================
# 7. Final prior summary
# ============================================================
cat("\n=== PRIORS FOR SAP ===\n")
cat("Weibull shape k  ~ Gamma(", round(k_alpha, 3), ",", round(k_beta, 3), ")\n")
cat("Intercept        ~ Normal(", round(logmu_mean, 3), ",", round(logmu_sd_inflated, 3), ")\n")
cat("Compound effects ~ Normal(0,", b_prior_sd, ")\n")
cat("Species effect   ~ Normal(0,", b_prior_sd, ")\n")
cat("Garden SD        ~ HalfNormal(0,", garden_sd_prior, ")\n")

# ══════════════════════════════════════════════════════════
#  PLOT 1 — log(mu) estimates with prior overlay
# ══════════════════════════════════════════════════════════
n_compounds     <- nlevels(pilot_estimates$compound)
compound_levels <- sort(unique(as.integer(as.character(pilot_estimates$compound))))

prior_ribbon <- tibble(
  x  = seq(0.5, n_compounds + 0.5, length.out = 200),
  lo = qnorm(0.05, logmu_mean, logmu_sd_inflated),
  hi = qnorm(0.95, logmu_mean, logmu_sd_inflated),
  mu = logmu_mean
)

p1 <- ggplot(
  pilot_estimates,
  aes(x = as.numeric(compound), y = log_mu_hat)
) +
  geom_ribbon(
    data = prior_ribbon,
    aes(x = x, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    fill = "#378ADD", alpha = 0.12
  ) +
  geom_hline(
    data = prior_ribbon,
    aes(yintercept = mu),
    colour = "#378ADD", linewidth = 0.6,
    linetype = "solid", inherit.aes = FALSE
  ) +
  geom_point(aes(shape = species), size = 2.5, alpha = 0.85) +
  scale_shape_manual(
    values = c("1" = 16, "2" = 17),
    labels = c("1" = "Species 1", "2" = "Species 2"),
    name   = NULL
  ) +
  scale_x_continuous(
    breaks = seq_along(compound_levels),
    labels = compound_levels
  ) +
  labs(
    title    = "Plot 1 — scale estimates",
    subtitle = "log(μ) per cell with 90% prior interval (blue)",
    x = "Compound", y = "log(μ̂)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ══════════════════════════════════════════════════════════
#  PLOT 2 — garden effect (log mu diff per compound x species)
# ══════════════════════════════════════════════════════════
garden_diffs_plot <- garden_diffs |>
  mutate(
    label    = paste0("C", compound, " S", species),
    positive = diff >= 0
  )

p2 <- ggplot(
  garden_diffs_plot,
  aes(x = reorder(label, diff), y = diff, fill = positive)
) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_hline(
    yintercept =  sd(garden_diffs$diff),
    linetype = "dashed", colour = "#854F0B", linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -sd(garden_diffs$diff),
    linetype = "dashed", colour = "#854F0B", linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#5DCAA5", "FALSE" = "#F0997B"),
    labels = c("TRUE" = "Garden 2 > Garden 1", "FALSE" = "Garden 1 > Garden 2"),
    name   = NULL
  ) +
  coord_flip() +
  labs(
    title    = "Plot 2 — garden effect",
    subtitle = "log(μ) difference garden 2 − garden 1 | dashed = ±1 SD",
    x = NULL, y = "Δ log(μ)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ══════════════════════════════════════════════════════════
#  PLOT 3 — prior predictive distribution (histogram)
# ══════════════════════════════════════════════════════════
pred_df <- tibble(t = pmin(t_sim, 60))

p3 <- ggplot(pred_df, aes(x = t)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40, fill = "#7F77DD", colour = "white", linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 30, colour = "#E24B4A",
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_vline(
    xintercept = median(t_sim), colour = "#1D9E75",
    linetype = "dotted", linewidth = 0.8
  ) +
  annotate("text",
    x = 31, y = Inf, label = "Day 30\n(study end)",
    hjust = 0, vjust = 1.4, size = 3, colour = "#E24B4A"
  ) +
  annotate("text",
    x = median(t_sim) - 1, y = Inf,
    label = paste0("Median\n", round(median(t_sim), 1), "d"),
    hjust = 1, vjust = 1.4, size = 3, colour = "#1D9E75"
  ) +
  labs(
    title    = "Plot 3 — prior predictive check",
    subtitle = paste0(
      round(mean(t_sim > 30) * 100, 1), "% survive past day 30  |  ",
      round(mean(t_sim < 5)  * 100, 1), "% die within day 5"
    ),
    x = "Simulated survival time (days, capped at 60)", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ══════════════════════════════════════════════════════════
#  PLOT 4 — k prior density  →  TeX/plots/prior_k_density.pdf
# ══════════════════════════════════════════════════════════
k_df <- data.frame(k = seq(0.01, 10, length.out = 500)) |>
  mutate(density = dgamma(k, shape = k_alpha, rate = k_beta))

prior_density_plot <- ggplot(k_df, aes(x = k, y = density)) +
  geom_area(fill = "slateblue2", alpha = 0.15) +
  geom_line(color = "slateblue2", linewidth = 1.2) +
  geom_vline(xintercept = k_hat, color = "tomato", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 1,     color = "grey40", linetype = "dotted", linewidth = 1) +
  annotate("label",
    x = k_hat, y = max(k_df$density) * 0.7,
    label = paste0("hat(k)[MLE] == ", round(k_hat, 2)), parse = TRUE,
    color = "tomato", size = 3.5, hjust = -0.05
  ) +
  annotate("label",
    x = 1, y = max(k_df$density) * 0.1,
    label = "k == 1~(stały~hazard)", parse = TRUE,
    color = "grey40", size = 3, hjust = -0.05
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(
    title    = "Gęstość prawdopodobieństwa",
    subtitle = paste0("k ~ Gamma(", round(k_alpha, 3), ", ", round(k_beta, 3), ")"),
    x = "k", y = ""
  ) +
  theme_minimal()
prior_density_plot

# ══════════════════════════════════════════════════════════
#  PLOT 5 — prior predictive survival curves  →  TeX/plots/prior_pred_surv.pdf
# ══════════════════════════════════════════════════════════
set.seed(42)
n_draws      <- 300
k_draws      <- rgamma(n_draws, shape = k_alpha, rate = k_beta)
lambda_draws <- exp(rnorm(n_draws, mean = logmu_mean, sd = logmu_sd_inflated))
t_seq        <- seq(0, 30, by = 0.5)

surv_prior <- map_dfr(1:n_draws, function(i) {
  tibble(
    t    = t_seq,
    surv = exp(-(t_seq / lambda_draws[i])^k_draws[i]),
    draw = i
  )
})

surv_bands <- surv_prior |>
  group_by(t) |>
  summarise(
    q025 = quantile(surv, 0.025),
    q25  = quantile(surv, 0.25),
    q50  = quantile(surv, 0.50),
    q75  = quantile(surv, 0.75),
    q975 = quantile(surv, 0.975),
    .groups = "drop"
  )

prior_pred_plot <- ggplot(surv_bands, aes(x = t)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "slateblue2", alpha = 0.20) +
  geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = "slateblue2", alpha = 0.35) +
  geom_line(aes(y = q50), color = "slateblue4", linewidth = 1.2) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title    = "Prior predictive check krzywych przeżycia",
    subtitle = sprintf(
      "Mediana oraz przedziały 50%% i 95%% (n = 300)\nPrior: k ~ Gamma(%.3f, %.3f), log(λ) ~ N(%.2f, %.2f)",
      k_alpha, k_beta, logmu_mean, logmu_sd_inflated
    ),
    x = "Czas (dni)",
    y = "Prawdopodobieństwo przeżycia"
  ) +
  theme_minimal()
prior_pred_plot

# ══════════════════════════════════════════════════════════
#  SAVE plots 4 and 5
# ══════════════════════════════════════════════════════════
ggsave(prior_density_plot,
  filename = "TeX/plots/prior_k_density.pdf",
  device = cairo_pdf, width = 6, height = 4
)
ggsave(prior_pred_plot,
  filename = "TeX/plots/prior_pred_surv.pdf",
  device = cairo_pdf, width = 7, height = 5
)

# ══════════════════════════════════════════════════════════
#  COMBINE plots 1-3 with patchwork
# ══════════════════════════════════════════════════════════
(p1 + p2) / p3 +
  plot_annotation(
    title   = "Prior elicitation — pilot data summary",
    caption = paste0(
      "Shape k ~ Gamma(", round(k_alpha, 3), ", ", round(k_beta, 3), ")  |  ",
      "Intercept ~ Normal(", round(logmu_mean, 3), ", ", round(logmu_sd_inflated, 3), ")  |  ",
      "b ~ Normal(0, ", b_prior_sd, ")  |  ",
      "Garden SD ~ HalfNormal(0, ", garden_sd_prior, ")"
    ),
    theme = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 9, colour = "grey50")
    )
  )
