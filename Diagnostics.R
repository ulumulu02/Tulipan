# Diagnostics.R
# Generates diagnostic and validation plots for the final report.
# Outputs: ppc_km.pdf, trace_plots.pdf, acf_plots.pdf, rank_plots.pdf
# Run from project root after fit_main_fixed.rds is available.

library(brms)
library(dplyr)
library(ggplot2)
library(purrr)
library(survival)
library(bayesplot)


bayesplot::color_scheme_set("mix-pink-purple")
bayesplot::color_scheme_view()
options(brms.backend = "cmdstanr")


# ── Setup ──────────────────────────────────────────────────────────────────────
fit <- readRDS("models/fit_main_fixed.rds")

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

compound_labels <- c(
  "1"  = "1 Distilled Water",
  "4"  = "4 Concentrate of Caducues",
  "5"  = "5 Distillate of Discovery",
  "6"  = "6 Essence of Epiphaneia",
  "7"  = "7 Four in December",
  "8"  = "8 Granules of Geheref",
  "9"  = "9 Kar-Hamel Mooh",
  "11" = "11 Noospherol",
  "12" = "12 Oil of John's Son",
  "13" = "13 Power of Perlimpinpin",
  "14" = "14 Spirit of Scienza",
  "15" = "15 Zest of Zen"
)
compound_lvls <- names(compound_labels)

draws <- as_draws_df(fit)

get_lambda <- function(c_lvl, draws) {
  k   <- draws$shape
  b0  <- draws$b_Intercept
  a_c <- if (c_lvl == "1") 0 else draws[[paste0("b_compound", c_lvl)]]
  exp(b0 + a_c) / gamma(1 + 1/k)
}

# ── 1. Posterior Predictive Check: KM overlay ─────────────────────────────────
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

km_df <- map_dfr(compound_lvls, function(c_lvl) {
  sub <- surv_final[surv_final$compound == c_lvl, ]
  km  <- survfit(Surv(time, event) ~ 1, data = sub)
  tibble(
    compound = as.integer(c_lvl),
    label    = compound_labels[c_lvl],
    t        = c(0, km$time),
    surv     = c(1, km$surv)
  )
})

# Fix factor levels so facets appear in ranking order
rank_order <- c("6 Essence of Epiphaneia", "5 Distillate of Discovery",
                "11 Noospherol", "8 Granules of Geheref",
                "12 Oil of John's Son", "13 Power of Perlimpinpin",
                "4 Concentrate of Caducues", "14 Spirit of Scienza",
                "15 Zest of Zen", "7 Four in December",
                "9 Kar-Hamel Mooh", "1 Distilled Water")

surv_curves <- surv_curves |>
  mutate(label = factor(label, levels = rank_order))
km_df <- km_df |>
  mutate(label = factor(label, levels = rank_order))

p_ppc <- ggplot() +
  geom_ribbon(data = surv_curves,
              aes(x = t, ymin = q025, ymax = q975),
              fill = "#C868C8", alpha = 0.15) +
  geom_ribbon(data = surv_curves,
              aes(x = t, ymin = q25, ymax = q75),
              fill = "#C868C8", alpha = 0.35) +
  geom_line(data = surv_curves,
            aes(x = t, y = q50),
            colour = "#7B2D87", linewidth = 0.8) +
  geom_step(data = km_df,
            aes(x = t, y = surv),
            colour = "black", linewidth = 0.55) +
  geom_vline(xintercept = c(10, 22), linetype = "dashed",
             colour = "tomato", linewidth = 0.4) +
  facet_wrap(~ label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 30, by = 10)) +
  labs(
    title    = "Posterior predictive check: krzywe przeżycia a posteriori vs Kaplan-Meier",
    subtitle = "Pasma: 50% i 95% CrI  |  czarna linia schodkowa: estymator KM  |  przerywane: t = 10 i 22 dni",
    x        = "Czas (dni)",
    y        = "P(przeżycie)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text    = element_text(size = 7.5),
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, colour = "grey40")
  )

ggsave(p_ppc,
  filename = "TeX/plots/ppc_km.pdf",
  device = cairo_pdf, width = 10, height = 12)
cat("Saved: ppc_km.pdf\n")

# ── 2–4. MCMC diagnostics (trace, ACF, rank) ──────────────────────────────────
post_array <- as.array(fit)

key_pars <- c("shape", "b_Intercept", "b_compound6", "b_compound9")

par_labels <- c(
  shape        = "k  (parametr kształtu)",
  b_Intercept  = "β₀  (wyraz wolny)",
  b_compound6  = "α₆  (Essence of Epiphaneia)",
  b_compound9  = "α₉  (Kar-Hamel Mooh)"
)

# Trace plots
p_trace <- mcmc_trace(post_array, pars = key_pars,
                      facet_args = list(labeller = as_labeller(par_labels))) +
  scale_x_continuous(breaks = c(0, 500, 1000)) +
  labs(title = "Trace plots — wybrane parametry modelu") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(p_trace,
  filename = "TeX/plots/trace_plots.pdf",
  device = cairo_pdf, width = 10, height = 8)
cat("Saved: trace_plots.pdf\n")

# ACF plots
p_acf <- mcmc_acf_bar(post_array, pars = key_pars, lags = 20,
                      facet_args = list(labeller = as_labeller(par_labels))) +
  labs(title = "Autokorelacja — wybrane parametry modelu") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(p_acf,
  filename = "TeX/plots/acf_plots.pdf",
  device = cairo_pdf, width = 10, height = 6)
cat("Saved: acf_plots.pdf\n")

# Rank overlay plots
p_rank <- mcmc_rank_overlay(post_array, pars = key_pars,
                            facet_args = list(labeller = as_labeller(par_labels))) +
  labs(title = "Wykresy rang — wybrane parametry modelu") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

ggsave(p_rank,
  filename = "TeX/plots/rank_plots.pdf",
  device = cairo_pdf, width = 10, height = 6)
cat("Saved: rank_plots.pdf\n")

cat("\nAll diagnostic plots saved to TeX/plots/\n")
