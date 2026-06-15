# Demo.R
# Główny model bayesowski — pełna wersja do prezentacji.
# 4 łańcuchy × 2000 iteracji (1000 warmup) → 4000 próbek a posteriori.

library(brms)
library(dplyr)
library(ggplot2)
library(purrr)

# ============================================================
# 1. Dane
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

surv_final  |> 
group_by(compound)  |> 
summarise(median= median(time))  |> 
arrange(compound)
cat("n flowers:", nrow(surv_final), "\n")
cat("Censoring rate:", round(mean(surv_final$event == 0) * 100, 1), "%\n")

surv_final |>
  group_by(compound) |>
  summarise(
    n = n(),
    median_time = median(time, na.rm = TRUE),
    censored_pct = sum(can.be.presented == 1, na.rm = TRUE) / n() * 100,
    events = sum(can.be.presented == 1, na.rm = TRUE),
    .groups = "drop"
  )
# ============================================================
# 2. Priors (z Prior.R, dane pilotażowe)
# ============================================================
priors_main <- c(
  prior(gamma(15.466, 3.544), class = shape),
  prior(normal(2.918, 0.362), class = Intercept),
  prior(normal(0, 0.24),      class = b)
)

# ============================================================
# 3. Fit
# ============================================================
t_start <- Sys.time()

fit_main <- brm(
  formula = bf(
    time | cens(1 - event) ~ compound + species + garden,
    family = weibull()
  ),
  data    = surv_final,
  prior   = priors_main,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  file    = "fit_main_fixed"
)

t_end     <- Sys.time()
t_elapsed <- difftime(t_end, t_start, units = "mins")
timing_msg <- sprintf("fit_main_fixed — %s | czas: %.1f min\n",
                      format(t_end, "%Y-%m-%d %H:%M:%S"),
                      as.numeric(t_elapsed))
cat(timing_msg)
write(timing_msg, file = "fit_timings.txt", append = TRUE)

# ============================================================
# 4. Diagnostyka zbieżności
# ============================================================
print(summary(fit_main))

cat("\nMax R-hat:    ", round(max(rhat(fit_main), na.rm = TRUE), 4), "\n")
cat("Min Bulk-ESS:", round(min(neff_ratio(fit_main), na.rm = TRUE) * 4000, 0), "\n")
cat("N divergent:  ", sum(nuts_params(fit_main)$Value[nuts_params(fit_main)$Parameter == "divergent__"]), "\n")

fit_main <- readRDS("fit_main_fixed.rds")
# Trace plots
mcmc_plot(fit_main, type = "trace",
  variable = c("b_Intercept", "shape", "b_compound6"))

# Density overlay
mcmc_plot(fit_main, type = "dens_overlay",
  variable = c("b_Intercept", "shape", "b_compound6"))

# Autokorelacja
mcmc_plot(fit_main, type = "acf",
  variable = c("b_Intercept", "shape", "b_compound6"))

# Efekty związków α_c
mcmc_plot(fit_main, type = "intervals",
  variable = "^b_compound", regex = TRUE) +
  labs(title = "Efekty związków α_c", x = "Wartość", y = "Związek") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50")

# ============================================================
# 5. Posterior median survival per compound
#    λ_c = exp(β₀ + α_c) / Γ(1+1/k)
#    median_c = λ_c · (ln 2)^{1/k}
# ============================================================
draws         <- as_draws_df(fit_main)
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

cat("\n=== RANKING — MEDIANA PRZEŻYCIA (DNI) ===\n")
print(median_table, n = Inf)

# ============================================================
# 6. Prawdopodobieństwa przeżycia t = 10 i t = 22 dni
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

cat("\n=== PRAWDOPODOBIEŃSTWA PRZEŻYCIA t=10 i t=22 ===\n")
print(surv_probs, n = Inf)

# ============================================================
# 7. Wykres rankingowy
# ============================================================
compound_labels <- c(
  "1"  = "1 Distilled Water",         "4"  = "4 Concentrate of Caducues",
  "5"  = "5 Distillate of Discovery",  "6"  = "6 Essence of Epiphaneia",
  "7"  = "7 Four in December",         "8"  = "8 Granules of Geheref",
  "9"  = "9 Kar-Hamel Mooh",           "11" = "11 Noospherol",
  "12" = "12 Oil of John's Son",       "13" = "13 Power of Perlimpinpin",
  "14" = "14 Spirit of Scienza",       "15" = "15 Zest of Zen"
)

ctrl_median <- median_table$post_median[median_table$compound == 1]

rank_df <- median_table |>
  mutate(
    label      = compound_labels[as.character(compound)],
    label      = factor(label, levels = rev(compound_labels[as.character(compound)])),
    is_control = compound == 1
  )

p_rank <- ggplot(rank_df, aes(x = post_median, y = label, colour = is_control)) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.35) +
  geom_point(size = 3) +
  geom_vline(xintercept = ctrl_median, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  scale_colour_manual(values = c("FALSE" = "slateblue4", "TRUE" = "grey40"),
                      guide = "none") +
  labs(
    title    = "Ranking związków według mediany czasu przeżycia",
    subtitle = paste0("Mediana a posteriori z 95% CrI  |  linia przerywana = kontrola (",
                      round(ctrl_median, 1), " dni)"),
    x = "Mediana przeżycia (dni)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_rank

ggsave(p_rank,
  filename = "TeX/plots/compound_ranking.pdf",
  device = cairo_pdf, width = 8, height = 6
)

# ============================================================
# 8. Krzywe przeżycia a posteriori
# ============================================================
t_seq <- seq(0, 30, by = 0.5)

surv_curves <- map_dfr(compound_lvls, function(c_lvl) {
  k        <- draws$shape
  lambda_c <- get_lambda(c_lvl, draws)

  surv_mat <- matrix(0, nrow = length(k), ncol = length(t_seq))
  for (j in seq_along(t_seq)) {
    surv_mat[, j] <- exp(-(t_seq[j] / lambda_c)^k)
  }

  tibble(
    label = compound_labels[c_lvl],
    t     = t_seq,
    q025  = apply(surv_mat, 2, quantile, 0.025),
    q25   = apply(surv_mat, 2, quantile, 0.25),
    q50   = apply(surv_mat, 2, quantile, 0.50),
    q75   = apply(surv_mat, 2, quantile, 0.75),
    q975  = apply(surv_mat, 2, quantile, 0.975)
  )
})

p_surv <- ggplot(surv_curves, aes(x = t)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "slateblue2", alpha = 0.15) +
  geom_ribbon(aes(ymin = q25,  ymax = q75),  fill = "slateblue2", alpha = 0.30) +
  geom_line(aes(y = q50), colour = "slateblue4", linewidth = 0.8) +
  geom_vline(xintercept = c(10, 22), linetype = "dashed",
             colour = "tomato", linewidth = 0.5) +
  facet_wrap(~ label, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 30, by = 10)) +
  labs(
    title    = "Krzywe przeżycia a posteriori według związku",
    subtitle = "Mediana z przedziałami 50% i 95% CrI  |  linie przerywane: dni 10 i 22",
    x = "Czas (dni)", y = "P(przeżycia)"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text    = element_text(size = 7.5),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"))

p_surv

ggsave(p_surv,
  filename = "TeX/plots/posterior_surv_curves.pdf",
  device = cairo_pdf, width = 10, height = 12
)


fit_main <- readRDS("fit_main_fixed.rds")
summary(fit)
pp_check(fit)
