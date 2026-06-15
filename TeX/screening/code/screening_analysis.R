library(survival)
library(ggplot2)
library(survminer) 
library(dplyr)
# Survival Curves - plot ---------------------
control_fit <- survfit(Surv(time, can.be.presented) ~ 1,
  data = filter(surv_rose, compound == 1)
)

control_df <- data.frame(
  time   = control_fit$time,
  surv   = control_fit$surv,
  upper  = control_fit$upper,
  lower  = control_fit$lower
)

# KM for all compounds
km_all <- survfit(Surv(time, can.be.presented) ~ compound, data = surv_rose)
# Podział na grupy - 0 (woda), 1 szkodliwe, 2 umiarkowane, 3 korzystne
print(km_all)

# Extract KM data for all compounds using survminer
km_df <- surv_summary(km_all, data = surv_rose) %>%
  mutate(
    compound = as.integer(gsub("compound=", "", strata)),
    compound_label = compound_names[compound]
  )

# Add control label column for facet titles
km_df$compound_label <- factor(km_df$compound_label,
  levels = compound_names
)

# Plot
ggplot() +
  # Grey control ribbon + line on every facet
  geom_ribbon(
    data = control_df,
    aes(x = time, ymin = lower, ymax = upper),
    fill = "grey70", alpha = 0.3, na.rm = TRUE
  ) +
  geom_step(
    data = control_df,
    aes(x = time, y = surv),
    color = "grey50", linewidth = 0.8, linetype = "dashed"
  ) +
  # Colored compound curves
  geom_ribbon(
    data = km_df,
    aes(x = time, ymin = lower, ymax = upper, fill = compound_label),
    alpha = 0.25, na.rm = TRUE
  ) +
  geom_step(
    data = km_df,
    aes(x = time, y = surv, color = compound_label),
    linewidth = 0.8
  ) +
  facet_wrap(~compound_label, ncol = 3) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(
    title = "Krzywe przeżycia Kaplana-Meiera według związku",
    subtitle = "Szara linia przerywana = Związek 1 (Woda destylowana) -- grupa kontrolna",
    x = "Czas (dni)",
    y = "Prawdopodobieństwo przeżycia"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8),
    plot.subtitle = element_text(color = "grey50")
  )

ggsave("TeX/pre_analysis/plots/km_compounds_facet.pdf",
  device = cairo_pdf, width = 10, height = 14, units = "in"
)


# One-sided log-rank test: each compound vs control (compound 1) ---------------
# H0: S_c(t) = S_1(t)
# H1: S_c(t) > S_1(t)  <=>  compound c has shorter survival than control

compounds <- setdiff(sort(unique(surv_rose$compound)), 1)
med_ctrl  <- median(surv_rose$time[surv_rose$compound == 1])

results <- lapply(compounds, function(c) {
  df  <- surv_rose %>% filter(compound %in% c(1, c))
  fit <- survdiff(Surv(time, can.be.presented) ~ compound, data = df)

  # Direction: compound c is worse if it has more observed events than expected
  idx_c       <- which(names(fit$n) == paste0("compound=", c))
  worse       <- fit$obs[idx_c] > fit$exp[idx_c]
  p_two_sided <- pchisq(fit$chisq, df = 1, lower.tail = FALSE)
  p_one_sided <- ifelse(worse, p_two_sided / 2, 1 - p_two_sided / 2)

  data.frame(
    compound    = c,
    median      = median(surv_rose$time[surv_rose$compound == c]),
    p_one_sided = p_one_sided
  )
})
# Poprawka Holma (FWER)
results_df <- bind_rows(results) %>%
  mutate(
    p_holm  = p.adjust(p_one_sided, method = "holm"),
    crit1   = p_holm < 0.05,            # statistically significant after correction
    exclude = crit1 
  )

print(results_df)
cat("\nKontrola - mediana:", med_ctrl, "dni\n")
cat("Próg praktyczny (50% mediany kontroli):", 0.5 * med_ctrl, "dni\n\n")
cat("Związki wykluczone z badania głównego:\n")
print(results_df %>% filter(exclude) %>% select(compound, median, p_holm))
