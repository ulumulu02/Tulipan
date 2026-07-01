
rose <- read.csv2("roses_analysis2.csv")
compound_names <- c(
  `1`  = "1 Distilled Water",
  `4`  = "4 Concentrate of Caducues",
  `5`  = "5 Distillate of Discovery",
  `6`  = "6 Essence of Epiphaneia",
  `7`  = "7 Four in December",
  `8`  = "8 Granules of Geheref",
  `9`  = "9 Kar-Hamel Mooh",
  `11` = "11 Noospherol",
  `12` = "12 Oil of John's Son",
  `13` = "13 Power of Perlimpinpin",
  `14` = "14 Spirit of Scienza",
  `15` = "15 Zest of Zen"
)

# Data prep -----------------
surv_rose <- as.data.frame(rose) %>%
  left_join(
    rose %>%
      filter(time == 0) %>%
      select(id, baseline_diam = diam),
    by = "id"
  ) %>%
  group_by(id) %>%
  slice_max(order_by = time, n = 1, with_ties = FALSE) %>%
  mutate(can.be.presented = ifelse(can.be.presented == FALSE, 1, 0))  |> 
  filter(!(compound %in% c(2,3,10))) |>
  ungroup()
rm(rose)

head(surv_rose)
# Exploratory Data Analysis (EDA) ----------------------------------------------

## 1 - Summary table of compounds ------------------------------------------

sum_table_time <- surv_rose %>%
  group_by(compound) %>%
  summarise(
    no_flowers = n(),
    mean_surv_time = mean(time),
    med_surv_time = median(time),
    sd_surv_time = sd(time),
    q1 = quantile(time, 0.25, na.rm = TRUE),
    q3 = quantile(time, 0.75, na.rm = TRUE),
    min = min(time),
    max=max(time)
  )
sum_table_time
latex_table <- xtable(sum_table_time,
  align = c("c", "c", "c", "c", "c", "c","c","c","c","c")
)
print(
  latex_table,
  type = "latex",
  include.rownames = FALSE
)

## 2 - Boxplot of survival time by compound ------------------------------
compound_boxplot <- surv_rose %>%
  ggplot(aes(x = factor(compound), y = time)) +
  geom_boxplot(alpha = 0.3, fill = "#C868C8", color= "#7B2D87") +
  coord_flip() +
  scale_x_discrete(labels = compound_names) +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  labs(
    title = "Rozkład czasu przeżycia kwiatów według zastosowanego związku",
    x = NULL,
    y = "Czas przeżycia (dni)"
  ) +
  theme_minimal() +
theme(
  plot.title = element_markdown(
    hjust = 0,
    margin = margin(t = 10, r = 10, b = 10, l = 10)
  ),
  plot.title.position = "plot",
  legend.position = "none"
)

compound_boxplot

ggsave(compound_boxplot,
  filename = "TeX/SAP/plots/boxplot_compounds.pdf",
  device = cairo_pdf,
  width = 9,
  height = 6
)

## 3 - KM Survival curves -----------------------------------------------------
km_all <- survfit(Surv(time, can.be.presented) ~ compound, data = surv_rose)

km_df <- surv_summary(km_all, data = surv_rose) %>%
  mutate(
    compound_raw = as.character(gsub("compound=", "", as.character(strata))),
    compound_label = factor(
      compound_names[compound_raw],
      levels = compound_names
    )
  )

# --- kontrola (compound == 1) ---
control_fit <- survfit(
  Surv(time, can.be.presented) ~ 1,
  data = filter(surv_rose, compound == 1)
)
unique(gsub("compound=", "", km_df$strata))
control_df <- data.frame(
  time  = control_fit$time,
  surv  = control_fit$surv,
  lower = control_fit$lower,
  upper = control_fit$upper
)

# --- powiel kontrolę na wszystkie panele ---
control_df_all <- km_df %>%
  distinct(compound_label) %>%
  tidyr::crossing(control_df)

# --- wykres ---
km_plot <- ggplot() +
  # kontrola (we wszystkich panelach)
  geom_ribbon(
    data = control_df_all,
    aes(x = time, ymin = lower, ymax = upper),
    fill = "grey70", alpha = 0.25
  ) +
  geom_step(
    data = control_df_all,
    aes(x = time, y = surv),
    color = "grey40", linetype = "dashed", linewidth = 0.7
  ) +

  # krzywe compoundów
  geom_ribbon(
    data = km_df,
    aes(x = time, ymin = lower, ymax = upper, fill = compound_label),
    alpha = 0.25
  ) + 
  geom_step(
    data = km_df,
    aes(x = time, y = surv, color = compound_label),
    linewidth = 0.8
  ) +

  facet_wrap(~compound_label, ncol = 3) +

  scale_color_manual(values = colorRampPalette(c("#FFB3DE", "#C868C8", "#7B2D87", "#4A0072"))(12)) +
  scale_fill_manual(values = colorRampPalette(c("#FFB3DE", "#C868C8", "#7B2D87", "#4A0072"))(12)) +

  labs(
    title = "Krzywe przeżycia Kaplana-Meiera według związku",
    subtitle = "Przerywana szara linia = kontrola (compound 1) w każdym panelu",
    x = "Czas (dni)",
    y = "Prawdopodobieństwo przeżycia"
  ) +

  theme_minimal() +
  theme(
    plot.title = element_markdown(margin = margin(t = 10, r = 10, b = 10, l = 0)),
    legend.position = "none",
    strip.text = element_text(size = 8),
    plot.subtitle = element_text(color = "grey50")
  )

ggsave(km_plot,
  filename = "TeX/SAP/plots/km_compounds_facet.pdf",
  device = cairo_pdf, width = 10, height = 14, units = "in"
)

## 4 - Histogram of survival time --------------------------------------------

surv_time_hist <- ggplot(surv_rose, aes(x = time)) +
  geom_histogram(
    color = "#7B2D87",
    fill = "#C868C8",
    binwidth = 2,
    alpha = 0.6
  ) +
  geom_vline(
    xintercept = median(surv_rose$time, na.rm = TRUE),
    linetype = "solid",
    color = "tomato",
    linewidth = 1.25
  ) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(breaks = seq(0, 200, by = 25)) +
  labs(
    title = "Rozkład czasu przeżycia róż (<span style='font-family:mono'>time</span>)",
    subtitle = "Dane pilotażowe n = 720",
    x = "Czas przeżycia (dni)",
    y = "Liczebność"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_text(color = "grey30")
  )

surv_time_hist
ggsave(surv_time_hist,
  filename = "TeX/SAP/plots/surv_time_hist.pdf",
  device = cairo_pdf,
  height = 5, width = 7, units = "in"
)

## 5 - KM Survival curves -> garden, species --------------
fit_g <- survfit(Surv(time, can.be.presented) ~ garden, data = surv_rose)
garden_km <- ggsurvplot(fit_g,
  data = surv_rose, conf.int = TRUE,
  xlab = "Czas (dni)",
  ylab = "Prawdopodobieństwo przeżycia",
  legend.title = "Ogród",
  legend.labs = c("Ogród południowy ", "Ogród północny")
)
aggregate(time ~ garden, data = surv_rose, FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))
ggsave(
  filename = "TeX/SAP/plots/garden_km.pdf",
  plot = garden_km$plot,
  device = cairo_pdf,
  height = 5, width = 7, units = "in"
)
fit_g <- survfit(Surv(time, can.be.presented) ~ species, data = surv_rose)
species_km <- ggsurvplot(fit_g,
  data = surv_rose, conf.int = TRUE, 
  xlab = "Czas (dni)",
  ylab = "Prawdopodobieństwo przeżycia",
  legend.title = "Gatunek",
  legend.labs = c("T. floribunda ", "T. hybrid")
)
aggregate(time ~ species, data = surv_rose, FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))
ggsave(
  filename = "TeX/SAP/plots/species_km.pdf",
  plot = species_km$plot,
  device = cairo_pdf,
  height = 5, width = 7, units = "in"
)
## 6 - CLogLog ----------
fit_clog <- survfit(Surv(time, can.be.presented) ~ compound, data = surv_rose)

clog_df <- surv_summary(fit_clog, data = surv_rose) %>%
  mutate(
    compound_raw   = as.character(gsub("compound=", "", as.character(strata))),
    compound_label = factor(
      unname(compound_names[compound_raw]),
      levels = unname(compound_names)
    ),
    cloglog = log(-log(surv))
  ) %>%
  filter(surv > 0 & surv < 1 & time > 0)

cloglog_plot <- ggplot(clog_df, aes(x = log(time), y = cloglog, color = compound_label)) +
  geom_step(linewidth = 0.8) +
  scale_color_manual(
    values = colorRampPalette(c("#FFB3DE", "#C868C8", "#7B2D87", "#4A0072"))(12)
  ) +
  labs(
    title = "Wykres Complementary Log-Log według związku",
    x = "log(czas) [dni]",
    y = "log(–log(S(t)))",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13),
    legend.position = "bottom",
    legend.text     = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 3))

cloglog_plot

ggsave(cloglog_plot,
  filename = "TeX/plots/cloglog.pdf",
  device = cairo_pdf,
  width = 9, height = 6
)

m_adj <- survreg(Surv(time, can.be.presented) ~ factor(compound),
  data = filter(surv_rose, time > 0, !compound %in% c(2, 3, 10)),
  dist = "weibull"
)
cat("Adjusted (shared k, compound covariate) k =", round(1 / m_adj$scale, 2), "\n"

## 6.5 Test parametru baseline_diam ------------------
ggplot(data = surv_rose, aes(x = baseline_diam, y = time, color = compound)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

m_ch <- survreg(Surv(time, can.be.presented) ~ factor(compound) + scale(baseline_diam),
  data = filter(surv_rose, time > 0, !compound %in% c(2, 3, 10)), dist = "weibull")
summary(m_ch)$table