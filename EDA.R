
rose <- read.csv2("roses_analysis2.csv")
compound_names <- c(
  "1 Distilled Water", "2 Apathic Acid", "3 Beerse Brew",
  "4 Concentrate of Caducues", "5 Distillate of Discovery",
  "6 Essence of Epiphaneia", "7 Four in December",
  "8 Granules of Geheref", "9 Kar-Hamel Mooh", "10 Lucifer's Liquid",
  "11 Noospherol", "12 Oil of John's Son", "13 Power of Perlimpinpin",
  "14 Spirit of Scienza", "15 Zest of Zen"
)

rose %>%
  group_by(compound) %>%
  summarise(
    median = median(time),
    mean = mean(time)
  ) %>%
  arrange(mean)
# Data prep -----------------
rose <- as.data.frame(rose)
surv_rose <- rose %>%
  left_join(
    rose %>%
      filter(time == 0) %>%
      select(id, baseline_diam = diam),
    by = "id"
  ) %>%
  group_by(id) %>%
  slice_max(order_by = time, n = 1, with_ties = FALSE) %>%
  ungroup()

# Usunięcie szkodliwych związków
surv_rose <- surv_rose %>%
  mutate(can.be.presented = ifelse(can.be.presented == FALSE, 1, 0))  |> 
  filter(!(compound %in% c(2,3,10)))

head(surv_rose)

# Categories of compounds -> control 0, harmful 1, moderate 2, promising 3
# surv_rose <- surv_rose %>%
#   mutate(category = case_when(
#     compound == 1 ~ 0,
#     compound %in% c(2, 3, 10) ~ 1,
#     compound %in% c(7, 9, 14) ~ 2,
#     compound %in% c(4, 5, 6, 8, 11, 12, 13, 15) ~ 3
#   ))

# Exploratory Data Analysis (EDA) ----------------------------------------------

## 1 - Summary table of compounds ------------------------------------------

sum_table_time <- surv_rose %>%
  group_by(compound) %>%
  summarise(
    no_flowers = n(),
    mean_surv_time = mean(time),
    med_surv_time = median(time),
    sd_surv_time = sd(time)
  )
sum_table_time
latex_table <- xtable(sum_table_time,
  align = c("c", "c", "c", "c", "c", "c")
)
print(
  latex_table,
  type = "latex",
  include.rownames = FALSE
)

## 2 - Boxplot of survival time by compound ------------------------------
compound_boxplot <- surv_rose %>%
  mutate(Grupa = case_when(
    compound %in% c(2, 3, 10) ~ "Związki szkodliwe",
    compound %in% c(7, 9, 14) ~ "Związki umiarkowane",
    compound == 1 ~ "Kontrola - woda destylowana",
    TRUE ~ "Związki obiecujące"
  )) %>%
  ggplot(aes(x = factor(compound), y = time, fill = Grupa)) +
  geom_boxplot(alpha = 0.3) +
  coord_flip() +
  scale_x_discrete(labels = compound_names) +
  scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  labs(
    title = "Rozkład przeżycia kwiatów według zastosowanego związku",
    x = NULL,
    y = "Czas przeżycia (dni)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

compound_boxplot

ggsave(compound_boxplot,
  filename = "TeX/plots/boxplot_compounds.pdf",
  device = cairo_pdf,
  width = 9,
  height = 6
)

## 3 - KM Survival curves -----------------------------------------------------
# KM for control (compound 1) — will be overlaid on every facet
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

ggsave("TeX/plots/km_compounds_facet.pdf",
  device = cairo_pdf, width = 10, height = 14, units = "in"
)

## 4 - Histogram of survival time --------------------------------------------

surv_time_hist <- ggplot(surv_rose, aes(x = time)) +
  geom_histogram(color = "slateblue4", bins = 15, fill = "slateblue2", alpha = 0.6) +
  labs(
    title = "Histogram czasu przeżycia róż",
    x = "Czas przeżycia (dni)",
    y = "Liczebność"
  ) +
  theme_minimal()
surv_time_hist

surv_time_hist <- ggplot(surv_rose, aes(x = time)) +
  geom_histogram(
    color = "slateblue4",
    fill = "slateblue2",
    bins = 15,
    alpha = 0.6
  ) +
  geom_vline(
    xintercept = median(surv_rose$time, na.rm = TRUE),
    linetype = "solid",
    color = "yellow",
    linewidth = 1.25
  ) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(breaks = seq(0, 200, by = 25)) +
  labs(
    title = "Rozkład czasu przeżycia róż (<span style='font-family:mono'>time</span>)",
    subtitle = "Dane pilotażowe n = 900",
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
  filename = "TeX/plots/surv_time_hist.pdf",
  device = cairo_pdf,
  height = 5, width = 7, units = "in"
)

## 5 - KM Survival curves -> garden, species --------------
fit_g <- survfit(Surv(time, can.be.presented) ~ garden, data = surv_rose)
garden_km <- ggsurvplot(fit_g,
  data = surv_rose, conf.int = TRUE, pval = TRUE,
  xlab = "Czas (dni)",
  ylab = "Prawdopodobieństwo przeżycia",
  legend.title = "Ogród",
  legend.labs = c("Ogród południowy ", "Ogród północny")
)
aggregate(time ~ garden, data = surv_rose, FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))
ggsave(
  filename = "TeX/plots/garden_km.pdf",
  plot = garden_km$plot,
  device = cairo_pdf,
  height = 5, width = 7, units = "in"
)

fit_g <- survfit(Surv(time, can.be.presented) ~ species, data = surv_rose)
species_km <- ggsurvplot(fit_g,
  data = surv_rose, conf.int = TRUE, pval = TRUE,
  xlab = "Czas (dni)",
  ylab = "Prawdopodobieństwo przeżycia",
  legend.title = "Gatunek",
  legend.labs = c("T. floribunda ", "T. hybrid")
)
aggregate(time ~ species, data = surv_rose, FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))
ggsave(
  filename = "TeX/plots/species_km.pdf",
  plot = species_km$plot,
  device = cairo_pdf,
  height = 5, width = 7, units = "in"
)
## 6 - CLogLog ----------
fit_clog <- survfit(Surv(time, can.be.presented) ~ category, data = surv_rose)

clog_df <- data.frame(
  time  = fit_clog$time,
  surv  = fit_clog$surv,
  group = rep(names(fit_clog$strata), fit_clog$strata)
) %>%
  filter(surv > 0 & surv < 1 & time > 0) %>%
  mutate(
    cloglog = log(-log(surv)),
    Grupa = case_when(
      group == "category=0" ~ "Kontrola - woda destylowana",
      group == "category=1" ~ "Związki szkodliwe",
      group == "category=2" ~ "Związki umiarkowane",
      group == "category=3" ~ "Związki obiecujące"
    )
  )

cloglog_plot <- ggplot(clog_df, aes(x = log(time), y = cloglog, color = Grupa)) +
  geom_step(linewidth = 1) +
  labs(
    title = "Wykres Complementary Log-Log według grupy związków",
    x = "log(czas) [dni]",
    y = "log(-log(S(t)))",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13),
    plot.subtitle   = element_text(color = "grey50"),
    legend.position = "bottom"
  )

cloglog_plot

ggsave(cloglog_plot,
  filename = "TeX/plots/cloglog.pdf",
  device = cairo_pdf,
  width = 7, height = 5
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