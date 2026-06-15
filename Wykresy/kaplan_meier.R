library(patchwork)
# Compound name lookup

# Kaplan-Meier z podziałem na 3 grupy
compound_names <- c(
  "1 Distilled Water", "2 Apathic Acid", "3 Beerse Brew",
  "4 Concentrate of Caducues", "5 Distillate of Discovery",
  "6 Essence of Epiphaneia", "7 Four in December",
  "8 Granules of Geheref", "9 Kar-Hamel Mooh", "10 Lucifer's Liquid",
  "11 Noospherol", "12 Oil of John's Son", "13 Power of Perlimpinpin",
  "14 Spirit of Scienza", "15 Zest of Zen"
)

get_km_df <- function(compounds) {
  df <- surv_rose %>% filter(compound %in% compounds)
  fit <- survfit(Surv(time, can.be.presented) ~ compound, data = df)
  surv_summary(fit, data = df) %>%
    mutate(
      compound_id    = as.integer(gsub("compound=", "", strata)),
      compound_label = factor(compound_names[compound_id],
                              levels = compound_names[compounds])
    )
}

control_fit <- survfit(Surv(time, can.be.presented) ~ 1,
                       data = filter(surv_rose, compound == 1))
control_df <- data.frame(
  time  = control_fit$time,
  surv  = control_fit$surv,
  upper = control_fit$upper,
  lower = control_fit$lower
)

groups <- list(
  "Plot 1 — Harmful (median ≤ 3 days)"         = c(2, 3, 10),
  "Plot 2 — Moderate (median 13-15 days)"       = c(7, 9, 14),
  "Plot 3a — Promising I (median ~17-18 days)"  = c(4, 6, 8, 15),
  "Plot 3b — Promising II (median ~19-21 days)" = c(5, 11, 12, 13)
)

group_colors <- list(
  "Plot 1 — Harmful (median ≤ 3 days)"         = c("#E41A1C", "#FF7F00", "#984EA3"),
  "Plot 2 — Moderate (median 13-15 days)"       = c("#4DAF4A", "#377EB8", "#A65628"),
  "Plot 3a — Promising I (median ~17-18 days)"  = c("#1B9E77", "#7570B3", "#D95F02", "#E7298A"),
  "Plot 3b — Promising II (median ~19-21 days)" = c("#66A61E", "#E6AB02", "#A6761D", "#666666")
)

# Force same x axis across all plots
make_plot <- function(group_name, compounds, colors) {
  km_df <- get_km_df(compounds)
  
  ggplot() +
    geom_ribbon(data = km_df,
                aes(x = time, ymin = lower, ymax = upper,
                    fill = compound_label),
                alpha = 0.2) +
    geom_step(data = km_df,
              aes(x = time, y = surv, color = compound_label),
              linewidth = 0.9) +
    geom_ribbon(data = control_df,
                aes(x = time, ymin = lower, ymax = upper),
                fill = "grey50", alpha = 0.15) +
    geom_step(data = control_df,
              aes(x = time, y = surv),
              color = "grey30", linewidth = 0.9, linetype = "dashed") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values  = colors) +
    scale_x_continuous(limits = c(0, 30)) +  # <-- same x range for all
    scale_y_continuous(limits = c(0, 1)) +   # <-- same y range for all
    labs(title = group_name, x = "Days", y = "Survival Probability",
         color = "Compound", fill = "Compound") +
    theme_minimal() +
    theme(
      plot.title       = element_text(size = 11, face = "bold"),
      legend.position  = "right",
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm")
    )
}


plots <- mapply(make_plot,
                names(groups), groups, group_colors,
                SIMPLIFY = FALSE)

final_plot <- plots[[1]] / plots[[2]] / plots[[3]] / plots[[4]] +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(
    title    = "Survival Curves by Compound Group",
    subtitle = "Dashed grey = Compound 1 (Distilled Water) reference in each panel",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(color = "grey50")
    )
  )
final_plot
ggsave("km_compounds.png", plot = final_plot, width = 10, height = 16, dpi = 300)
