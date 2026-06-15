# 2. Create the faceted histogram
surv_rose_labeled <- surv_rose %>%
  mutate(compound_name = factor(compound, levels = 1:15, labels = compound_names))
facet_histograms_density <- ggplot(surv_rose_labeled, aes(x = time)) +
  # 1. Use after_stat(density) so the histogram and line use the same scale
  geom_histogram(
    aes(y = after_stat(density)), 
    bins = 20, 
    fill = "slateblue2", 
    color = "slateblue4", 
    alpha = 0.4, # Lower alpha so the line stands out
    boundary = 0
  ) +
  # 2. Add the density line
  geom_density(color = "red", linewidth = 1) + 
  
  facet_wrap(~compound_name, ncol = 5) + 
  scale_x_continuous(breaks = seq(0, 30, by = 10)) + 
  
  labs(
    title = "Survival Density Distributions",
    subtitle = "Red lines indicate density; look for 'two humps' for bimodality",
    x = "Survival Time (Days)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank()
  )

facet_histograms_density

ggplot(surv_rose[surv_rose$compound == 10,],
       aes(x=time)) +
  geom_histogram(bins = 5) +
  facet_wrap(~species)
