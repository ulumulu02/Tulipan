# Histogramy czasu przeżycia
# Histogram czasu przeżycia dla wszystkich obserwacji, oprocz zwiazkow szkodliwych
tempik <- surv_rose %>%  
  filter(!(compound %in% c(2,3,10)))
surv_time_hist <- ggplot(tempik, aes(x = time)) +
  geom_histogram(
    binwidth = 2,          
    fill = "slateblue2", 
    color = "slateblue4", 
    alpha = 0.6,
    boundary = 0            
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = seq(0, max(surv_rose$time, na.rm = TRUE), by = 5)) +
  labs(
    title = "Distribution of Survival Times",
    subtitle = "Counts of observations grouped by time intervals",
    x = "Survival Time (Days)",
    y = "Number of Subjects"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank() # Cleans up the background
  )
surv_time_hist

# Histogramy czasu przeżycia z podziałem na cztery grupy
 group_surv_time_hist <- surv_rose %>%
   mutate(Grupa = case_when(
     compound %in% c(2, 3, 10) ~ "Związki szkodliwe",
     compound %in% c(7, 9, 14) ~ "Związki umiarkowane",
     compound == 1             ~ "Kontrola - woda destylowana",
     TRUE                      ~ "Związki obiecujące"
   )) %>%
   ggplot(aes(x = time, fill = Grupa)) +
   geom_histogram(binwidth = 1, color = "white", alpha = 0.8) +
   facet_wrap(~ Grupa, ncol = 2) +
   labs(title = "Rozkład czasu przeżycia róż dla czterech podgrup (<span style='font-family:mono'>category</span>)",
        x = "Czas przeżycia (dni)",
        y = "Liczebność") +      
   theme_minimal() +
    theme(
    plot.title = element_markdown(),
    legend.position = "none")
    
ggsave(group_surv_time_hist,
       filename = "TeX/plots/group_surv_time_hist.pdf",
       device = cairo_pdf,
       width = 7,
       height = 5
       )
