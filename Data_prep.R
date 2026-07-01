# =====================================================
# Statystyki opisowe dla głównych danych
# =====================================================
library(dplyr)
rose <- read.csv2("roses_analysis2_final.csv")
rose <- as.data.frame(rose)
# Data prep -----------------
main_rose <- rose %>%
    filter(!(compound %in% c(2,3,10)))  |> 
  left_join(
    rose %>%
      filter(time == 0) %>%
      select(id, baseline_diam = diam),
    by = "id"
  ) %>%
  group_by(id) %>%
  slice_max(order_by = time, n = 1, with_ties = FALSE) %>%
  ungroup()
main_rose <- main_rose %>%
  mutate(can.be.presented = ifelse(can.be.presented == FALSE, 1, 0))  |> 

