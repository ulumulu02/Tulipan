library(dplyr)
rm(list=ls())
rose <- read.csv2("TeX/screening/code/roses_analysis1.csv")
compound_names <- c(
  "1 Distilled Water", "2 Apathic Acid", "3 Beerse Brew",
  "4 Concentrate of Caducues", "5 Distillate of Discovery",
  "6 Essence of Epiphaneia", "7 Four in December",
  "8 Granules of Geheref", "9 Kar-Hamel Mooh", "10 Lucifer's Liquid",
  "11 Noospherol", "12 Oil of John's Son", "13 Power of Perlimpinpin",
  "14 Spirit of Scienza", "15 Zest of Zen"
)
rose <- as.data.frame(rose)

rose %>%
  group_by(compound) %>%
  summarise(
    median = median(time),
    mean = mean(time)
  ) %>%
  arrange(median)
# Data prep -----------------
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
surv_rose <- surv_rose %>%
  mutate(can.be.presented = ifelse(can.be.presented == FALSE, 1, 0))

surv_rose  |> 
  group_by(compound)  |> 
  summarise(median = median(time))
# head(surv_rose)

write.csv(surv_rose, "TeX/screening/code/surv_rose.csv")
# surv_rose  |> 
#   group_by(compound) %>%
#   summarise(
#     median = median(time),
#     mean = mean(time)
#   )  |> 
#   arrange(median)

