# =====================================================
# Statystyki opisowe dla głównych danych
# =====================================================
library(dplyr)
library(survival)
d <- read.csv2("roses_analysis2_final.csv")
s <- d |>
  group_by(id) |>
  slice_max(order_by = time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    event    = ifelse(can.be.presented == FALSE, 1L, 0L),
    time     = pmax(as.numeric(time), 0.5),
    compound = factor(compound)
  ) |>
  filter(!(compound %in% c(2, 3, 10))) |>
  droplevels()

lvls <- c("1","4","5","6","7","8","9","11","12","13","14","15")
for (c_lvl in lvls) {
  sub <- s[s$compound == c_lvl, ]
  km  <- survfit(Surv(time, event) ~ 1, data = sub)
  med <- summary(km)$table["median"]
  cat(sprintf("%-2s | n=%d | ev=%d | cens=%d (%.1f%%) | KM_med=%.1f\n",
              c_lvl, nrow(sub), sum(sub$event==1), sum(sub$event==0),
              mean(sub$event==0)*100, med))
}
