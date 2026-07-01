# extract_results.R — pull tables from fit_main_fixed.rds for the final report

library(brms)
library(dplyr)
library(purrr)

fit <- readRDS("models/fit_main_fixed.rds")

# ── 1. Convergence diagnostics ───────────────────────────────────────────────
s <- summary(fit)
params_fixed  <- rownames(s$fixed)
params_spec   <- rownames(s$spec_pars)  # shape

rhat_df <- data.frame(
  param     = c(params_fixed, params_spec),
  rhat      = round(c(s$fixed[, "Rhat"],       s$spec_pars[, "Rhat"]),       4),
  bulk_ess  = round(c(s$fixed[, "Bulk_ESS"],   s$spec_pars[, "Bulk_ESS"]),   0),
  tail_ess  = round(c(s$fixed[, "Tail_ESS"],   s$spec_pars[, "Tail_ESS"]),   0)
)
cat("\n=== CONVERGENCE ===\n"); print(rhat_df, row.names = FALSE)

# ── 2. Posterior summary (fixed effects + shape) ─────────────────────────────
post_df <- data.frame(
  param  = c(params_fixed, params_spec),
  median = round(c(s$fixed[, "Estimate"],     s$spec_pars[, "Estimate"]),     3),
  q025   = round(c(s$fixed[, "l-95% CI"],     s$spec_pars[, "l-95% CI"]),     3),
  q975   = round(c(s$fixed[, "u-95% CI"],     s$spec_pars[, "u-95% CI"]),     3)
)
cat("\n=== POSTERIOR SUMMARY ===\n"); print(post_df, row.names = FALSE)

# ── 3. Compound ranking + survival probs ─────────────────────────────────────
draws <- as_draws_df(fit)

compound_lvls <- c("1","4","5","6","7","8","9","11","12","13","14","15")

get_lambda <- function(c_lvl, draws) {
  k   <- draws$shape
  b0  <- draws$b_Intercept
  a_c <- if (c_lvl == "1") 0 else draws[[paste0("b_compound", c_lvl)]]
  exp(b0 + a_c) / gamma(1 + 1/k)
}

results <- map_dfr(compound_lvls, function(c_lvl) {
  k        <- draws$shape
  lambda_c <- get_lambda(c_lvl, draws)
  median_c <- lambda_c * log(2)^(1/k)
  tibble(
    compound    = as.integer(c_lvl),
    post_median = round(median(median_c), 2),
    q025        = round(quantile(median_c, 0.025), 2),
    q975        = round(quantile(median_c, 0.975), 2),
    P_surv_10   = round(mean(exp(-(10 / lambda_c)^k)), 3),
    P_surv_22   = round(mean(exp(-(22 / lambda_c)^k)), 3)
  )
}) |> arrange(desc(post_median))

cat("\n=== COMPOUND RANKING ===\n"); print(results, n = Inf)
