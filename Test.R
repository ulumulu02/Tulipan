## ============================================================
## Compare Weibull vs Log-logistic vs Log-normal AFT models
## per compound, using AIC/BIC + diagnostic plots
## ============================================================

# install.packages(c("survival", "flexsurv", "dplyr", "ggplot2"))
library(survival)
library(flexsurv)
library(dplyr)
library(ggplot2)
data <- surv_rose
## ------------------------------------------------------------
## 1. Data format expected
## ------------------------------------------------------------
## A data frame `df` with columns:
##   compound : factor/character, compound name (e.g. "1 Distilled Water")
##   time     : numeric, time-to-event or censoring time (days)
##   can.be.presented   : 1 = event observed, 0 = right-censored
##
## Example placeholder (REPLACE with your real data, e.g. via read.csv):
## df <- read.csv("your_survival_data.csv")

## ------------------------------------------------------------
## 2. Fit all three AFT models for a single compound
## ------------------------------------------------------------
fit_aft_models <- function(data) {
  surv_obj <- Surv(data$time, data$can.be.presented)

  fit_weibull <- flexsurvreg(surv_obj ~ 1, data = data, dist = "weibull")
  fit_loglogistic <- flexsurvreg(surv_obj ~ 1, data = data, dist = "llogis")
  fit_lognormal <- flexsurvreg(surv_obj ~ 1, data = data, dist = "lnorm")

  list(
    weibull = fit_weibull,
    loglogistic = fit_loglogistic,
    lognormal = fit_lognormal
  )
}

## ------------------------------------------------------------
## 3. Extract AIC/BIC into a tidy comparison table
## ------------------------------------------------------------
compare_models <- function(fits, n_obs) {
  get_stats <- function(fit, model_name) {
    ll <- fit$loglik
    k <- length(fit$coefficients)
    aic <- AIC(fit)
    bic <- -2 * ll + k * log(n_obs)
    data.frame(
      model = model_name,
      logLik = ll,
      k = k,
      AIC = aic,
      BIC = bic
    )
  }

  bind_rows(
    get_stats(fits$weibull, "Weibull"),
    get_stats(fits$loglogistic, "Log-logistic"),
    get_stats(fits$lognormal, "Log-normal")
  ) %>%
    arrange(AIC)
}

## ------------------------------------------------------------
## 4. Run per-compound comparison across the whole dataset
## ------------------------------------------------------------
run_all_compounds <- function(df) {
  compounds <- unique(df$compound)
  results <- list()

  for (cmpd in compounds) {
    sub <- df %>% filter(compound == cmpd)

    # skip compounds with too few events to fit 2-parameter models reliably
    if (sum(sub$can.be.presented) < 3) {
      message("Skipping ", cmpd, ": fewer than 3 observed events")
      next
    }

    fits <- tryCatch(fit_aft_models(sub), error = function(e) {
      message("Fit failed for ", cmpd, ": ", e$message)
      NULL
    })
    if (is.null(fits)) next

    tab <- compare_models(fits, nrow(sub))
    tab$compound <- cmpd
    tab$best_model <- tab$model[which.min(tab$AIC)]
    results[[cmpd]] <- tab
  }

  bind_rows(results)
}

## ------------------------------------------------------------
## 5. Summarize: which distribution wins most often
## ------------------------------------------------------------
summarize_winners <- function(results) {
  results %>%
    filter(model == best_model) %>%
    count(model, name = "n_compounds_best") %>%
    arrange(desc(n_compounds_best))
}

## ------------------------------------------------------------
## 6. Diagnostic overlay plot: KM curve + 3 fitted survival curves
## ------------------------------------------------------------
plot_compound_fit <- function(data, fits, compound_name) {
  km <- survfit(Surv(time, can.be.presented) ~ 1, data = data)

  time_grid <- seq(min(data$time), max(data$time), length.out = 200)

  pred_df <- bind_rows(
    data.frame(time = time_grid,
               surv = summary(fits$weibull, t = time_grid)[[1]]$est,
               model = "Weibull"),
    data.frame(time = time_grid,
               surv = summary(fits$loglogistic, t = time_grid)[[1]]$est,
               model = "Log-logistic"),
    data.frame(time = time_grid,
               surv = summary(fits$lognormal, t = time_grid)[[1]]$est,
               model = "Log-normal")
  )

  km_df <- data.frame(time = km$time, surv = km$surv)

  ggplot() +
    geom_step(data = km_df, aes(time, surv), color = "black", linewidth = 0.8) +
    geom_line(data = pred_df, aes(time, surv, color = model), linewidth = 0.7) +
    labs(title = paste("KM vs fitted AFT models:", compound_name),
         x = "Time (days)", y = "Survival probability", color = "Model") +
    theme_minimal()
}

## ------------------------------------------------------------
## 7. Example usage
## ------------------------------------------------------------
results <- run_all_compounds(data)
print(results)
print(summarize_winners(results))

# Plot for one specific compound
sub <- data %>% filter(compound == "7 Four in December")
fits <- fit_aft_models(sub)
plot_compound_fit(sub, fits, "7 Four in December")

# Save full comparison table
write.csv(results, "aft_model_comparison.csv", row.names = FALSE)