ggplot(pilot_estimates, aes(x = mu_hat)) +
  geom_histogram(bins = 30)

ggplot(pilot_estimates, aes(x = factor(compound), y = mu_hat)) +
  geom_boxplot()

k_sim  <- rgamma(5000, shape = k_alpha, rate = k_beta)
mu_sim <- exp(rnorm(5000, logmu_mean, logmu_sd_inflated))

t_sim <- rweibull(5000, shape = k_sim, scale = mu_sim)

ggplot(data.frame(t_sim), aes(x = t_sim)) +
  geom_histogram(bins = 50) +
  labs(title = "Prior predictive distribution of survival time")

mean(k_sim)
mean(mu_sim)


weibull_data <- rweibull(n = 1000, shape = mean(k_sim), scale = mean(mu_sim))


density_est <- density(t_sim)

ggplot(data.frame(t_sim), aes(x = t_sim)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 color = "black", fill = "skyblue") +
  geom_line(data = data.frame(x = density_est$x, y = density_est$y),
            aes(x = x, y = y), color = "red", linewidth = 1) +
  labs(title = "Mixture Distribution (Empirical Density)",
       x = "Values", y = "Density") +
  theme_minimal()



