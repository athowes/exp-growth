library(rstan)
library(tidyverse)
library(bayesplot)

cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

compiled <- rstan::stan_model("4pl.stan")
cycles <- 20

sim <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 0, n = rep(0, cycles))
)

sim_summary <- summary(sim)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  #' I've hard coded this, better to be automatic
  #' Somewhere in sim@par_dims
  mutate(t = c(1, 1, 1, 1, 1:cycles, 1:cycles, 1))

lambda_prior <- sim_summary %>%
  filter(substr(rowname, 1, 6) == "lambda")

pdf("lambda-prior.pdf", h = 5, w = 6.25)

lambda_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(col = cbpalette[1], linetype = "dashed") +
  geom_ribbon(alpha = 0.5, fill = cbpalette[1]) +
  theme_minimal() +
  lims(y = c(0, 10)) +
  labs(x = "t", y = "fluorescence")

dev.off()

n_sim <- rstan::extract(sim, )$n_sim %>%
  as.matrix()

#' Example of one simulated data from the model
n_sim[1, ]
df_one <- data.frame(n = n_sim[1, ], t = 1:cycles)

pdf("lambda-prior-data.pdf", h = 5, w = 6.25)

lambda_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(col = cbpalette[1], linetype = "dashed") +
  geom_ribbon(alpha = 0.5, fill = cbpalette[1]) +
  geom_point(data = df_one, aes(x = t, y = n), inherit.aes = FALSE) +
  theme_minimal() +
  lims(y = c(0, 10)) +
  labs(x = "t", y = "fluorescence")

dev.off()

#' Let's fit the model to that data and see what it thinks the posterior would be
fit <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 1, n = df_one$n)
)

fit_summary <- summary(fit)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  #' I've hard coded this, better to be automatic
  #' Somewhere in sim@par_dims
  mutate(t = c(1, 1, 1, 1, 1:cycles, 1:cycles, 1))

lambda_posterior <- fit_summary %>%
  filter(substr(rowname, 1, 6) == "lambda")

pdf("lambda-prior-data-posterior.pdf", h = 5, w = 6.25)

lambda_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(linetype = "dashed") +
  geom_ribbon(alpha = 0.5, fill = cbpalette[1]) +
  geom_point(data = df_one, aes(x = t, y = n), inherit.aes = FALSE) +
  geom_line(
    data = lambda_posterior, aes(x = t, y = mean), inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = lambda_posterior, aes(x = t, ymin = `2.5%`, ymax = `97.5%`), inherit.aes = FALSE,
    fill = cbpalette[3], alpha = 0.5
  ) +
  theme_minimal() +
  lims(y = c(0, 20)) +
  labs(x = "t", y = "fluorescence")

dev.off()
