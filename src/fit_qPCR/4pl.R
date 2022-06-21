cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

compiled <- rstan::stan_model("4pl.stan")
cycles <- 20

sim <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 0, f = rep(0, cycles))
)

sim_summary <- summary(sim)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  #' I've hard coded this, better to be automatic
  #' Somewhere in sim@par_dims
  mutate(t = c(1, 1, 1, 1, 1, 1, 1, 1:cycles, 1:cycles, 1))

n_prior <- sim_summary %>%
  filter(substr(rowname, 1, 1) == "n")

pdf("n-prior.pdf", h = 5, w = 6.25)

n_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(col = cbpalette[1], linetype = "dashed") +
  geom_ribbon(alpha = 0.5, fill = cbpalette[1]) +
  theme_minimal() +
  lims(y = c(0, 12)) +
  labs(x = "t", y = "output")

dev.off()

f_sim <- rstan::extract(sim, )$f_sim %>%
  as.matrix()

#' Example of one simulated data from the model
f_sim[1, ]
df_one <- data.frame(f = f_sim[1, ], t = 1:cycles)

pdf("n-prior-data.pdf", h = 5, w = 6.25)

n_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(col = cbpalette[1], linetype = "dashed") +
  geom_ribbon(alpha = 0.5, fill = cbpalette[1]) +
  geom_point(data = df_one, aes(x = t, y = f), inherit.aes = FALSE) +
  theme_minimal() +
  lims(y = c(0, 20)) +
  labs(x = "t", y = "fluorescence")

dev.off()

#' Let's fit the model to that data and see what it thinks the posterior would be
fit <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 1, f = df_one$f)
)

fit_summary <- summary(fit)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  #' I've hard coded this, better to be automatic
  #' Somewhere in sim@par_dims
  mutate(t = c(1, 1, 1, 1, 1, 1, 1, 1:cycles, 1:cycles, 1))

n_posterior <- fit_summary %>%
  filter(substr(rowname, 1, 1) == "n")

pdf("n-prior-data-posterior.pdf", h = 5, w = 6.25)

n_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(linetype = "dashed") +
  geom_ribbon(alpha = 0.5, fill = cbpalette[1]) +
  geom_point(data = df_one, aes(x = t, y = f), inherit.aes = FALSE) +
  geom_line(
    data = n_posterior, aes(x = t, y = mean), inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = n_posterior, aes(x = t, ymin = `2.5%`, ymax = `97.5%`), inherit.aes = FALSE,
    fill = cbpalette[3], alpha = 0.5
  ) +
  theme_minimal() +
  labs(x = "t", y = "fluorescence")

dev.off()
