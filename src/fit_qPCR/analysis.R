set.seed(3)

#' Colour-blind friendly colours
cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' 4PL (four parameter logistic)

#' Based on amplification curve data I make the following judgements
#' * The right-hand asymptote is around 8 on the "Delta Rn" scale
#' * The number of cycles is around 40
#' * The curves are relatively steep when they start taking off
#' * We could fix the left-hand asymptote to zero

#' Simulate data from the prior
compiled <- rstan::stan_model("4pl/4pl.stan")
cycles <- 40

sim <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 0, f = rep(0, cycles))
)

stan_par_index <- function(fit) {
  fit@par_dims %>%
    lapply(function(x) rep(1:ifelse(length(x) == 0, 1, x[1]))) %>%
    unlist()
}

sim_summary <- rstan::summary(sim)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  mutate(t = stan_par_index(sim))

f_prior <- sim_summary %>%
  filter(str_starts(rowname, "f_sim"))

pdf("4pl/prior.pdf", h = 3, w = 6.25)

plot_prior <- f_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(size = 1, col = cbpalette[1]) +
  geom_ribbon(alpha = 0.2, fill = cbpalette[1]) +
  theme_minimal() +
  labs(x = "Cycle number", y = "Fluoresence (Delta Rn)")

plot_prior

dev.off()

pdf("4pl/prior-stan.pdf", h = 8, w = 6.25)

rstan::stan_plot(sim, pars = sim@model_pars) + theme_minimal()
rstan::stan_trace(sim, pars = sim@model_pars) + theme_minimal()

dev.off()

f_sim <- rstan::extract(sim)$f_sim %>%
  as.matrix()

#' Simulated curves from the prior
n_curves <- 10
chosen_curves <- seq(from = 1, to = (100 * (n_curves - 1)) + 1, by = 100)
f_sim_curves <- f_sim[chosen_curves, 1:40] %>%
  t() %>%
  as.data.frame() %>%
  gather(
    key = "sim",
    value = "f"
  ) %>%
  mutate(
    sim = as.numeric(str_remove(sim, "V")),
    t = rep(1:40, length(chosen_curves)),
    observed = as.factor(c(rep(1, 40), rep(0, 40 * (n_curves - 1))))
  )

pdf("4pl/prior-data.pdf", h = 3, w = 6.25)

plot_prior_data <- plot_prior +
  geom_line(
    data = f_sim_curves, aes(x = t, y = f, group = sim, col = observed),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  scale_color_manual(values = c("grey80", "black")) +
  labs(col = "Observed?")

plot_prior_data

dev.off()

#' Suppose we observed one of the simulated data
obs <- f_sim_curves %>%
  filter(observed == 1)

fit <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 1, f = obs$f)
)

fit_summary <- rstan::summary(fit)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  mutate(t = stan_par_index(fit))

f_posterior <- fit_summary %>%
  filter(str_starts(rowname, "f_sim"))

pdf("4pl/prior-data-posterior.pdf", h = 3, w = 6.25)

plot_prior_data_posterior <- plot_prior_data +
  geom_line(
    data = f_posterior, aes(x = t, y = mean), inherit.aes = FALSE,
    col = cbpalette[3], size = 1
  ) +
  geom_ribbon(
    data = f_posterior, aes(x = t, ymin = `2.5%`, ymax = `97.5%`), inherit.aes = FALSE,
    fill = cbpalette[3], alpha = 0.5
  )

plot_prior_data_posterior

dev.off()

pdf("4pl/posterior-stan.pdf", h = 8, w = 6.25)

rstan::stan_plot(fit, pars = fit@model_pars) + theme_minimal()
rstan::stan_trace(fit, pars = fit@model_pars) + theme_minimal()

dev.off()


#' 4PL with processing step

#' Follow the same steps as above
compiled <- rstan::stan_model("4pl-processing/4pl-processing.stan")
cycles <- 40

sim <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 0, f = rep(0, cycles))
)

sim_summary <- rstan::summary(sim)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  mutate(t = stan_par_index(sim))

f_prior <- sim_summary %>%
  filter(str_starts(rowname, "f_sim"))

pdf("4pl-processing/prior.pdf", h = 3, w = 6.25)

plot_prior <- f_prior %>%
  ggplot(aes(x = t, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line(size = 1, col = cbpalette[1]) +
  geom_ribbon(alpha = 0.2, fill = cbpalette[1]) +
  theme_minimal() +
  labs(x = "Cycle number", y = "Fluoresence (Delta Rn)")

plot_prior

dev.off()

pdf("4pl-processing/prior-stan.pdf", h = 8, w = 6.25)

rstan::stan_plot(sim, pars = sim@model_pars) + theme_minimal()
rstan::stan_trace(sim, pars = sim@model_pars) + theme_minimal()

dev.off()

f_sim <- rstan::extract(sim)$f_sim %>%
  as.matrix()

#' Simulated curves from the prior
n_curves <- 10
chosen_curves <- seq(from = 1, to = (100 * (n_curves - 1)) + 1, by = 100)
f_sim_curves <- f_sim[chosen_curves, 1:40] %>%
  t() %>%
  as.data.frame() %>%
  gather(
    key = "sim",
    value = "f"
  ) %>%
  mutate(
    sim = as.numeric(str_remove(sim, "V")),
    t = rep(1:40, length(chosen_curves)),
    observed = as.factor(c(rep(1, 40), rep(0, 40 * (n_curves - 1))))
  )

pdf("4pl-processing/prior-data.pdf", h = 3, w = 6.25)

plot_prior_data <- plot_prior +
  geom_line(
    data = f_sim_curves, aes(x = t, y = f, group = sim, col = observed),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  scale_color_manual(values = c("grey80", "black")) +
  labs(col = "Observed?")

plot_prior_data

dev.off()

obs <- f_sim_curves %>%
  filter(observed == 1)

fit <- rstan::sampling(
  compiled,
  data = list(T = cycles, flag_run_estimation = 1, f = obs$f)
)

fit_summary <- rstan::summary(fit)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  mutate(t = stan_par_index(fit))

f_posterior <- fit_summary %>%
  filter(str_starts(rowname, "f_sim"))

pdf("4pl-processing/prior-data-posterior.pdf", h = 3, w = 6.25)

plot_prior_data_posterior <- plot_prior_data +
  geom_line(
    data = f_posterior, aes(x = t, y = mean), inherit.aes = FALSE,
    col = cbpalette[3], size = 1
  ) +
  geom_ribbon(
    data = f_posterior, aes(x = t, ymin = `2.5%`, ymax = `97.5%`), inherit.aes = FALSE,
    fill = cbpalette[3], alpha = 0.5
  )

plot_prior_data_posterior

dev.off()

pdf("4pl-processing/posterior-stan.pdf", h = 8, w = 6.25)

rstan::stan_plot(fit, pars = fit@model_pars) + theme_minimal()
rstan::stan_trace(fit, pars = fit@model_pars) + theme_minimal()

dev.off()
