library(rstan)
library(tidyverse)

compiled <- rstan::stan_model("4pl.stan")
sim <- rstan::sampling(compiled, data = list(T = 10, flag_run_estimation = 0))

summary_lambda <- summary(sim)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  filter(substr(rowname, 1, 6) == "lambda") %>%
  mutate(id = 1:n())

ggplot(summary_lambda, aes(x = id, y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_line() +
  geom_ribbon(alpha = 0.5) +
  theme_minimal() +
  labs(x = "t", y = "latent quantity")
