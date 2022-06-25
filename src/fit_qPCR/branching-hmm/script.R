#' Model code based on suggestion by Perry de Valpine
library(nimble)

branching_hmm <- nimbleCode({
  N_cont ~ dunif(0, 1000)
  N <- round(N_cont)
  r ~ dunif(0, 1)
  E ~ dunif(0, 1)
  sigma_f ~ T(dnorm(0, sd = 1), 0, Inf)
  alpha ~ dunif(0, 100)
  n[1]~ dbinom(N_0, r)
  for(i in 2:C) {
    n_new[i] ~ dbinom(n[i - 1], E)
    n[i] <- n[i - 1] + n_new[i]
    f[i] ~ dnorm(alpha * n[i], sd = sigma_f)
  }
})

sim_model <- nimbleModel(
  code = branching_hmm,
  constants = list(C = 40)
)

sim_model$calculate("f")

nodes_to_sim <- sim_model$getDependencies(
  c("f"),
  self = F,
  downstream = T
)

sim_model$simulate(nodes_to_sim)
