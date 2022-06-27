#' Model code based on suggestion by Perry de Valpine
library(nimble)

branching_hmm <- nimbleCode({
  N_cont ~ dunif(0, 1000)
  N <- round(N_cont)
  r ~ dunif(0, 1)
  E ~ dunif(0, 1)
  sigma_f ~ T(dnorm(0, sd = 1), 0, Inf)
  alpha ~ dunif(0, 100)
  n[1] ~ dbinom(prob = r, size = N)
  for(i in 2:C) {
    n_new[i] ~ dbinom(prob = E, size = n[i - 1])
    n[i] <- n[i - 1] + n_new[i]
    f[i] ~ dnorm(alpha * n[i], sd = sigma_f)
  }
})

# inits <- list(N_cont = 500, r = 0.9, E = 0.9, sigma_f = 0.1, alpha = 10)

sim_model <- nimbleModel(
  code = branching_hmm,
  constants = list(C = 20)
)

sim_model$initializeInfo()
#' [Note] Missing values (NAs) or non-finite values were found in model variables: N_cont, N, r, E, sigma_f, alpha, n, n_new, lifted_alpha_times_n_oBi_cB_L11, f.
#' [Note] This is not an error, but some or all variables may need to be initialized for certain algorithms to operate properly.
#' [Note] For more information on model initialization, see help(modelInitialization).

#' "MCMC will auto-initialize but will do so from the prior distribution.
#' This can cause slow convergence, especially in the case of diffuse priors."
#' This is exactly what we want for prior sampling, so should leave the model without inits above

#' The variable names are
sim_model$getVarNames()
#' "N_cont", "N", "r", "E", "sigma_f", "alpha", "n", "n_new", "lifted_alpha_times_n_oBi_cB_L11", "f"
#' The "lifted" variable is to do with the multiplication between alpha and n in line 11

#' The node names are
sim_model$getNodeNames()
#' Similar to the variable names, but including all of the indices

sim_model$getNodeNames(determOnly = TRUE)
sim_model$getNodeNames(stochOnly = TRUE)

#' There is no data currently (I just want to simulate from the model)
sim_model$getNodeNames(dataOnly = TRUE)

#' I think these are the nodes which are children of this node (up to the next stochastic node)
#' You can get not just the direct children by using downstream = TRUE to give all descendants
#' Parent nodes can be determined with $getParents()
sim_model$getDependencies("N_cont")
sim_model$getDependencies("N")
sim_model$getDependencies("r")
sim_model$getDependencies("E")
sim_model$getDependencies("sigma_f")
sim_model$getDependencies("alpha")

#' iGraph
plot(sim_model$getGraph())

#' The fluorescence values are end nodes
sim_model$isEndNode("f")

sim_model$getDistribution("N_cont")
sim_model$getDistribution("N") #' For discrete parameters this should just give NA I suppose
sim_model$getDistribution("r")

#' Simulate all variables from the prior
all_vars <- sim_model$getVarNames()
sim_model$simulate(all_vars)

#' Looks like exponential growth to me!
plot(sim_model$n)
plot(log(sim_model$n))
plot(sim_model$f)
plot(log(sim_model$f))
