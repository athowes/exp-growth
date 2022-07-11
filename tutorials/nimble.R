#' Following https://r-nimble.org/html_manual/cha-lightning-intro.html
#' Other resources which may be useful:
#' * https://r-nimble.org/cheatsheets/NimbleCheatSheet.pdf
#' * https://groups.google.com/g/nimble-users

library(nimble)
library(tidyverse)

pumpCode <- nimbleCode({
  for (i in 1:N){
    theta[i] ~ dgamma(alpha,beta)
    lambda[i] <- theta[i]*t[i]
    x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

pumpConsts <- list(
  N = 10,
  t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
)

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1, theta = rep(0.1, pumpConsts$N))

#' x[i] is the number of failures recorded during t[i] for pump i
#' theta[i] is a failure rate
#' Want to estimate alpha and beta

pump <- nimbleModel(
  code = pumpCode,
  name = "pump",
  constants = pumpConsts,
  data = pumpData,
  inits = pumpInits
)

pump$getNodeNames()

pump$x
pump$logProb_x
pump$alpha
pump$theta

pump$plotGraph()

pump$getDependencies(c("alpha", "beta"))
pump$getDependencies(c("alpha", "beta"), determOnly = TRUE)

pump[["lifted_d1_over_beta"]]

set.seed(1)

pump$simulate("theta")
pump$theta

pump$lambda
pump$logProb_x
pump$getLogProb("x")

pump$calculate(pump$getDependencies("theta"))

pump$lambda

Cpump <- compileNimble(pump)

Cpump$theta

#' N.B. variables included in "monitors" have their posterior sampples recorded
#' The default variables to be monitored are all top-level stochastic nodes
mcmc_out <- nimbleMCMC(
  code = pumpCode,
  constants = pumpConsts,
  data = pumpData,
  inits = pumpInits,
  nchains = 2,
  niter = 1000,
  summary = TRUE,
  WAIC = TRUE,
  monitors = c("alpha", "beta", "theta")
)

names(mcmc_out)

mcmc_out$summary
mcmc_out$WAIC

pumpConf <- configureMCMC(pump, print = TRUE)

#' Conjugate relationships are detected for all nodes except alpha
pumpConf$addMonitors(c("alpha", "beta", "theta"))

pumpMCMC <- buildMCMC(pumpConf)
CpumpMCMC <- compileNimble(pumpMCMC, project = pump)

niter <- 1000
set.seed(1)
samples <- runMCMC(CpumpMCMC, niter = niter)

plot(samples[, "alpha"], type = "l")
plot(samples[, "beta"], type = "l")
plot(samples[, "alpha"], samples[, "beta"]) #' You can see that the alpha and beta samples are correlated

acf(samples[, "alpha"]) #' Correlation in the chains up to lag 6 or 7
acf(samples[, "beta"])

#' Customising the MCMC
#' NIMBLE has lots of options for customising the MCMC algorithm that is used in order to improve sampling performance
#' Here are some resources I've quickly found about it:
#' * https://www.stat.berkeley.edu/~paciorek/presentations/paciorek-isba16.pdf
#' * https://mmeredith.net/blog/2021/NIMBLE_adaptation.htm
#' Going to add an adaptive block sampler on alpha and beta jointly
pumpConf$addSampler(
  target = c("alpha", "beta"),
  type = "RW_block",
  control = list(adaptInterval = 100)
)

#' I get this note:
#' [Note] Assigning an RW_block sampler to nodes with very different scales can result in low MCMC efficiency.
#' If all nodes assigned to RW_block are not on a similar scale, we recommend providing an informed value for
#' the "propCov" control list argument, or using the AFSS sampler instead.
#' I can't remember details on this but

pumpMCMC2 <- buildMCMC(pumpConf)

# need to reset the nimbleFunctions in order to add the new MCMC
CpumpNewMCMC <- compileNimble(pumpMCMC2, project  = pump,
                              resetFunctions = TRUE)

CpumpNewMCMC$run(niter)

samples_new <- as.matrix(CpumpNewMCMC$mvSamples)

plot(samples_new[ , "alpha"], type = "l", xlab = "iteration", ylab = expression(alpha))
plot(samples_new[ , "beta"], type = "l", xlab = "iteration", ylab = expression(beta))
plot(samples_new[ , "alpha"], samples_new[ , "beta"], xlab = expression(alpha), ylab = expression(beta))
plot(samples_new[ , "theta[1]"], type = "l", xlab = "iteration", ylab = expression(theta[1]))

#' Not convinced that this looks substantially better, but ok, the point is that you can customise the sampler

#' You can also run other algorithms, like Monte Carlo EM
pump2 <- pump$newModel()
box <- list(list(c("alpha", "beta"), c(0, Inf)))
pump_mcem <- buildMCEM(model = pump2, latentNodes = "theta[1:10]", boxConstraints = box)
pump_mle <- pump_mcem$run()
pump_mle

#' You can also create your own functions in NIMBLE

sim_nodes_many <- nimbleFunction(
  setup <- function(model, nodes) {
    mv <- modelValues(model)
    deps <- model$getDependencies(nodes)
    all_nodes <- model$getNodeNames()
  },
  run <- function(n = integer()) {
    resize(mv, n)
    for(i in 1:n) {
      model$simulate(nodes)
      model$calculate(deps)
      copy(from = model, nodes = all_nodes, to = mv, rowTo = i, logProb = TRUE)
    }
  }
)

sim_nodes_theta1to5 <- sim_nodes_many(pump, "theta[1:5]")
sim_nodes_theta6to10 <- sim_nodes_many(pump, "theta[6:10]")

#' setup is written in R
#' run is written in NIMBLE
#' run requires type information e.g. n = integer()
#' for loop looks like R, but only sequential integer iteration is allowed
#' calculate and simulate can be used in NIMBLE
#' copy is a special function which can be used to record values from the modelinto modelValues object
#' multiple instances / specialisations can be made by calling sim_nodes_many with different arguments
