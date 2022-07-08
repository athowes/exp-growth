#' Following https://jsilve24.github.io/fido/articles/introduction-to-fido.html

#' Bayesian multinomial logistic-normal models
#' Pibble is one type of fido model (for multivariate linear regression)
#'
#' * Y is a [D x N] matrix of counts
#' * j-th column of Y as Y_j
#' * X is a [Q x N] covariate matrix
#'
#' Then the pibble model is
#'
#' Y_j ~ Multinomial(pi_j)
#' pi_j = phi^-1(eta_j)
#' eta_j ~ N(Lambda X_j, Sigma)
#' Sigma ~ MN(Theta, Sigma, Gamma)
#' Sigma ~ W^-1()
#'
#' where MN is a matrix normal, W^-1 is the inverse Wishart, phi(pi_j) = eta_j is some
#' transformation from reals to the simplex e.g. ALR (a.k.a. identified softmax)
#' Note that the particular transform used is not important.
#'
#' Justin highlights that the main modelling aspect of pibble is the MN distribution.

#' Let's analyze some data!
library(MicrobeDS)
library(phyloseq)
library(dplyr)
library(fido)

set.seed(899)

data("RISK_CCFA")

#' Drop low abundant taxa and samples
dat <- RISK_CCFA %>%
  subset_samples(
    disease_stat != "missing",
    immunosup != "missing"
  ) %>%
  subset_samples(diagnosis %in% c("no", "CD")) %>%
  subset_samples(steroids=="false") %>%
  subset_samples(antibiotics=="false") %>%
  subset_samples(biologics=="false") %>%
  subset_samples(biopsy_location=="Terminal ileum") %>%
  tax_glom("Family") %>%
  prune_samples(sample_sums(.) >= 5000, .) %>%
  filter_taxa(function(x) sum(x > 3) > 0.10*length(x), TRUE)

sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>%
  mutate(
    age = as.numeric(as.character(age)),
    diagnosis = relevel(factor(diagnosis, ordered = FALSE), ref = "no"),
    disease_stat = relevel(factor(disease_stat, ordered = FALSE), ref = "non-inflamed")
  )

X <- t(model.matrix(~diagnosis + disease_stat + age, data = sample_dat))
Y <- otu_table(dat)
useful::corner(X)
useful::corner(Y)

#' Now we're going to specify priors
#' Sigma is the covariance between log-ratios
#' We will get to Sigma via Omega, the covariance on log-absolute abundances
#'
#' Sigma = G Omega G^T
#'
#' where G is a [D - 1 x D] matrix given by G = [I_{D - 1}; -1_{D - 1}] -- i.e the ALR_D
#' contrast matrix.

D <- ntaxa(dat)
upsilon <- D + 3 #' D + 3
Omega <- diag(D)
G <- cbind(diag(D - 1), -1)
Xi <- (upsilon - D) * G %*% Omega %*% t(G)

Theta <- matrix(0, D - 1, nrow(X))
Gamma <- diag(nrow(X))

#' Prior predictive checks recommended!
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
print(priors)
ppc(priors)

#' fido uses ALR coordinate system by default
#' But it designed to work with many different coordinate systems (ALR, CLR, ILR, proportions)

priors <- to_clr(priors)
summary(priors, pars = "Lambda", gather_prob = TRUE, as_factor = TRUE, use_names = TRUE)

names_covariates(priors) <- rownames(X)
p <- plot(priors, par = "Lambda") +
  ggplot2::xlim(c(-10, 10))

priors$Y <- Y
posterior <- refit(priors, optim_method = "adam")

tax <- tax_table(dat)[, c("Class", "Family")]
tax <- apply(tax, 1, paste, collapse = "_")
names_categories(posterior) <- tax

ppc(posterior) + ggplot2::coord_cartesian(ylim = c(0, 30000))
ppc_summary(posterior)

posterior_summary <- summary(posterior, pars = "Lambda")$Lambda

#' Ones where the 95% CI doesn't contain zero
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5), ]
focus <- unique(focus$coord)
plot(posterior, par = "Lambda", focus.coord = focus, focus.cov = rownames(X)[2:4])

#' age has non-zero effect, but it is very small

posterior_summary <- filter(posterior_summary, covariate == "diagnosisCD")
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)

tax_table(dat)[taxa_names(dat)[which(names_coords(posterior) %in% focus)]]

plot(posterior, par = "Lambda", focus.coord = focus, focus.cov = rownames(X)[2])
