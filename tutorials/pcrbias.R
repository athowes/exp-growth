#' Notes on https://www.biorxiv.org/content/10.1101/604025v1.full.pdf
#' Following https://jsilve24.github.io/fido/articles/mitigating-pcrbias.html

#' Preferential amplificaiton of over 3.5x for specific templates
#' Mismatch between primer and template causes up to 10x error
#'
#' Computational method called "alpine" requies mocks
#' They're claiming that they have a calibration curve which doesn't
#'
#' w_ij = a_j b_j ^ x_i
#' log(w_i) = log(a) = log(b) x_i

library(fido)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(5903)
data(pcrbias_mock)

useful::corner(Y)
head(metadata)

X <- t(model.matrix(~ cycle_num + sample_num + machine - 1, data = metadata))
X[, 1:5]

fit <- pibble(Y = Y, X = X, Gamma = 10 * diag(nrow(X)))
fit <- to_clr(fit)

focus_covariate <- rownames(X)[which(grepl("sample_num", rownames(X)))]

plot(fit, par="Lambda", focus.cov = focus_covariate) +
  theme(strip.text.y=element_text(angle=0, hjust=1)) +
  facet_grid(.data$covariate ~ .)
