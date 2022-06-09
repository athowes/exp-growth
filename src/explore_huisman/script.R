#' https://github.com/JSHuisman/wastewaterRe
#' https://ehp.niehs.nih.gov/doi/10.1289/EHP10050

#' Relation of wastewater to incidence or transmision rates is driven by assumptions about
#' virus excretion rates into sewer system.
#' Excretion varies by individual, and through time after infection, summarised by the
#' shedding load distribution (SLD).
#'
#' Pepper mild mottle virus (PMMoV) is a plant virus found in wastewater in high concentrations
#' and in fairly constant loads and serves to detect anomalies in the collected sample.
#'
#' Concentrations of RNA targets multiplied by daily flow rate to estimate the total number of
#' genome copies (gc) shed by people within the catchment per day.
#' Samples with PMMoV outside the mean plus or minus three standard deviations were discarded.

#' i, j: day index
#' w_j: the SLD such that \sum_j w_j = 1 (i.e. it's a distribution)
#' C_i: viral RNA on day i
#' I_i: infections on day i
#' M: normalisation factor, accounting for details of the sewer system, sampling point
#' within wastewater treatment plant, choice of sample matrix and processing pipeline
#' (to what extent does this need to be determined if modelling multiple systems together?)
#' N: total amount of virus shed by infected ind

#' Used interpolation and smoothing to fill gaps in wastewater data
#' Used EM algorithm for deconvolving the time series

#' Gastrointestinal SLD from Benefield was Gamma with mean 6.7 and SD 7.0
gamma_mean <- 6.7
gamma_sd <- 7.0

#' gamma_mean = alpha / beta
#' gamma_variance = alpha / beta^2
beta <- gamma_mean / gamma_sd^2 #'rate
alpha <- gamma_mean * beta #' shape

pdf("gastro-sld-benefield.pdf", h = 5, w = 6.25)

data.frame(x = 1:14, y = dgamma(1:14, shape = alpha, rate = beta)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  labs(x = "Day", y = "", title = "Gastrointestinal SLD from Benefield ") +
  theme_minimal()

dev.off()

#' Tried subsampling to see how important it has to have high frequency data
