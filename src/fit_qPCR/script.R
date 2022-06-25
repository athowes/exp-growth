# orderly::orderly_develop_start("fit_qPCR")
# setwd("src/fit_qPCR")

#' Run models
source("4pl/script.R")
source("4pl-processing/script.R")

#' Generate report
rmarkdown::render("qPCR.Rmd")
