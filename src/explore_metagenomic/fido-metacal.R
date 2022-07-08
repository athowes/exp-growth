#' Following https://microbiomemeasurement.org/posts/bias-estimation-from-mocks-with-fido/

#' Advantages over metacal:
#' * Model log relative efficiencies as linear functions
#' * Random counting error inherent to sequencing
#' * Regularise efficiency estimates using a prior (or model relatedness between taxa)

library(tidyverse)
library(metacal)
library(fido)
library(ggdist)
library(cowplot)
library(patchwork)

sam <- system.file(
  "extdata", "brooks2015-sample-data.csv", package = "metacal"
) %>%
  read_csv(col_types = "cffcic") %>%
  dplyr::rename_with(str_to_lower)

observed <- system.file(
  "extdata", "brooks2015-observed.csv", package = "metacal"
) %>%
  read_csv %>%
  select(-Other) %>%
  column_to_rownames("Sample") %>%
  as("matrix")

actual <- system.file(
  "extdata", "brooks2015-actual.csv", package = "metacal"
) %>%
  read_csv %>%
  column_to_rownames("Sample") %>%
  as("matrix")

stopifnot(setequal(colnames(actual), colnames(observed)))
stopifnot(setequal(rownames(actual), rownames(observed)))
stopifnot(setequal(rownames(actual), sam$sample))

observed <- observed %>% t
observed <- observed[, sam$sample]
actual <- actual %>% t
actual <- actual[rownames(observed), colnames(observed)]

stopifnot(identical(dimnames(observed), dimnames(actual)))
stopifnot(identical(colnames(observed), sam$sample))

sam %>% glimpse()
observed %>% corner()
actual %>% corner() %>% round(2)

n_samples <- ncol(observed)
n_taxa <- nrow(observed)
