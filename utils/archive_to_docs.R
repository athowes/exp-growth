#' Set the working directory to the project root
setwd(rprojroot::find_rstudio_root_file())

#' Name of the report to move
report <- "explore_poisson-regression"

#' Artefacts to be moved
filenames <- yaml::read_yaml(file.path(paste0("src/", report, "/orderly.yml")))$artefacts[[1]]$data$filenames

#' Latest version in archive
latest <- orderly::orderly_latest(report)

#' Copy files over
files_from <- paste0("archive/", report, "/20220711-121610-2153998e/", filenames)
files_to <- paste0("docs/", filenames)

file.copy(from = files_from, to = files_to, overwrite = TRUE)
