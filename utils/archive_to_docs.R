#' Set the working directory to the project root
setwd(rprojroot::find_rstudio_root_file())

archive_to_docs <- function(report) {
  #' Artefacts to be moved
  filenames <- yaml::read_yaml(file.path(paste0("src/", report, "/orderly.yml")))$artefacts[[1]]$data$filenames

  #' Latest version in archive
  latest <- orderly::orderly_latest(report)

  #' Copy files over
  files_from <- paste0("archive/", report, "/", latest, "/", filenames)
  files_to <- paste0("docs/", filenames)

  file.copy(from = files_from, to = files_to, overwrite = TRUE)
}

#' Names of the reports to move
archive_to_docs("explore_poisson-regression")
archive_to_docs("docs_CoDA")
archive_to_docs("fit_qPCR")
