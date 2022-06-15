# orderly::orderly_develop_start("explore_qPCR")
# setwd("src/explore_qPCR/")

cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

files <- list.files(path = "data")

result_extracts <- lapply(files, function(file) {
  df <- readxl::read_excel(
    paste0("data/", file),
    sheet = "Results",
    skip = 42
  )

  #' Sometimes the files are empty
  if(nrow(df) == 0) {
    return(NULL)
  } else {
    return(mutate(df, assay = paste(file)))
  }
})

df <- bind_rows(result_extracts)

df <- df %>%
  #' Split well position into the row and column
  tidyr::separate(`Well Position`, into = c("row", "column"), sep = 1) %>%
  #' Split sample name into barcode and condition
  tidyr::separate(`Sample Name`, into = c("barcode", "condition_id"), sep = "_", remove = FALSE) %>%
  mutate(
    row = as.factor(row),
    column = as.numeric(column),
    condition = as.numeric(stringr::str_extract(condition_id, "\\d+")),
    Ct = ifelse(CT == "Undetermined", NA, as.numeric(CT))
  ) %>%
  #' Remove columns which are just NA
  select(where(~ any(!is.na(.x)))) %>%
  #' Replace T/F coding with TRUE/FALSE
  mutate(across(CQCONF:OUTLIERRG, ~ as.logical(ifelse(.x == "T", TRUE, FALSE))))

#' * `Well` The number of the well
#' * `row` The row that the well is in (A-H)
#' * `column` The column that the well is in (1-9)
#' * `Omit` Should the X be omitted?
#' * `Sample Name`: Format is Q[num1]_D[num2] where
#' [num1] is the tracer / pair identifier
#' [num2] is the experimental condition
#' * `barcode`: the tracer / pair identifier Q[num1]
#' * `condition_id`: the experimental condition D[num2]
#' * `condition`: the experimental condition [num2]
#' Note: NTC refers to "No Template Control"
#' * `Target Name`: These are mostly identical to Sample Name
#' * `Task`: All "UNKNOWN" here (?)
#' * `Reporter`: The type of reporter dye
#' * `Quencher`: The type of quencher
#' Note: for reporter and quencher, see:
#' https://www.ebi.ac.uk/training/online/courses/functional-genomics-ii-common-technologies-and-data-analysis-methods/real-time-pcr/
#' * `Ct`: PCR cycle number at which the reaction curve intersects the threshold (also called Cq)
#' * `Ct Mean`: Mean(Ct) over experiments with the same condition
#' * `Ct SD`: SD(Ct) over experiments with the same condition
#' * `Automatic Ct Threshold`: Is the (automatic) Ct Threshold used?
#' * `Ct Threshold`: The threshold (determined automatically?)
#' * `Automatic Baseline` Is the (automatic) baseline used?
#' * `Baseline Start`: (?)
#' * `Baseline End`: (?) Is diff(baseline) related to Ct?
#' * `Amp Status`: Whether amplification occurred
#' * `Cq Conf`: ?
#' * `Target Color`: ?
#' * `CQCONF`: ?
#' * `THOLDFAIL`: ?
#' * `EXPFAIL`: ?
#' * `HIGHSD`: ?
#' * `SPIKE`: ?
#' * `DRNMIN`: ?
#' * `NOISE`: ?
#' * `OUTLI`: ?

pdf("assay.pdf", h = 5, w = 6.25)

df %>%
  split(.$assay) %>%
  lapply(function(x) {
    ggplot(x, aes(x = column, y = forcats::fct_rev(row), label = `Sample Name`, fill = condition)) +
      geom_tile(alpha = 0.8) +
      geom_text(size = 2) +
      scale_x_continuous(breaks = 1:9) +
      scale_fill_viridis_c() +
      labs(title = paste0("File: ", x$assay[1]), x = "", y = "", fill = "log10 copies + 1 (per uL)") +
      theme_minimal() +
      theme(
        legend.position = "bottom"
      )
  })

dev.off()

pdf("ct.pdf", h = 5, w = 6.25)

df %>%
  filter(!(barcode %in% c("Blank", "NTC"))) %>%
  split(.$barcode) %>%
  lapply(function(x) {
    fit <- lm(Ct ~ condition, data = x)
    coeff <- summary(fit)$coefficients %>% as.data.frame()
    intercept <- coeff$Estimate[1]
    slope <- coeff$Estimate[2]
    E <- -1 + 10^(-1 / slope)

    ggplot(x, aes(x = condition, y = Ct)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = intercept, slope = slope, col = cbpalette[3]) +
      annotate(
        geom = "text", x = max(x$condition) - 2.5, y = max(x$Ct) - 1,
        label = paste0("Equation: y = ", round(slope, 2), "x + ", round(intercept, 2)),
        hjust = 0, vjust = 1, size = 4
      ) +
      annotate(
        geom = "text", x = max(x$condition) - 2.5, y = max(x$Ct) - 2,
        label = paste0("R^2 = ", round(summary(fit)$r.squared, 3)),
        hjust = 0, vjust = 1, size = 4
      ) +
      annotate(
        geom = "text", x = max(x$condition) - 2.5, y = max(x$Ct) - 3,
        label = paste0("Efficiency = ", 100 * round(E, 3), "%"),
        hjust = 0, vjust = 1, size = 4
      ) +
      labs(
        title = paste0("Barcode: ", x$barcode[1]),
        x = "log10 copies + 1 (per uL)"
      ) +
      theme_minimal()
  })

dev.off()

#' Moving to analysis of the amplification curves
#' Question: where is the threshold Rn value?
file <- files[1]

amp <- readxl::read_excel(
  paste0("data/", file),
  sheet = "Amplification Data",
  skip = 42
)

amp <- amp %>%
  select(Well, Cycle, Rn, "Delta Rn") %>%
  left_join(
    df %>%
      filter(assay == file) %>%
      select(Well, row, column, barcode, condition, Ct),
    by = "Well"
  ) %>%
  filter(!is.na(row), !is.na(column), column > 2)

pdf("amp.pdf", h = 8, w = 6.25)

amp %>%
  ggplot(aes(x = Cycle, y = Rn, col = condition)) +
  geom_line(size = 1) +
  geom_vline(
    data = select(amp, row, column, Ct), aes(xintercept = Ct),
    col = "grey40", linetype = "dashed", alpha = 0.8
  ) +
  facet_grid(row ~ column) +
  theme_minimal() +
  scale_color_viridis_c() +
  labs(col = "log10 copies + 1 (per uL)") +
  theme(
    legend.position = "bottom"
  )

dev.off()
