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

df <- bind_rows(result_extracts) %>%
  janitor::clean_names()

df <- df %>%
  #' Split well position into the row and column
  tidyr::separate(well_position, into = c("row", "column"), sep = 1) %>%
  #' Split sample name into barcode and condition
  tidyr::separate(sample_name, into = c("barcode", "condition_id"), sep = "_", remove = FALSE) %>%
  mutate(
    row = as.factor(row),
    column = as.numeric(column),
    condition = as.numeric(stringr::str_extract(condition_id, "\\d+")),
    ct = ifelse(ct == "Undetermined", NA, as.numeric(ct))
  ) %>%
  #' Remove columns which are just NA
  select(where(~ any(!is.na(.x)))) %>%
  #' Replace T/F coding with TRUE/FALSE
  mutate(across(cqconf:outlierrg, ~ as.logical(ifelse(.x == "T", TRUE, FALSE))))

#' * `well` The number of the well
#' * `row` The row that the well is in (A-H)
#' * `column` The column that the well is in (1-9)
#' * `mit` Should the X be omitted?
#' * `sample_name`: Format is Q[num1]_D[num2] where
#' [num1] is the tracer / pair identifier
#' [num2] is the experimental condition
#' * `barcode`: the tracer / pair identifier Q[num1]
#' * `condition_id`: the experimental condition D[num2]
#' * `condition`: the experimental condition [num2]
#' Note: NTC refers to "No Template Control"
#' * `target_name`: These are mostly identical to Sample Name
#' * `task`: All "UNKNOWN" here (?)
#' * `reporter`: The type of reporter dye
#' * `quencher`: The type of quencher
#' Note: for reporter and quencher, see:
#' https://www.ebi.ac.uk/training/online/courses/functional-genomics-ii-common-technologies-and-data-analysis-methods/real-time-pcr/
#' * `ct`: PCR cycle number at which the reaction curve intersects the threshold (also called Cq)
#' * `ct_mean`: Mean(Ct) over experiments with the same condition
#' * `ct_sd`: SD(Ct) over experiments with the same condition
#' * `automatic_ct_threshold`: Is the (automatic) Ct Threshold used?
#' * `ct_threshold`: The threshold (determined automatically)
#' * `automatic_baseline` Is the (automatic) baseline used?
#' * `baseline_start`: (?)
#' * `baseline_end`: (?) Is diff(baseline) related to Ct?
#' * `amp_status`: Whether amplification occurred
#' * `cq_conf`: ?
#' * `target_color`: ?
#' * `cqconf`: ?
#' * `tholdfail`: ?
#' * `expfail`: ?
#' * `highsd`: ?
#' * `spike`: ?
#' * `drnmin`: ?
#' * `noise`: ?
#' * `outlierrg`: ?
#' * `assay`: ?
#' * `prfdrop`: ?
#' * `badrox`: ?
#' * `prflow`: ?
#' * `ctfail`: ?
#' * `noamp`: ?

pdf("assay.pdf", h = 5, w = 6.25)

df %>%
  split(.$assay) %>%
  lapply(function(x) {
    ggplot(x, aes(x = column, y = forcats::fct_rev(row), label = sample_name, fill = condition)) +
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
    fit <- lm(ct ~ condition, data = x)
    coef <- broom::tidy(fit)
    intercept <- coef$estimate[1]
    slope <- coef$estimate[2]
    E <- -1 + 10^(-1 / slope)

    ggplot(x, aes(x = condition, y = ct)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = intercept, slope = slope, col = cbpalette[3]) +
      annotate(
        geom = "text", x = max(x$condition) - 2.5, y = max(x$ct) - 1,
        label = paste0("Equation: y = ", round(slope, 2), "x + ", round(intercept, 2)),
        hjust = 0, vjust = 1, size = 4
      ) +
      annotate(
        geom = "text", x = max(x$condition) - 2.5, y = max(x$ct) - 2,
        label = paste0("R^2 = ", round(summary(fit)$r.squared, 3)),
        hjust = 0, vjust = 1, size = 4
      ) +
      annotate(
        geom = "text", x = max(x$condition) - 2.5, y = max(x$ct) - 3,
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
amp_extracts <- lapply(files, function(file) {
  amp <- readxl::read_excel(
    paste0("data/", file),
    sheet = "Amplification Data",
    skip = 42
  )

  #' Sometimes the files are empty
  if(nrow(amp) == 0) {
    return(NULL)
  } else {
    return(mutate(amp, assay = paste(file)))
  }
})

amp <- bind_rows(amp_extracts) %>%
  janitor::clean_names()

amp <- amp %>%
  select(well, cycle, rn, delta_rn, assay) %>%
  left_join(
    df %>%
      select(well, row, column, barcode, condition, ct, ct_threshold, assay),
    by = c("well", "assay")
  ) %>%
  filter(!is.na(row), !is.na(column), column > 2)

pdf("amp.pdf", h = 8, w = 6.25)

amp %>%
  split(.$assay) %>%
  lapply(function(x) {
    ggplot(x, aes(x = cycle, y = delta_rn, col = condition)) +
    geom_line(size = 1) +
    geom_vline(
      data = select(x, row, column, ct), aes(xintercept = ct),
      col = "grey40", linetype = "dashed", alpha = 0.8
    ) +
    geom_hline(
      data = select(x, row, column, ct_threshold), aes(yintercept = ct_threshold),
      col = "grey40", linetype = "dashed", alpha = 0.8
    ) +
    facet_grid(row ~ column) +
    theme_minimal() +
    scale_color_viridis_c() +
    labs(col = "log10 copies + 1 (per uL)") +
    theme(
      legend.position = "bottom"
    )
  })

dev.off()

pdf("log-amp.pdf", h = 8, w = 6.25)

amp %>%
  split(.$assay) %>%
  lapply(function(x) {
    ggplot(x, aes(x = cycle, y = log(rn), col = condition)) +
      geom_line(size = 1) +
      geom_vline(
        data = select(x, row, column, ct), aes(xintercept = ct),
        col = "grey40", linetype = "dashed", alpha = 0.8
      ) +
      facet_grid(row ~ column) +
      theme_minimal() +
      scale_color_viridis_c() +
      labs(col = "log10 copies + 1 (per uL)") +
      theme(
        legend.position = "bottom"
      )
  })

dev.off()

#' Trying to determine the difference between rn and delta_rn
amp <- amp %>%
  mutate(rn_diff = rn - delta_rn)

pdf("rn-diff.pdf", h = 5, w = 6.25)

amp %>%
  split(.$assay) %>%
  lapply(function(x) {
    ggplot(x, aes(x = cycle, y = rn_diff)) +
    geom_line() +
    facet_grid(row ~ column) +
    theme_minimal()
  })

dev.off()

baselines <- df %>%
  filter(
    well == 3,
    assay == "2022-05-09_AG_barcodes006_008.xls"
  ) %>%
  select(baseline_start, baseline_end)

rn_subset <- amp %>%
  filter(
    well == 3,
    assay == "2022-05-09_AG_barcodes006_008.xls",
    cycle %in% baselines$baseline_start:baselines$baseline_end
  )

rn_fit <- lm(rn ~ 1 + cycle, data = rn_subset)
rn_coef <- broom::tidy(rn_fit)

pdf("rn-diff-test.pdf", h = 5, w = 6.25)

amp %>%
  filter(
    well == 3,
    assay == "2022-05-09_AG_barcodes006_008.xls"
  ) %>%
  ggplot(aes(x = cycle, y = rn_diff)) +
    geom_line() +
    geom_abline(
      intercept = rn_coef$estimate[1],
      slope = rn_coef$estimate[2],
      col = "red",
      linetype = "dashed"
    ) +
  theme_minimal()

dev.off()
