# orderly::orderly_develop_start("explore_qPCR")
# setwd("src/explore_qPCR/")

df <- read_excel("data/2022-05-09_AG_barcodes006_008.xls", sheet = "Results", skip = 42)

df <- df %>%
  #' Split well position into the row and column
  tidyr::separate(`Well Position`, into = c("row", "column"), sep = 1) %>%
  mutate(
    row = as.factor(row),
    column = as.numeric(column)
  ) %>%
  #' Remove columns which are just NA
  select(where(~ any(!is.na(.x)))) %>%
  #' Replace T/F coding with TRUE/FALSE
  mutate(across(CQCONF:OUTLIERRG, ~ as.logical(ifelse(.x == "T", TRUE, FALSE))))

pdf("assay.pdf", h = 5, w = 6.25)

ggplot(df, aes(x = column, y = forcats::fct_rev(row), label = `Sample Name`)) +
  geom_tile(colour = "black", alpha = 0.5) +
  geom_text(size = 2) +
  scale_x_continuous(breaks = 1:9) +
  labs(x = "", y = "") +
  theme_minimal()

dev.off()
