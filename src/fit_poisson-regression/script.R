# orderly::orderly_develop_start("fit_poisson-regression")
# setwd("src/fit_poisson-regression")

cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Import simulated data
df <- readRDS("depends/sample0.rds") %>%
  as_tibble()

#' Check what the data looks like
pdf("data-sample0.pdf", h = 4, w = 6.25)

df %>%
  ggplot(aes(x = day, y = count, col = regime, group = id)) +
  geom_line(alpha = 0.1) +
  theme_minimal() +
  scale_color_manual(values = cbpalette) +
  labs(x = "Day", y = "Count", col = "Regime")

dev.off()

#' Fit to each unique id
#' Note that kmer is not unique, and represents the kmer within a given organism
fits <- tibble(id = unique(df$id)) %>%
  mutate(
    fit = map(id, ~glm(count ~ 1 + day, family = "poisson", data = filter(df, id == .))),
    fit = map(fit, broom::tidy, conf.int = TRUE)
  ) %>%
  unnest(fit)

pdf("fits-sample0.pdf", h = 4, w = 6.25)

fits %>%
  left_join(
    select(df, id, regime),
    by = "id"
  ) %>%
  filter(term == "day") %>%
  ggplot(aes(x = id, y = estimate, col = regime)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    scale_color_manual(values = cbpalette) +
    labs(x = "k-mer ID", y = "Slope point estimate", col = "Regime")

dev.off()

saveRDS(fits, "fits-sample0.rds")

bm <- bench::mark(
  glm(count ~ 1 + day, family = "poisson", data = filter(df, id == 1))
)

saveRDS(bm, "benchmark-sample0.rds")
