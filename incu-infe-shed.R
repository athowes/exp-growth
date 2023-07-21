library(readr)
library(tidyverse)

cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

url <- "https://raw.githubusercontent.com/lennijusten/Epi-curves/master/pathogens_of_concern.csv?token=GHSAT0AAAAAABUU3LCZOUJFPJUD56V4524QYYRGGZA"
poc <- read_csv(url)

names(poc)

df <- poc %>%
  select(pathogen, incubation_period_mean, infectious_period_duration_mean, infectious_before_symptoms_mean, shedding_duration_mean, shedding_before_symptoms_mean) %>%
  mutate(
    incubation_start = ifelse(!is.na(incubation_period_mean), 0, NA),
    incubation_end = incubation_start + incubation_period_mean,
    infectious_start = incubation_end + infectious_before_symptoms_mean,
    infectious_end = infectious_start + infectious_period_duration_mean,
    shedding_start = incubation_end + shedding_before_symptoms_mean,
    shedding_end = shedding_start + shedding_duration_mean
  ) %>%
  select(pathogen, ends_with("start"), ends_with("end")) %>%
  pivot_longer(
    cols = -pathogen,
    names_to = "x",
    values_to = "time"
  ) %>%
  separate(x, c("status", "starting"))

pdf("incu-infe-shed.pdf", h = 7, w = 6.25)

df %>%
  mutate(
    status = str_to_title(status),
    starting = str_to_title(starting)
  ) %>%
  filter(!(pathogen %in% c("Tuberculosis", "HIV"))) %>%
  mutate(pathogen_status = interaction(pathogen, status)) %>%
  ggplot(aes(x = fct_rev(pathogen), y = time, group = status, col = status, shape = starting)) +
  geom_point(position = position_dodge(width = 0.6), size = 2, alpha = 0.8) +
  geom_line(aes(group = pathogen_status), position = position_dodge(width = 0.6), alpha = 0.8) +
  scale_shape_manual(values = c(15, 19)) +
  scale_colour_manual(values = cbpalette) +
  coord_flip() +
  labs(x = "Pathogen", y = "Day", col = "", shape = "") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
  )

dev.off()
