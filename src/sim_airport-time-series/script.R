# orderly::orderly_develop_start("sim_airport-time-series")
# setwd("src/sim_airport-time-series")

set.seed(2)

# Specifying and running the model

# Time-step duration
dt <- 0.05

# Final time point
end <- 50

mod <- stoch_seir_dust$new(
  beta = 2,
  gamma = 1,
  sigma = 1,
  population_size = 10^6,
  start_infections = 10,
  capacity_per_flight = 2000,
  num_flights = 50,
  num_flightsAB = 25,
  dt = dt,
  virus_shed = 10,
  non_virus_shed = 100000,
  volume_wastewater = 500 * 10^6, # 500 litres = 500 * 1000 ml
  shedding_freq = 1,
  bias = 1,
  m_tot = 10000000
)

output <- mod$run(1:(end/dt))
transformed_output <- mod$transform_variables(output)

extract_vars <- function(output, vars) {
  state_df <- do.call(cbind, output[vars]) %>%
    as.data.frame() %>%
    pivot_longer(
      -time,
      names_to = "var",
      values_to = "value"
    )
}

state_df <- extract_vars(
  output = transformed_output,
  vars = c("time", "S", "E", "I", "R")
)

ggplot(state_df, aes(x = time, y = value, colour = var)) +
  geom_line() +
  labs(x = "Time (days)", y = "Number of individuals in each state", col = "State") +
  theme_minimal()

# Number of infections on each flight per day
num_flightsAB <- mod$contents()[["num_flightsAB"]]
indiv_flight_infections <- data.frame(time = transformed_output$time, transformed_output$n_inf_specific_flightABOut)
colnames(indiv_flight_infections) <- c("time", paste0("flight", 1:num_flights))

indiv_flight_infections

airplane_infections_df <- indiv_flight_infections %>%
  pivot_longer(
    -time,
    names_to = "flight",
    values_to = "infections",
    names_prefix = "flight",
    names_transform = as.numeric
  ) %>%
  mutate(
    day_interval = cut(time, breaks = max(time)),
    day = cut_limit(day_interval, type = "upper")
  ) %>%
  group_by(flight, day) %>%
  summarise(infections = sum(infections))

final_day <- 25

pdf("airplane-infections.pdf", h = 5, w = 6.25)

ggplot(airplane_infections_df, aes(x = day, y = infections, group = day)) +
  geom_boxplot(fill = NA, outlier.colour = NA, alpha = 0.1) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.5) +
  theme(legend.position = "none") +
  lims(x = c(0, final_day), y = c(0, 20)) +
  labs(x = "Time (days)", y = "Number of infections on each flight") +
  theme_minimal()

ggplot(airplane_infections_df, aes(x = day, y = infections, group = flight)) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  lims(x = c(0, final_day), y = c(0, 20)) +
  labs(x = "Time (days)", y = "Number of infections on each flight") +
  theme_minimal()

dev.off()

# Number of mapped reads on each flight per day (deterministic model)
indiv_flight_reads <- data.frame(time = transformed_output$time, transformed_output$m_i_Out)
colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flights))

airplane_reads_df <- indiv_flight_reads %>%
  pivot_longer(
    -time,
    names_to = "flight",
    values_to = "reads",
    names_prefix = "flight",
    names_transform = as.numeric
  ) %>%
  mutate(
    day_interval = cut(time, breaks = max(time)),
    day = cut_limit(day_interval, type = "upper")
  ) %>%
  group_by(flight, day) %>%
  summarise(reads = sum(reads))

pdf("airplane-reads-deterministic.pdf", h = 5, w = 6.25)

ggplot(airplane_reads_df, aes(x = day, y = reads, group = day)) +
  geom_boxplot(fill = NA, outlier.colour = NA, alpha = 0.1) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.5) +
  theme(legend.position = "none") +
  lims(x = c(0, final_day), y = c(0, 12.5)) +
  labs(x = "Time (days)", y = "Number of reads from each flight's wastewater") +
  theme_minimal()

ggplot(airplane_reads_df, aes(x = day, y = reads, group = flight)) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  lims(x = c(0, final_day), y = c(0, 12.5)) +
  labs(x = "Time (days)", y = "Number of reads from each flight's wastewater") +
  theme_minimal()

dev.off()

# Number of mapped reads on each flight per day (stochastic model)
indiv_flight_reads <- data.frame(time = transformed_output$time, transformed_output$stoch_m_i_pois_Out)
colnames(indiv_flight_reads) <- c("time", paste0("flight", 1:num_flights))

airplane_reads_df <- indiv_flight_reads %>%
  pivot_longer(
    -time,
    names_to = "flight",
    values_to = "reads",
    names_prefix = "flight",
    names_transform = as.numeric
  ) %>%
  mutate(
    day_interval = cut(time, breaks = max(time)),
    day = cut_limit(day_interval, type = "upper")
  ) %>%
  group_by(flight, day) %>%
  summarise(reads = sum(reads))

pdf("airplane-reads-stochastic.pdf", h = 5, w = 6.25)

ggplot(airplane_reads_df, aes(x = day, y = reads, group = day)) +
  geom_boxplot(fill = NA, outlier.colour = NA, alpha = 0.1) +
  geom_jitter(size = 0.5, width = 0.15, alpha = 0.5) +
  theme(legend.position = "none") +
  lims(x = c(0, final_day), y = c(0, 12.5)) +
  labs(x = "Time (days)", y = "Number of reads from each flight's wastewater") +
  theme_minimal()

ggplot(airplane_reads_df, aes(x = day, y = reads, group = flight)) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  lims(x = c(0, final_day), y = c(0, 12.5)) +
  labs(x = "Time (days)", y = "Number of reads from each flight's wastewater") +
  theme_minimal()

dev.off()
