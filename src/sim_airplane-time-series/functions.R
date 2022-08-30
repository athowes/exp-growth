stoch_seir_dust <- odin::odin({

  # Epidemiological module ----

  # Probability of a contact successfully transmitting the disease
  beta <- user()

  # Rate of transition from E to I (incubation period)
  gamma <- user()

  # Rate of transition from I to R (rate of recovery)
  sigma <- user()

  # Overall size of population
  population_size <- user()

  # Starting number of infections (in Exposed compartment)
  start_infections <- user()

  # Force of infection parameter
  lambda <- ((beta * I) / N)

  # Convert from epidemiological rates to probabilities of leaving each compartment
  p_SE <- 1 - exp(-lambda * dt) # S to E
  p_EI <- 1 - exp(-gamma * dt)  # E to I
  p_IR <- 1 - exp(-sigma * dt)  # I to R

  # Binomial draws for number of people moving between compartments at each timestep
  n_SE <- rbinom(S, p_SE) # S to E
  n_EI <- rbinom(E, p_EI) # E to I
  n_IR <- rbinom(I, p_IR) # I to R

  # Stochastic model updates for epidemiological states (S, E, I, R) and quantities of interest (new infectious)
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR - n_inf_flight
  update(R) <- R + n_IR # TODO: return n_inf_flight to here
  update(N) <- S + E + I + R
  update(new_infectious) <- n_EI

  # Note: odin.dust doesn't let you output calculated quantities like n_EI & n_IR without making
  # them tracked variables, which requires having initial() and update() calls for them
  update(n_EI_Output) <- n_EI
  update(n_IR_Output) <- n_IR

  # Initial values for epidemiological states and quantities of interest
  initial(S) <- population_size - start_infections
  initial(E) <- start_infections
  initial(I) <- 0
  initial(R) <- 0
  initial(N) <- S + E + I + R
  initial(new_infectious) <- 0
  initial(n_EI_Output) <- 0
  initial(n_IR_Output) <- 0

  # Surveillance module ----

  # Capacity of a single flight
  capacity_per_flight <- user()

  # Number of flights per day
  num_flights <- user()

  # Number of flights per day from Location A to Location B (NAO Location)
  num_flightsAB <- user()

  # Calculate the number of infected people taking a flight based on number of infections and airport / airplane parameters
  # Note: the rhyper function has inputs and outputs as follows
  # * 1st arg = # white balls
  # * 2nd arg = # black balls
  # * 3rd arg = # balls drawn
  # * output = # white balls that get drawn in this context
  n_inf_flight <- rhyper(I, population_size - I, dt * num_flights * capacity_per_flight)

  # Note: still figuring out whether dt is required here: 2nd arg divided by 3rd arg is ~equivalent to p for binomial,
  # so it doesn't make any huge practical difference (especially when 1st arg << 2nd & 3rd args)
  n_inf_flightAB <- rhyper(n_inf_flight, num_flights * capacity_per_flight * dt, num_flightsAB * capacity_per_flight * dt)

  # Distributing (detectable) infections that take a flight (n_inf_take_flight) across all the possible flights they could get on
  capacity_individual_flight[] <- capacity_per_flight
  n_inf_specific_flightAB[] <- rmhyper(n_inf_flightAB, capacity_individual_flight)

  # Stochastic model updates for outputs relevant to surveillance (number infected on flights etc.)
  update(n_inf_flightOut) <- n_inf_flight
  update(n_inf_flightABOut) <- n_inf_flightAB
  update(n_inf_specific_flightABOut[]) <- n_inf_specific_flightAB[i]

  # Initial values for outputs relevant to surveillance (number infected on flights etc.)
  initial(n_inf_flightOut) <- 0
  initial(n_inf_flightABOut) <- 0
  initial(n_inf_specific_flightABOut[]) <- 0

  # Metagenomic sequencing module ----

  # Metagenomic and sequencing parameters

  # Average number of defecation events per person per flight
  shedding_freq <- user()

  # Average amount of material shed per event for our virus of interest (defecation)
  virus_shed <- user()

  # Average amount of other nucleic acid (i.e. not virus of interest) shed per event (defecation)
  non_virus_shed <- user()

  # Volume of wastewater into which individuals shed (to convert abundance to concentration)
  volume_wastewater <- user()

  # Bias term for the metagenomic model
  bias <- user()

  # Total amount of sequencing that is done
  m_tot <- user()

  # Calculating the number of shedding events from infected and uninfected individuals
  infected_indiv_shedding_events[] <- rpois(n_inf_specific_flightAB[i] * shedding_freq)
  uninfected_indiv_shedding_events[] <- rpois((capacity_individual_flight[i] - n_inf_specific_flightAB[i]) * shedding_freq)

  # Calculating amount of nucleic acid shed into wastewater on each flight (i.e. abundance)
  amount_virus_indiv_flight[] <- infected_indiv_shedding_events[i] * virus_shed
  amount_non_virus_indiv_flight[] <- (uninfected_indiv_shedding_events[i] + infected_indiv_shedding_events[i]) * non_virus_shed

  # Converting the abundance of nucleic acid on each flight into the concentration
  conc_virus_indiv_flight[]  <- amount_virus_indiv_flight[i] / volume_wastewater
  conc_non_virus_indiv_flight[]  <- amount_non_virus_indiv_flight[i] / volume_wastewater

  # Converting nucleic acid concentration into sequencing reads
  m_i[] <- m_tot * (conc_virus_indiv_flight[i] * bias)/((conc_virus_indiv_flight[i] * bias) + conc_non_virus_indiv_flight[i])

  # Note: # nbinom to be implemented in odin shortly (not yet there / need to ask Rich about parameterising)
  stoch_m_i_pois[] <- rpois(m_i[i])
  stoch_m_not_i[] <- m_tot - stoch_m_i_pois[i]

  # Stochastic model updates for outputs relevant to metagenomic sequencing
  update(amount_virus_indiv_flight_Out[]) <- amount_virus_indiv_flight[i]
  update(amount_non_virus_indiv_flight_Out[]) <- amount_non_virus_indiv_flight[i]
  update(conc_virus_indiv_flight_Out[]) <- conc_virus_indiv_flight[i]
  update(conc_non_virus_indiv_flight_Out[]) <- conc_non_virus_indiv_flight[i]
  update(m_i_Out[]) <- m_i[i]
  update(stoch_m_i_pois_Out[]) <- stoch_m_i_pois[i] # replace with negbinom when implemented
  update(stoch_m_not_i_Out[]) <- stoch_m_not_i[i] # replace with negbinom when implemented

  # Initial values
  initial(amount_virus_indiv_flight_Out[]) <- 0
  initial(amount_non_virus_indiv_flight_Out[]) <- 0
  initial(conc_virus_indiv_flight_Out[]) <- 0
  initial(conc_non_virus_indiv_flight_Out[]) <- 0
  initial(m_i_Out[]) <- 0
  initial(stoch_m_i_pois_Out[]) <- 0
  initial(stoch_m_not_i_Out[]) <- 0

  # Miscellaneous model requirements ----

  # Definition of the time-step and output as "time"
  dt <- user(0.05)
  initial(time) <- 0
  update(time) <- (step + 1) * dt

  # Specifying the dimensions of the different vectors in the model (required for C compilation)
  dim(capacity_individual_flight) <- num_flightsAB
  dim(n_inf_specific_flightAB) <- num_flightsAB
  dim(n_inf_specific_flightABOut) <- num_flightsAB
  dim(infected_indiv_shedding_events) <- num_flightsAB
  dim(uninfected_indiv_shedding_events) <- num_flightsAB
  dim(amount_virus_indiv_flight) <- num_flightsAB
  dim(amount_non_virus_indiv_flight) <- num_flightsAB
  dim(conc_virus_indiv_flight) <- num_flightsAB
  dim(conc_non_virus_indiv_flight) <- num_flightsAB
  dim(m_i) <- num_flightsAB
  dim(stoch_m_i_pois) <- num_flightsAB
  dim(stoch_m_not_i) <- num_flightsAB
  dim(amount_virus_indiv_flight_Out) <- num_flightsAB
  dim(amount_non_virus_indiv_flight_Out) <- num_flightsAB
  dim(conc_virus_indiv_flight_Out) <- num_flightsAB
  dim(conc_non_virus_indiv_flight_Out) <- num_flightsAB
  dim(m_i_Out) <- num_flightsAB
  dim(stoch_m_i_pois_Out) <- num_flightsAB
  dim(stoch_m_not_i_Out) <- num_flightsAB

})

cut_limit <- function(x, type) {
  if(type == "lower") return(as.numeric(gsub(",.*","", gsub("\\(|\\[|\\)|\\]","", x))))
  if(type == "upper") return (as.numeric(gsub(".*,","", gsub("\\(|\\[|\\)|\\]","", x))))
  else return(warning('type must be either "lower" or "upper"'))
}

midpoints <- function(x, dp = 2){
  lower <- cut_limit(x, type = "lower")
  upper <- cut_limit(x, type = "uppper")
  return(round(lower + (upper - lower) / 2, dp))
}
