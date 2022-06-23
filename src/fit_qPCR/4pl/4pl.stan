// 4pl.stan

functions {
  real four_parameter_logistic(real t, real n_0, real n_inf, real m, real b) {
    return(n_inf +  ((n_0 - n_inf) / (1 + pow(t / m, b))));
  }
}

data {
  int T;                                     // Total number of cycles
  int<lower=0, upper=1> flag_run_estimation; // Should the likelihood be used?
  real f[T];                                  // Observed fluorescence data
  // real<lower=0> threshold;                    // Threshold value
}

parameters {
  // 4PL
  // real<lower=0> n_0;    // Asymptote as t -> 0 (assume this is zero)
  real<lower=0> n_inf;     // Asymptote as t -> inf
  real<lower=0> m;         // Midpoint
  real<lower=0> b;         // Slope

  // Link to fluoresence
  // real<lower=0> beta_0; // Background fluoresence (assume this is removed)
  real<lower=0> beta;      // Slope (assume counts proportional to fluoresence)
  real<lower=0> sigma_f;   // Standard deviation of fluoresence
}

transformed parameters {
  vector[T] n; // Counts

  // 4PL logistic
  for(t in 1:T) {
    n[t] = four_parameter_logistic(t, 0, n_inf, m, b);
  }
}

model {
  // Prior
  // n_0 ~ normal(1, 0.5); // A plausible range for the initial amount
  n_inf ~ normal(8, 0.5);  // A plausible range for the final amount
  m ~ normal(20, 5);       // The midpoint at around 20 with high uncertainty
  b ~ normal(10, 1);       // And the slope at around 10

  // beta_0 ~ normal(0, 0.5);
  beta ~ normal(1, 0.1);     // Assume to be a 1-to-1 relationship approximately
  sigma_f ~ normal(0, 0.1); // Some amount of signal noise in fluoresence measurements

  // Likelihood (evaluated if flag_run_estimation is TRUE)
  if(flag_run_estimation == 1){
    f ~ normal(beta * n, sigma_f);
  }
}

generated quantities {
  real f_sim[T];
  f_sim = normal_rng(beta * n, sigma_f);
}
