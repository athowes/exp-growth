// 4pl.stan

data {
  int T;                                     // Total number of cycles
  int<lower=0, upper=1> flag_run_estimation; // Should the likelihood be used?
  int n[T];                              // Observed counts
}

parameters {
  real<lower=0> n_0;   // Asymptote as t -> 0
  real<lower=0> n_inf; // Asymptote as t -> inf
  real<lower=0> m;     // Midpoint
  real<lower=0> beta;  // Slope
}

transformed parameters {
  vector[T] lambda;    // Intensities

  // 4PL logistic
  for(t in 1:T) {
    lambda[t] = n_inf +  ((n_0 - n_inf) / (1 + pow(t / m, beta)));
  }
}

model {
  // Prior
  n_0 ~ normal(1, 0.5);    // Give a plausible range for the initial amount
  n_inf ~ normal(10, 0.5); // Give a plausible range for the final amount
  m ~ normal(10, 0.5);
  beta ~ normal(5, 0.5);

  // Likelihood (evaluated if flag_run_estimation is TRUE)
  if(flag_run_estimation == 1){
    n ~ poisson(lambda);
  }
}

generated quantities {
  vector[T] n_sim;
  for(t in 1:T) {
    n_sim[t] = poisson_rng(lambda[t]);
  }
}
