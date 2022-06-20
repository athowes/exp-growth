// 4pl.stan

data {
  int T;                                     // Total number of cycles
  int<lower=0, upper=1> flag_run_estimation; // Should the likelihood be used?
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
    lambda[t] = n_inf +  ((n_0 - n_inf) / pow(1 + (t / m), beta));
  }
}

model {
  n_0 ~ normal(1, 0.1);
  n_inf ~ normal(10, 0.1);
  m ~ normal(10, 0.1);
  beta ~ normal(5, 0.1);
}

generated quantities {
  vector[T] n_sim;
  for(t in 1:T) {
    n_sim[t] = poisson_rng(lambda[t]);
  }
}
