// 4pl-processing.stan

functions {
  real four_parameter_logistic(real t, real n_0, real n_inf, real m, real b) {
    return(n_inf +  ((n_0 - n_inf) / (1 + pow(t / m, b))));
  }

  real xbinomial_logit_lpdf(real y, real m, real eta) {
    real dens = lchoose(m, y) + y * log(inv_logit(eta)) + (m - y) * log(1 - inv_logit(eta));
    real jacobian = log(inv_logit(eta)) + log(inv_logit(-eta)); // Alternative: eta - 2 * log(1 + exp(eta))
    return dens + jacobian;
  }
}

data {
  int T;                                     // Total number of cycles
  int<lower=0, upper=1> flag_run_estimation; // Should the likelihood be used?
  real f[T];                                 // Observed fluorescence data
}

parameters {
  // Processing
  real<lower=0> N_0;     // Pre-processing initial copies
  real r_logit;          // Proportion DNA successfully released during processing

  // 4PL
  real<lower=0> n_0;     // Asymptote as t -> 0
  real<lower=0> n_inf;   // Asymptote as t -> inf
  real<lower=0> m;       // Midpoint
  real<lower=0> b;       // Slope

  // Link to fluoresence
  real<lower=0> beta_0;  // Intercept (assume positive intercept)
  real<lower=0> beta;    // Slope (assume count fluoresence relationship is positive)
  real<lower=0> sigma_f; // Standard deviation
}

transformed parameters {
  vector[T] n; // Counts
  for(t in 1:T) {
    n[t] = four_parameter_logistic(t, n_0, n_inf, m, b);
  }
}

model {
  // Prior

  // Processing
  N_0 ~ normal(1, 0.1);
  r_logit ~ normal(-3, 0.1);
  n_0 ~ xbinomial_logit(N_0, r_logit);

  // 4PL logistic
  n_inf ~ normal(10, 0.5); // Give a plausible range for the final amount
  m ~ normal(10, 0.5);
  b ~ normal(5, 0.5);

  // Link to fluoresence
  beta_0 ~ normal(0, 0.5);
  beta ~ normal(1, 0.1);
  sigma_f ~ normal(0, 0.25);

  // Likelihood (evaluated if flag_run_estimation is TRUE)

  if(flag_run_estimation == 1){
    f ~ normal(beta_0 + beta * n, sigma_f);
  }
}

generated quantities {
  real<lower=0, upper=1> r = inv_logit(r_logit);

  real f_sim[T];
  f_sim = normal_rng(beta_0 + beta * n, sigma_f);
}
