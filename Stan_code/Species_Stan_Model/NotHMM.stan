data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> N_stands;

  array[N] int<lower=0> y;
  array[N] int<lower=1, upper=S> sp;
  array[N] int<lower=1, upper=N_stands> stand_id;

  vector<lower=0>[N] area;
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  // Global intercept
  real grand_mean;

  // Species effects
  vector[S] alpha_species;
  real<lower=0> sigma_species;

  // Stand effects
  vector[N_stands] alpha_stand;
  real<lower=0> sigma_stand;

  // Dispersion
  real<lower=0> phi;
}

model {
  // Priors
  grand_mean ~ normal(2.6, 1.0);

  alpha_species ~ normal(0, sigma_species);
  sigma_species ~ normal(0, 0.5);

  alpha_stand ~ normal(0, sigma_stand);
  sigma_stand ~ normal(0, 0.5);

  phi ~ gamma(4.0, 0.6);

  // Likelihood
  for (t in 1:N) {
    real log_mu = grand_mean
                  + alpha_species[sp[t]]
                  + alpha_stand[stand_id[t]]
                  + log_area_ratio[t];

    y[t] ~ neg_binomial_2_log(log_mu, phi);
  }
}

generated quantities {
  array[N] int y_rep;

  for (t in 1:N) {
    real log_mu = grand_mean
                  + alpha_species[sp[t]]
                  + alpha_stand[stand_id[t]]
                  + log_area_ratio[t];

    y_rep[t] = neg_binomial_2_log_rng(log_mu, phi);
  }
}
