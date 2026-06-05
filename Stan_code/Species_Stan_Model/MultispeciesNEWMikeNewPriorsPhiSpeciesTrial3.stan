//HMM 2 state - Multi species
//Transition matrix partially pooled by species and stand (non-centered)
//Distribution: two NB for each state partially pooled by species and stand (non-centered)
//Overdispersion per state and per species
//Rho: per species (simplex array over S)
//Forced state 2 above state 1 with log1p_exp(delta) — aligned with single-species model

data {
  int<lower=1> N;
  int<lower=1> F;         // number of species x stand series
  int<lower=1> S;         // number of species
  int<lower=1> N_stands;  // number of unique stands

  array[N] int<lower=0> y;
  array[N] int<lower=1, upper=S> sp;

  array[F] int<lower=1>              start_idxs;
  array[F] int<upper=N>              end_idxs;
  array[F] int<lower=1, upper=N_stands> stand_id;

  vector<lower=0>[N] area;
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  // Initial state distribution 
  simplex[2] rho;

  // Transitions (pooled by species and stand) 
  real grand_logit_theta1;
  real grand_logit_theta2;

  vector[S]       alpha_theta1_species_nc;
  vector[S]       alpha_theta2_species_nc;
  real<lower=0>   sigma_theta1_species;
  real<lower=0>   sigma_theta2_species;

  vector[N_stands] alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand_nc;
  real<lower=0>    sigma_theta1_stand;
  real<lower=0>    sigma_theta2_stand;

  // Emission means (pooled by species and stand) 
  // Low state: grand mean + species effect + stand effect
  real             mu_log_low;
  vector[S]        alpha_low_species_nc;
  real<lower=0>    sigma_low_species;
  vector[N_stands] alpha_low_stand_nc;
  real<lower=0>    sigma_low_stand;

  // Delta (log-scale gap between states): grand mean + species  + stand 
  //real             mu_log_delta;
  //vector[S]        log_delta_species_nc;
  //real<lower=0>    sigma_log_delta_species;
  //vector[N_stands] log_delta_stand_nc;
  //real<lower=0>    sigma_log_delta_stand;
  real <lower= 0> delta;

  // Dispersion: one per state
  vector [S] log_phi_species;
  vector [2] log_phi_state;
}

transformed parameters {
  // Non-centered species and stand effects 
  vector[S]        alpha_theta1_species = sigma_theta1_species * alpha_theta1_species_nc;
  vector[S]        alpha_theta2_species = sigma_theta2_species * alpha_theta2_species_nc;
  vector[N_stands] alpha_theta1_stand   = sigma_theta1_stand   * alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand   = sigma_theta2_stand   * alpha_theta2_stand_nc;

  vector[S]        alpha_low_species    = sigma_low_species     * alpha_low_species_nc;
  vector[N_stands] alpha_low_stand      = sigma_low_stand       * alpha_low_stand_nc;

  //vector[S]        log_delta_species    = sigma_log_delta_species * log_delta_species_nc;
  //vector[N_stands] log_delta_stand      = sigma_log_delta_stand   * log_delta_stand_nc;

  //Per-series emission means 
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;

  for (f in 1:F) {
    int s  = sp[start_idxs[f]]; 
    int st = stand_id[f];

    log_alpha_low[f]  = mu_log_low
                      + alpha_low_species[s]
                      + alpha_low_stand[st];

    log_alpha_high[f] = log_alpha_low[f]
                  + delta;
                  //+ log1p_exp(mu_log_delta)
                  //+ log1p_exp(log_delta_species[s])
                  //+ log1p_exp(log_delta_stand[st])
                  ;
  }

  // --- Per-series transition matrices ---
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;

  for (f in 1:F) {
    int s  = sp[start_idxs[f]];
    int st = stand_id[f];

    theta1[f] = inv_logit(grand_logit_theta1
                         + alpha_theta1_species[s]
                         + alpha_theta1_stand[st]);
    theta2[f] = inv_logit(grand_logit_theta2
                         + alpha_theta2_species[s]
                         + alpha_theta2_stand[st]);

    Gamma[f][1, 1] = theta1[f];
    Gamma[f][1, 2] = 1 - theta1[f];
    Gamma[f][2, 1] = 1 - theta2[f];
    Gamma[f][2, 2] = theta2[f];
  }

  // --- Emission log-likelihoods ---
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int s = sp[start_idxs[f]];  
    for (t in start_idxs[f]:end_idxs[f]) {
      log_omega[1, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1]));
      log_omega[2, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2]));
    }
  }
}

model {
  // --- Initial state ---
  rho ~ dirichlet(rep_vector(2.0, 2));

  // --- Transitions ---
  grand_logit_theta1 ~ normal(2,   0.5); //Changed from the previous version MultispeciesNEWMike.stan 
  grand_logit_theta2 ~ normal(1,   0.5); //Changed from the previous version MultispeciesNEWMike.stan 

  alpha_theta1_species_nc ~ normal(0, 1);
  alpha_theta2_species_nc ~ normal(0, 1);
  sigma_theta1_species    ~ normal(0, 1);//Changed from the previous version MultispeciesNEWMike.stan
  sigma_theta2_species    ~ normal(0, 1);//Changed from the previous version MultispeciesNEWMike.stan

  alpha_theta1_stand_nc ~ normal(0, 1);
  alpha_theta2_stand_nc ~ normal(0, 1);
  sigma_theta1_stand    ~ normal(0, 0.7);//Changed from the previous version MultispeciesNEWMike.stan
  sigma_theta2_stand    ~ normal(0, 0.7);//Changed from the previous version MultispeciesNEWMike.stan

  // Emission means
  mu_log_low           ~ normal(2.8, 1.0);
  alpha_low_species_nc  ~ normal(0, 1);
  sigma_low_species    ~ normal(0,2);//Changed from the previous version MultispeciesNEWMike.stan
  alpha_low_stand_nc   ~ normal(0, 1);
  sigma_low_stand      ~ normal(0.5, 1);//Changed from the previous version MultispeciesNEWMike.stan

  //mu_log_delta           ~ normal(1.5, 1.5);
  //log_delta_species_nc   ~ normal(0, 1);
  //sigma_log_delta_species ~ normal(0, 1);
  //log_delta_stand_nc     ~ normal(0, 1);
  //sigma_log_delta_stand  ~ normal(0, 1);
  delta ~ normal (1,0.3);
  
  // Dispersion
  log_phi_species ~ normal(log(4), 0.6); 
  log_phi_state ~ normal(0, 0.3); 

  // HMM marginal likelihood
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]],
        Gamma[f],
        rho);
}

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];

    state[start_id:end_id] = hmm_latent_rng(
        log_omega[, start_id:end_id],
        Gamma[f],
        rho);

    for (t in start_id:end_id) {
  int s = sp[t];  
  if (state[t] == 1)
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1]));
  else
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2]));
    }
  }
}


