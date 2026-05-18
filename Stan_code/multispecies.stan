

data {
  int<lower=1> N;
  int<lower=1> F;
  
  int<lower=1> N_species;       // number of species
  int<lower=1> N_stands;

  array[N] int<lower=0> y;
  array[N] int<lower=1, upper=N_species> sp;

  array[F] int<lower=1> start_idxs;
  array[F] int<upper=N> end_idxs;
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
  simplex[2] rho;

  // Transitions
  real grand_logit_theta1;
  real grand_logit_theta2;
  
  vector[N_species] alpha_theta1_species_tilde;   
  vector[N_species] alpha_theta2_species_tilde;  
  real<lower=0> sigma_theta1_species;
  real<lower=0> sigma_theta2_species;

  vector[N_stands] alpha_theta1_stand_tilde;
  vector[N_stands] alpha_theta2_stand_tilde;
  real<lower=0> sigma_theta1_stand;
  real<lower=0> sigma_theta2_stand;

  // Emission means; 
  real mu_log_low;
  real mu_log_delta;
  

  vector [N_species] log_alpha_low_species;
  vector [N_stands] log_alpha_low_stand_tilde;
  real<lower=0> sigma_low_species;
  real<lower=0> sigma_low_stand;

  vector[N_species] log_delta_species_tilde;
  vector[N_stands] log_delta_stand_tilde;
  real<lower=0> sigma_delta_species;
  real<lower=0> sigma_delta_stand;

  // Dispersion per state
  real<lower=0> phi_low;
  real<lower=0> phi_high;
}

transformed parameters {
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;
  
  vector[N_species] alpha_theta1_species = sigma_theta1_species * alpha_theta1_species_tilde;
  vector[N_stands] alpha_theta1_stand = sigma_theta1_stand * alpha_theta1_stand_tilde;
  
  vector[N_species] alpha_theta2_species = sigma_theta2_species * alpha_theta2_species_tilde;
  vector[N_stands] alpha_theta2_stand = sigma_theta2_stand * alpha_theta2_stand_tilde;
  
  vector[N_stands] log_alpha_low_stand = sigma_low_stand * log_alpha_low_stand_tilde;
  
  vector[N_species] log_delta_species = sigma_delta_species * log_delta_species_tilde;
  vector[N_stands] log_delta_stand = sigma_delta_stand * log_delta_stand_tilde;
  
  for (f in 1:F) {
    
    int s = sp[start_idxs[f]];
    int st = stand_id[f];

    log_alpha_low[f] = mu_log_low + log_alpha_low_species[s] + log_alpha_low_stand[st];
    log_alpha_high[f] = log_alpha_low[f] + log1p_exp(mu_log_delta + log_delta_species[s] + log_delta_stand[st]);

    theta1[f] = inv_logit(grand_logit_theta1 + alpha_theta1_species[s] + alpha_theta1_stand[st]);
    theta2[f] = inv_logit(grand_logit_theta2 + alpha_theta2_species[s] + alpha_theta2_stand[st]);

    Gamma[f][1, 1] = theta1[f];
    Gamma[f][1, 2] = 1 - theta1[f];
    Gamma[f][2, 1] = 1 - theta2[f];
    Gamma[f][2, 2] = theta2[f];
  }

  matrix[2, N] log_omega;
  for (f in 1:F) {
    for (t in start_idxs[f]:end_idxs[f]) {
      log_omega[1, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_low[f]  + log_area_ratio[t], phi_low);
      log_omega[2, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_high[f] + log_area_ratio[t], phi_high);
    }
  }
}

model {
  rho ~ dirichlet(rep_vector(2.0, 2));

  // Transitions
  grand_logit_theta1 ~ normal(1, 0.7);//mean of 73% chance of staying in low state
  grand_logit_theta2 ~ normal(0, 0.7); 

  alpha_theta1_species_tilde ~ normal(0, 1);  
  alpha_theta1_stand_tilde ~ normal(0, 1);   
  sigma_theta1_species ~ normal(0, 0.7);
  sigma_theta1_stand ~ normal(0, 0.7);
  
  alpha_theta2_species_tilde ~ normal(0, 1);  
  alpha_theta2_stand_tilde ~ normal(0, 1);  
  sigma_theta2_species ~ normal(0, 0.7);
  sigma_theta2_stand ~ normal(0, 0.7);

  mu_log_low               ~ normal(2, 2.0);
  log_alpha_low_species           ~ normal(0, sigma_low_species);
  log_alpha_low_stand_tilde           ~ normal(0, 1);
  sigma_low_species           ~ normal(0, 3);
  sigma_low_stand           ~ normal(0, 1);

  mu_log_delta              ~ normal(1.5, 1.5);
  log_delta_species_tilde           ~ normal(0, 1);
  log_delta_stand_tilde           ~ normal(0, 1);
  sigma_delta_species           ~ normal (0, 3);
  sigma_delta_stand           ~ normal (0,1);

  // Dispersion
  phi_low ~ normal(0, 1);
  phi_high ~ normal(0, 1);

  // Likelihood
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]], Gamma[f], rho);
}

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];

    state[start_id:end_id] = hmm_latent_rng(
        log_omega[, start_id:end_id], Gamma[f], rho);

    for (t in start_id:end_id) {
      if (state[t] == 1)
        y_rep[t] = neg_binomial_2_log_rng(
            log_alpha_low[f]  + log_area_ratio[t], phi_low);
      else
        y_rep[t] = neg_binomial_2_log_rng(
            log_alpha_high[f] + log_area_ratio[t], phi_high);
    }
  }
}

