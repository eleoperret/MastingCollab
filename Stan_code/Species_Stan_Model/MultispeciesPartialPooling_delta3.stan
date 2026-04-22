//HMM 2 state
//partial pooling on the distribution for the high and low state for species and stand
//NC high stands + theta 1 and thetha 2 stands

data {
  int<lower=1> N;
  int<lower=1> F;       // number of species x stand series
  int<lower=1> S;       // number of species
  int<lower=1> N_stands; // NEW : number of unique stands

  array[N] int<lower=0> y;
  array[N] int<lower=1, upper=S> sp;

  array[F] int<lower=1> start_idxs;
  array[F] int<upper=N> end_idxs;
  array[F] int<lower=1, upper=N_stands> stand_id; // NEW: stand index per series

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

  // Transitions (partially pooled by species and stand)
  real grand_logit_theta1; //grand mean for theta
  real grand_logit_theta2;
  //Non-centered now
  vector[S] alpha_theta1_species_nc;
  vector[S] alpha_theta2_species_nc;
  real<lower=0> sigma_theta1_species;
  real<lower=0> sigma_theta2_species;

  vector[N_stands] alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand_nc;
  real<lower=0> sigma_theta1_stand;
  real<lower=0> sigma_theta2_stand;

  // Grand means (intercepts) — one per state, replaces ordered log_means
  real grand_mean_low;
  real log_delta_high_grand_mean;

  // Species random effects 
  vector[S] alpha_low_species;
  vector[S] log_delta_high_species;
  real<lower=0> sigma_low_species;
  real<lower=0> sigma_log_delta_high_species;

  // Stand random effects 
  vector[N_stands] alpha_low_stand;
  vector[N_stands] log_delta_high_stand_nc;
  real<lower=0> sigma_low_stand;
  real<lower=0> sigma_log_delta_high_stand;

  // Dispersion (not partially pooled but per species)
  real<lower=0> phi1;
  real<lower=0> phi2;
}

transformed parameters {
  
  vector[N_stands] log_delta_high_stand = log_delta_high_stand_nc * sigma_log_delta_high_stand;
  
  // Random effects: grand_mean + species effect + stand effect
  vector[F] log_alpha_low;
  vector[F] log_delta_high;
  vector[F] log_alpha_high;
  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    int st = stand_id[f];
    log_alpha_low[f]  = grand_mean_low  + alpha_low_species[s]  + alpha_low_stand[st];
    log_delta_high[f] = log_delta_high_grand_mean + log_delta_high_species[s] 
    + log_delta_high_stand[st];
    log_alpha_high[f] = log_sum_exp(log_alpha_low[f], log_delta_high[f]);
  }

   //Non-centré Theta
  vector[S] alpha_theta1_species = sigma_theta1_species * alpha_theta1_species_nc;
  vector[S] alpha_theta2_species = sigma_theta2_species * alpha_theta2_species_nc;
  
  vector[N_stands] alpha_theta1_stand = sigma_theta1_stand * alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand = sigma_theta2_stand * alpha_theta2_stand_nc;
  
  
  // Random effects for transitions
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

  // Emission log-likelihoods 
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    for (t in start_id:end_id) {
      int s= sp[t];
      log_omega[1, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_low[f]  + log_area_ratio[t], phi1);
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_high[f] + log_area_ratio[t], phi2);
    }
  }
}

model {
  rho ~ dirichlet(rep_vector(8.0, 2));

  // Grand means
  grand_logit_theta1 ~ normal(1.4, 1); //changed those priors here too
  grand_logit_theta2 ~ normal(-1.4, 1);

  alpha_theta1_species_nc ~ normal(0,1);
  alpha_theta2_species_nc ~ normal(0,1);
  sigma_theta1_species ~ normal(0, 0.5);
  sigma_theta2_species ~ normal(0, 0.5);

  alpha_theta1_stand_nc ~ normal(0, 1);
  alpha_theta2_stand_nc ~ normal(0, 1);
  sigma_theta1_stand ~ normal(0, 0.5);
  sigma_theta2_stand ~ normal(0, 0.5);

  // Grand means — set priors to roughly the centre of your old per-species priors
  grand_mean_low  ~ normal(3, 2);//changed that

  // Species random effects
  alpha_low_species  ~ normal(0,sigma_low_species);
  sigma_low_species  ~ normal(0, 2);

  // Stand random effects
  alpha_low_stand  ~ normal(0,  sigma_low_stand);
  sigma_low_stand  ~ normal(0, 0.5);
  
  // 
  log_delta_high_grand_mean ~ normal(2.5, 1.0); //normal(3, 1);
  sigma_log_delta_high_species ~ normal(0, 1);
  log_delta_high_species ~ normal(0, sigma_log_delta_high_species);
  sigma_log_delta_high_stand ~ normal(0, 1);
  log_delta_high_stand_nc ~ normal(0, 1);

  // Dispersion
  phi1    ~ gamma(4.0, 0.6);
  phi2   ~ gamma(4.0, 0.6);

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    target += hmm_marginal(log_omega[, start_id:end_id], Gamma[f], rho);
  }
}

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;
  vector[N_stands]               stand_sync;

  // --- state sequence + y_rep ---
  for (f in 1:F) {
    int f_start = start_idxs[f];  // renamed
    int f_end   = end_idxs[f];    // renamed

    state[f_start:f_end] = hmm_latent_rng(
      log_omega[, f_start:f_end], Gamma[f], rho
    );

    for (t in f_start:f_end) {
      if (state[t] == 1)
        y_rep[t] = neg_binomial_2_log_rng(log_alpha_low[f]  + log_area_ratio[t], phi1);
      else
        y_rep[t] = neg_binomial_2_log_rng(log_alpha_high[f] + log_area_ratio[t], phi2);
    }
  }

  // --- synchrony ---
  {
    array[N_stands] real sum_prop;
    array[N_stands] int  stand_n_years;

    for (st in 1:N_stands) {
      sum_prop[st]      = 0.0;
      stand_n_years[st] = 0;
    }

    for (st in 1:N_stands) {
      int f_stand = 0;
      for (f in 1:F)
        if (stand_id[f] == st)
          f_stand = f;

      stand_n_years[st] = end_idxs[f_stand] - start_idxs[f_stand] + 1;

      for (yr in 1:stand_n_years[st]) {
        int count_low  = 0;
        int count_high = 0;
        for (f in 1:F) {
          if (stand_id[f] == st) {
            int t = start_idxs[f] + yr - 1;
            if (state[t] == 1) count_low  += 1;
            else                count_high += 1;
          }
        }
        sum_prop[st] += max(count_low, count_high) * 1.0 / (count_low + count_high);
      }
    }

    for (st in 1:N_stands)
      stand_sync[st] = sum_prop[st] / stand_n_years[st];
  }
}

