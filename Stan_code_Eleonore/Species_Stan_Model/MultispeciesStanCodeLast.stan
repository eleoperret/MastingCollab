//HMM model for seed production at Mt. Rainier (multispecies)
//Eléonore Perret + Victor Van der Meersch + Elizabeth Wolkovich + Janneke Hille Ris Lambers

//August 2026 (Last data modified)

// We have a species x stand combination of seed production over multiple years. This seed production switches from state 1 (low or base seed production) to state 2 (high/higher seed production or masting). The seed production is following a negative binomial (as by its nature of seed count and highly overdispersed).

//Model : 
// Two latent state

//Data specification: As not all stands were surveyed starting at the same year, we used an ragged vector to establish the start and end observation for a specific stand. We also added for SUNR missing year (2023), missing observations for the 3 present species and finally we also pooled together the seed counts from all traps present at a site that we then divided by the total area of seed trap (see transformed data + negative binomial model)

// Seed production:
// Species and stand effect partial pooling on seed production mean
// Low state partial pooling centered and = global low mean + species effect + stand effect
// High state partial pooling non-centered for the species and stands and to enforce that state 1 is always smaller than state 2, we use the delta between the two states seed production (log1p_exp(delta). = low state mean  + delta positif.
//for non-centering, we used visual plots to diagnose funelling and non-centered accordingly

//Transition matrix:

//Rho: probability of starting in state 1 or 2.  

//The transition matrix is defined with theta 1 (probability of remaining in low seed production) and theta 2 (probability of remaining in high seed production)
// It is partially pooled by species ONLY (non-centered)

//Overdispersion per state and per species


//DATA

data {
  int<lower=1> N;                                   // total number of observations
  int<lower=1> F;                                   // number of species x stand series
  int<lower=1> S;                                   // number of species
  int<lower=1> N_stands;                            // number of unique stands

  array[N] int<lower=0> y;                          //Seed counts

  array[F] int<lower=1> start_idxs;
  array[F] int<upper=N> end_idxs;
  array[F] int<lower=1, upper=N_stands> stand_id;
  array[F] int<lower=1, upper=S> species_id;

  vector<lower=0>[N] area;
  
  array [N] int <lower= 0, upper = 1> obs_missing;
}


//TRANSFORMED DATA

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}


//PARAMETERS

parameters {
  // Initial state probability  
  simplex[2] rho; // rho = 1 , probability of stating in the low state and rho= 2 probability of starting in the high state

  // Transitions (pooled by species and non-centered) 
  real                   grand_logit_theta1;
  real                   grand_logit_theta2;

  vector[S]              alpha_theta1_species_nc;// species deviation from the global transition probability for theta 1 and non-centered 
  vector[S]              alpha_theta2_species_nc;// same but for theta 2
  real<lower=0>          sigma_theta1_species;// variation among species in theta 1
  real<lower=0>          sigma_theta2_species; // same for theta 2

  // Emission means (pooled by species and stand) 
  // Low state: grand mean + species effect + stand effect (centered)
  real                   mu_log_low;//global mean of seed production
  vector[S]              alpha_low_species;// species deviation from the low state mean
  real<lower=0>          sigma_low_species;//variation among species
  matrix [S, N_stands]   alpha_stand;//stand deviation from the low state mean
  real <lower=0>         sigma_stand;//variation among stands

  // Delta (log-scale gap between states,): grand mean + species effect + stand effect
  real                   mu_log_delta;//difference between low and high state
  vector[S]              log_delta_species_nc;//deviation per species non-centered
  real<lower=0>          sigma_log_delta_species;//variation amon species
  matrix [S, N_stands]   log_delta_stand_nc;//deviation per stands non-centered
  real <lower=0>         sigma_delta_stand;//variation among stands

  // Dispersion: (for the negative binomial), per state and species
  vector [S]             log_phi_species;//species specific overdispersion
  vector [2]             log_phi_state;//state specific overdispersion
}


//TRANSFORMED PARAMETERS
//what will be actually be used by my likelihood

transformed parameters {
  // Non-centered species effects on transition probabilities
  vector[S] alpha_theta1_species = sigma_theta1_species * alpha_theta1_species_nc;
  vector[S] alpha_theta2_species = sigma_theta2_species * alpha_theta2_species_nc;

  // Non-centered species and stands effects on seed production for the delta (~high seed production state)
  vector[S] log_delta_species    = sigma_log_delta_species * log_delta_species_nc;
  matrix[S, N_stands] log_delta_stand   = sigma_delta_stand * log_delta_stand_nc;
  
  //Per-series emission means (seed production)
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;

  for (f in 1:F) {
    int s  = species_id[f]; 
    int st = stand_id[f];

    log_alpha_low[f]  = mu_log_low
                      + alpha_low_species[s]
                      + alpha_stand[s,st];

    log_alpha_high[f] = log_alpha_low[f]
                  + log1p_exp(mu_log_delta)
                  + log1p_exp(log_delta_species[s])
                  + log1p_exp(log_delta_stand[s,st]);
  }

  //  Per-series transition matrices
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;

  for (f in 1:F) {
    int s  = species_id[f];

    theta1[f] = inv_logit(grand_logit_theta1
                         + alpha_theta1_species[s])
                         ;
    theta2[f] = inv_logit(grand_logit_theta2
                         + alpha_theta2_species[s])
                      ;

    Gamma[f][1, 1] = theta1[f];
    Gamma[f][1, 2] = 1 - theta1[f];
    Gamma[f][2, 1] = 1 - theta2[f];
    Gamma[f][2, 2] = theta2[f];
  }

  // Emission log-likelihoods
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int s = species_id[f];  
    for (t in start_idxs[f]:end_idxs[f]) {
      if (obs_missing[t]){ //Missing observations do not provide information 
        log_omega[1,t]= 0;
        log_omega [2,t] =0;
      } else {
        log_omega[1, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1]));
      log_omega[2, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2]));
      }
    }
  }
}


//MODEL

model {
  // Initial state 
  rho                       ~ dirichlet(rep_vector(2.0, 2));

  //Transitions 
  grand_logit_theta1        ~ normal(0.5, 1); 
  grand_logit_theta2        ~ normal(-1, 0.5); 

  alpha_theta1_species_nc   ~ normal(0, 1);
  alpha_theta2_species_nc   ~ normal(0, 1);
  sigma_theta1_species      ~ normal(0, 0.7);
  sigma_theta2_species      ~ normal(0, 0.7);

  // Emission means
  mu_log_low                ~ normal(2, 1.5);
  alpha_low_species         ~ normal(0, sigma_low_species);
  sigma_low_species         ~ normal(0.5, 1);
  for (s in 1:S) 
    alpha_stand [s]         ~ normal(0, sigma_stand);
  sigma_stand               ~ normal(1, 0.5);

  mu_log_delta              ~ normal(1.5, 1.5);
  log_delta_species_nc      ~ normal(0, 1);
  sigma_log_delta_species   ~ normal(1.5, 0.7);
  for (s in 1:S)
    log_delta_stand_nc [s]  ~ normal(0, 1);
  sigma_delta_stand         ~ normal(1, 0.5); 

  // Dispersion
  log_phi_species           ~ normal(log(4), 0.6); 
  log_phi_state             ~ normal(0, 0.3); 

  // HMM marginal likelihood
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]],
        Gamma[f],
        rho);
}


//GENERATED QUANTILES

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
  int s = species_id[f];  
  
  if (obs_missing[t]) {
    y_rep[t] = 0;
  } else if (state[t] == 1) {
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1]));
  }else{
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2]));
      }
    }
  }
}

