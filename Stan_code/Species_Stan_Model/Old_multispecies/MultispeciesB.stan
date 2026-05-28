// 2-state HMM for seed production - multispecies model
// Now with stand-level random effects for BOTH states

data {
  int<lower=1> N;                          // total observations
  int<lower=1> F;                          // number of stands
  int<lower=1> S;                          // number of species
  
  array[N] int<lower=0> y;                 // seed counts
  array[N] int<lower=1, upper=S> sp;       // species index per observation
  
  array[F] int<lower=1> start_idxs;        // start index per stand
  array[F] int<upper=N> end_idxs;          // end index per stand
  
  vector<lower=0>[N] area;                 // trap area per observation
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  simplex[2] rho;                          // initial state probabilities
  
  vector<lower=0, upper=1>[S] theta1;      // prob stay in low state per species
  vector<lower=0, upper=1>[S] theta2;      // prob stay in high state per species
  
  array[S] ordered[2] log_means;           // [s][1]=low mean, [s][2]=high mean
  
  // Stand-level random effects for BOTH states
  vector[F] stand_effect_raw_low;          // NEW: for low state
  vector[F] stand_effect_raw_high;         // renamed from stand_effect_raw
  
  // Separate sigmas for each state
  real<lower=0> sigma_low;                 // NEW: SD for low state
  real<lower=0> sigma_high;                // renamed from sigma
  
  vector[S] log_phi1;
  vector[S] log_phi2;
}

transformed parameters {
  vector<lower=0>[S] phi1 = exp(log_phi1);
  vector<lower=0>[S] phi2 = exp(log_phi2);
  
  // Stand-level means for BOTH states
  vector[F] log_alpha_low;                 // NEW: low state with stand effect
  vector[F] log_alpha_high;                // renamed from log_alpha
  
  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    log_alpha_low[f]  = log_means[s][1] + stand_effect_raw_low[f] * sigma_low;
    log_alpha_high[f] = log_means[s][2] + stand_effect_raw_high[f] * sigma_high;
  }
  
  array[S] matrix[2, 2] Gamma;
  for (s in 1:S) {
    Gamma[s][1, 1] = theta1[s];
    Gamma[s][1, 2] = 1 - theta1[s];
    Gamma[s][2, 1] = 1 - theta2[s];
    Gamma[s][2, 2] = theta2[s];
  }
  
  // Log-likelihoods — now both states use stand effects
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    for (t in start_id:end_id) {
      int s = sp[t];
      // LOW state now uses log_alpha_low[f] instead of log_means[s][1]
      log_omega[1, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_low[f]  + log_area_ratio[t], phi1[s]);
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_high[f] + log_area_ratio[t], phi2[s]);
    }
  }
}

model {
  rho ~ dirichlet([8.0, 2.0]');

  theta1 ~ beta(5, 1);
  theta2 ~ beta(1, 4);

  // Species-specific priors (unchanged)
  log_means[1][1] ~ normal(0, 0.8); 
  log_means[1][2] ~ normal(5.2, 0.6); 

  log_means[2][1] ~ normal(2.3, 0.8);  
  log_means[2][2] ~ normal(5.8, 0.8);   

  log_means[3][1] ~ normal(0, log(30)/2.57);
  log_means[3][2] ~ normal(5.2, 0.5);   

  log_means[4][1] ~ normal(4, 0.5);   
  log_means[4][2] ~ normal(7, 0.5);   

  log_means[5][1] ~ normal(4, 0.5);   
  log_means[5][2] ~ normal(7, 0.5);  

  // Priors for BOTH stand effect SDs
  sigma_low  ~ normal(0, 0.5);
  sigma_high ~ normal(0, 0.8);
  
  // Non-centered parameterization for both
  stand_effect_raw_low  ~ normal(0, 1);
  stand_effect_raw_high ~ normal(0, 1);

  log_phi1 ~ normal(log(2.0), 0.5);
  log_phi2 ~ normal(log(6.5), 0.5);

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    int s        = sp[start_id];
    target += hmm_marginal(log_omega[, start_id:end_id], Gamma[s], rho);
  }
}

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;
  
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    int s_stand  = sp[start_id];
    
    state[start_id:end_id] = hmm_latent_rng(
      log_omega[, start_id:end_id], Gamma[s_stand], rho
    );
    
    for (t in start_id:end_id) {
      int s = sp[t];
      if (state[t] == 1) {
        y_rep[t] = neg_binomial_2_log_rng(log_alpha_low[f]  + log_area_ratio[t], phi1[s]);
      } else {
        y_rep[t] = neg_binomial_2_log_rng(log_alpha_high[f] + log_area_ratio[t], phi2[s]);
      }
    }
  }
}



