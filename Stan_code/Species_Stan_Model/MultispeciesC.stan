// 2-state HMM for seed production - multispecies model
// Expands single-species PSME model to multiple species
// Key changes from single-species:
  //   - log_lambda, log_mu, phi1, phi2 now vectors of length S
//   - theta1, theta2 now per species not per stand
//   - ordered[2] replaces lower=log_lambda constraint
//   - log_area_ratio precomputed in transformed data
//   - sp[] index added to data to identify species per observation

data {
  int<lower=1> N; // total observations
  int<lower=1> F; // number of stands x species combination
  int<lower=1> S; // number of species
  
  array[N] int<lower=0> y; // seed counts
  array[N] int<lower=1, upper=S> sp;// species index per observation
  
  array[F] int<lower=1> start_idxs; // start index per species x tand
  array[F] int<upper=N> end_idxs;   // end index per species x stand
  
  vector<lower=0>[N] area;// trap area per observation
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  simplex[2] rho;// initial state probabilities
  
  // Transitions now per species not per stand
  vector<lower=0, upper=1>[S] theta1; // prob stay in low  state per species
  vector<lower=0, upper=1>[S] theta2;   // prob stay in high state per species
  
  // Emission means now per species
  array[S] ordered[2] log_means;
  
  // Stand-level random effects — unchanged from single-species
  vector[F] stand_effect_low_raw;
  vector[F] stand_effect_high_raw;

  real<lower=0> sigma_low;
  real<lower=0> sigma_high;
  
  
  // Dispersion now per species,
  vector[S] log_phi1;
  vector[S] log_phi2;
}

transformed parameters {
  vector<lower=0>[S] phi1 = exp(log_phi1);
  vector<lower=0>[S] phi2 = exp(log_phi2);
  
  // Stand-level partial pooling for both state
  vector[F] log_alpha_high;
  vector[F] log_alpha_low;

  for (f in 1:F) {
    int s = sp[start_idxs[f]];
  
    log_alpha_low[f]  = log_means[s][1] + stand_effect_low_raw[f]  * sigma_low;
    log_alpha_high[f] = log_means[s][2] + stand_effect_high_raw[f] * sigma_high;
  }
  
  // Transition matrices per species not per stand
  array[S] matrix[2, 2] Gamma;
  for (s in 1:S) {
    Gamma[s][1, 1] = theta1[s];
    Gamma[s][1, 2] = 1 - theta1[s];
    Gamma[s][2, 1] = 1 - theta2[s];
    Gamma[s][2, 2] = theta2[s];
  }
  
  // Emission log-likelihoods
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    for (t in start_id:end_id) {
      int s = sp[t];    
      log_omega[1, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_low[f]  + log_area_ratio[t], phi1[s]);
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_high[f]    + log_area_ratio[t], phi2[s]);
    }
  }
}


model {
  rho ~ dirichlet(rep_vector(8.0, 2));

  theta1 ~ beta(4, 1);
  theta2 ~ beta(1, 4);

  // ABAM (s=1)
  log_means[1][1] ~ normal(1.6, 0.8);
  log_means[1][2] ~ normal(4.1, 0.8);

  // CANO (s=2)
  log_means[2][1] ~ normal(2.3, 0.8);
  log_means[2][2] ~ normal(5.7, 0.8);

  // PSME (s=3) 
  log_means[3][1] ~ normal(2.3, 1.0);
  log_means[3][2] ~ normal(5.0, 1.0);

  // THPL (s=4)
  log_means[4][1] ~ normal(3.0, 0.8);
  log_means[4][2] ~ normal(6.7, 0.8);

  // TSHE (s=5) 
  log_means[5][1] ~ normal(3.9, 1.0);
  log_means[5][2] ~ normal(6.9, 1.0);

  stand_effect_low_raw  ~ normal(0, 1);
  stand_effect_high_raw ~ normal(0, 1);

  sigma_low  ~ normal(0, 0.5);
  sigma_high ~ normal(0, 0.5);

  log_phi1 ~ normal(log(6.5), 0.5);
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
        y_rep[t] = neg_binomial_2_log_rng(log_alpha_low[f] + log_area_ratio[t], phi1[s]);
      } else {
        y_rep[t] = neg_binomial_2_log_rng(log_alpha_high[f]    + log_area_ratio[t], phi2[s]);
      }
    }
  }
}
