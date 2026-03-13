// 2-state HMM for seed production - multispecies model
// Expands single-species PSME model to multiple species
// Key changes from single-species:
//   - log_lambda, log_mu, phi1, phi2 now vectors of length S
//   - theta1, theta2 now per species not per stand
//   - ordered[2] replaces lower=log_lambda constraint
//   - log_area_ratio precomputed in transformed data
//   - sp[] index added to data to identify species per observation

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
  // Precomputed once — was inside transformed parameters before which meant it was recomputed at every leapfrog step
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  simplex[2] rho;                          // initial state probabilities
  
  // Transitions now per species not per stand
  // In single-species: theta1[F], theta2[F]
  // Now:               theta1[S], theta2[S]
  vector<lower=0, upper=1>[S] theta1;      // prob stay in low  state per species
  vector<lower=0, upper=1>[S] theta2;      // prob stay in high state per species
  
  // In single-species: real log_lambda; real<lower=log_lambda> log_mu
  // Now: ordered[2] — enforces log_means[s][1] < log_means[s][2]
  // This replaces the lower= constraint which only works for scalars?
  array[S] ordered[2] log_means;          // [s][1]=low mean, [s][2]=high mean
  
  // Stand-level random effects — unchanged from single-species
  vector[F] stand_effect_raw;
  real<lower=0> sigma;
  
  // Dispersion now per species, on log scale for numerical stability
  // In single-species: real<lower=0> phi1; real<lower=0> phi2
  // Now: log scale vectors, back-transformed in transformed parameters
  vector[S] log_phi1;
  vector[S] log_phi2;
}

transformed parameters {
  // Back-transform dispersion — always positive by construction
  vector<lower=0>[S] phi1 = exp(log_phi1);
  vector<lower=0>[S] phi2 = exp(log_phi2);
  
  // Stand-level high state mean
  // In single-species: log_alpha[f] = log_mu + stand_effect_raw[f] * sigma
  // Now: uses log_means[s][2] where s is the species of stand f
  vector[F] log_alpha;
  for (f in 1:F) {
    int s = sp[start_idxs[f]];            // species of this stand
    log_alpha[f] = log_means[s][2] + stand_effect_raw[f] * sigma;
  }
  
  // Transition matrices now per species not per stand
  // In single-species: matrix[2,2] Gamma[F]
  // Now: array[S] matrix[2,2] Gamma
  array[S] matrix[2, 2] Gamma;
  for (s in 1:S) {
    Gamma[s][1, 1] = theta1[s];
    Gamma[s][1, 2] = 1 - theta1[s];
    Gamma[s][2, 1] = 1 - theta2[s];
    Gamma[s][2, 2] = theta2[s];
  }
  
  // Log-likelihoods — now uses species index s
  // In single-species: log_lambda, phi1 
  // Now:               log_means[s][1], phi1[s] (species-indexed)
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    for (t in start_id:end_id) {
      int s = sp[t];                       // species of this observation
      log_omega[1, t] = neg_binomial_2_log_lpmf(y[t] | log_means[s][1] + log_area_ratio[t], phi1[s]);
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha[f]    + log_area_ratio[t], phi2[s]);
    }
  }
}




model {
  rho ~ dirichlet([8.0, 2.0]');  // strongly favor starting in low state

  theta1 ~ beta(5, 1);   // mean = 0.83 — strong tendency to stay in LOW (non-masting is the norm)
  theta2 ~ beta(1, 4);   // mean = 0.20 — masting rarely persists year-to-year

  // Species-specific priors 

  // ABAM: q90=88 → log(88) ≈ 4.5
  log_means[1][1] ~ normal(0, 0.8); 
  log_means[1][2] ~ normal(5.2, 0.6); 

  // CANO: q90=224 → log(224) ≈ 5.4
  log_means[2][1] ~ normal(2.3, 0.8);  
  log_means[2][2] ~ normal(5.8, 0.8);   

  // PSME: q90=123 → log(123) ≈ 4.8
  log_means[3][1] ~ normal(0, log(30)/2.57);
  log_means[3][2] ~ normal(5.2, 0.5);   

  // THPL: q90=707 → log(707) ≈ 6.6
  log_means[4][1] ~ normal(4,0.5);   
  log_means[4][2] ~ normal(7,0.5);   

  // TSHE: q90=1668 → log(1668) ≈ 7.4
  log_means[5][1] ~ normal(4,0.5);   
  log_means[5][2] ~ normal(7,0.5) ;  

  sigma            ~ normal(0, 0.5);
  stand_effect_raw ~ normal(0, 1);

  log_phi1 ~ normal(log(2.0), 0.5);   // more overdispersed low state
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
    
    // In single-species: hmm_latent_rng(..., Gamma[f], rho)
    // Now:               hmm_latent_rng(..., Gamma[s_stand], rho)
    state[start_id:end_id] = hmm_latent_rng(
      log_omega[, start_id:end_id], Gamma[s_stand], rho
    );
    
    for (t in start_id:end_id) {
      int s = sp[t];
      if (state[t] == 1) {
        y_rep[t] = neg_binomial_2_log_rng(log_means[s][1] + log_area_ratio[t], phi1[s]);
      } else {
        y_rep[t] = neg_binomial_2_log_rng(log_alpha[f]    + log_area_ratio[t], phi2[s]);
      }
    }
  }
}
