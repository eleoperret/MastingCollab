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
  // Precomputed once — was inside transformed parameters before
  // which meant it was recomputed at every leapfrog step
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
  
  // Emission means now per species
  // In single-species: real log_lambda; real<lower=log_lambda> log_mu
  // Now: ordered[2] — enforces log_means[s][1] < log_means[s][2]
  // This replaces the lower= constraint which only works for scalars
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
  // Now:               array[S] matrix[2,2] Gamma
  array[S] matrix[2, 2] Gamma;
  for (s in 1:S) {
    Gamma[s][1, 1] = theta1[s];
    Gamma[s][1, 2] = 1 - theta1[s];
    Gamma[s][2, 1] = 1 - theta2[s];
    Gamma[s][2, 2] = theta2[s];
  }
  
  // Emission log-likelihoods — now uses species index s
  // In single-species: log_lambda, phi1 (scalars)
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

//model {
  // Initial state
  rho ~ dirichlet(rep_vector(1.0, 2));
  
  // Transitions — now vectorized over species
  // In single-species: loop over F stands
  // Now: single vectorized statement over S species
  theta1 ~ beta(3, 1);
  theta2 ~ beta(2, 2);
  
  // Emission priors — now per species
  // In single-species: log_lambda ~ normal(0, log(30)/2.57)
  //                    log_mu     ~ normal(log(200), 0.5)
  // Now: species-specific informed priors from histogram inspection
  // Adjust these based on which species you include
  for (s in 1:S) {
    log_means[s][1] ~ normal(2.0, 1.5);   // weakly informative low  state
    log_means[s][2] ~ normal(5.0, 1.5);   // weakly informative high state
  }
  
  // Stand random effects — unchanged
  sigma            ~ normal(0, 0.5);
  stand_effect_raw ~ normal(0, 1);
  
  // Dispersion — now on log scale per species
  // In single-species: phi1 ~ normal(6.5, 3)
  // Now: log scale prior, exp(log(6.5)) = 6.5 same centre but more stable
  log_phi1 ~ normal(log(6.5), 0.5);
  log_phi2 ~ normal(log(6.5), 0.5);
  
  // Likelihood — now uses Gamma[s] not Gamma[f]
  // In single-species: hmm_marginal(..., Gamma[f], rho)
  // Now:               hmm_marginal(..., Gamma[s], rho)
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    int s        = sp[start_id];           // species of this stand
    target += hmm_marginal(log_omega[, start_id:end_id], Gamma[s], rho);
  }
}


model {
  rho ~ dirichlet(rep_vector(1.0, 2));

  theta1 ~ beta(3, 1);
  theta2 ~ beta(2, 2);

  // Species-specific priors derived from quantile analysis
  // Logic: low state centre from q10-q25 range on log scale
  //        high state centre from q75-q90 range on log scale
  // SD = 0.8 for species with clear separation (ABAM, CANO, THPL)
  // SD = 1.0 for species with unclear separation (PSME, TSHE)

  // ABAM (s=1)
  // Low  state: q25=2   → log(5)  ≈ 1.6
  // High state: q75=33  → log(60) ≈ 4.1
  log_means[1][1] ~ normal(1.6, 0.8);
  log_means[1][2] ~ normal(4.1, 0.8);

  // CANO (s=2)
  // Low  state: q25=3    → log(10)  ≈ 2.3
  // High state: q90=224  → log(300) ≈ 5.7
  log_means[2][1] ~ normal(2.3, 0.8);
  log_means[2][2] ~ normal(5.7, 0.8);

  // PSME (s=3) — unclear separation, wider priors
  // Low  state: q25=6   → log(10)  ≈ 2.3
  // High state: q90=123 → log(150) ≈ 5.0
  log_means[3][1] ~ normal(2.3, 1.0);
  log_means[3][2] ~ normal(5.0, 1.0);

  // THPL (s=4)
  // Low  state: q25=10  → log(20)  ≈ 3.0
  // High state: q90=707 → log(800) ≈ 6.7
  log_means[4][1] ~ normal(3.0, 0.8);
  log_means[4][2] ~ normal(6.7, 0.8);

  // TSHE (s=5) — both states high, wider priors
  // Low  state: q10=21  → log(50)   ≈ 3.9
  // High state: q90=1668→ log(1000) ≈ 6.9
  log_means[5][1] ~ normal(3.9, 1.0);
  log_means[5][2] ~ normal(6.9, 1.0);

  sigma            ~ normal(0, 0.5);
  stand_effect_raw ~ normal(0, 1);

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