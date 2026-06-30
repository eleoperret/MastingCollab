// =============================================================================
// 2-state HMM for seed production — multispecies model
// =============================================================================
// CHANGES FROM v1 (in order of implementation):
//
// [1] WIDENED PRIORS
//     - sigma_low  : normal(0, 0.5) → normal(0, 1.0)
//     - sigma_high : normal(0, 0.8) → normal(0, 1.0)
//     - log_phi1   : normal(log(2.0), 0.5) → normal(log(2.0), 1.0)
//     - log_phi2   : normal(log(6.5), 0.5) → normal(log(6.5), 1.0)
//     - theta1     : beta(5, 1) → beta(3, 1)   [less strong persistence prior]
//     - theta2     : beta(1, 4) → beta(1, 3)   [less strong persistence prior]
//
// [2] PARTIALLY POOLED OVERDISPERSION across species
//     - log_phi1[s] ~ normal(mu_log_phi1, sigma_log_phi1)  [low state]
//     - log_phi2[s] ~ normal(mu_log_phi2, sigma_log_phi2)  [high state]
//     - Hyperpriors on mu and sigma for each state
//     - Species with little data borrow strength from others
//
// [3] SPECIES-SPECIFIC INITIAL STATE PROBABILITIES (rho)
//     - rho was a single shared simplex[2] → now array[S] simplex[2]
//     - Partial pooling via Dirichlet hyperprior:
//         rho[s] ~ dirichlet(rho_conc)
//         rho_conc ~ gamma(2, 1)   [hyperprior on concentration]
//     - Each species can have its own baseline masting frequency
//
// TODO (not yet implemented, lower priority):
// [4] Correlated stand effects across states (bivariate normal / Cholesky)
//     — would capture that a productive stand tends to be high in both states
// =============================================================================

data {
  int<lower=1> N;                          // total observations
  int<lower=1> F;                          // number of series (spp x stand)
  int<lower=1> S;                          // number of species

  array[N] int<lower=0> y;                 // seed counts
  array[N] int<lower=1, upper=S> sp;       // species index per observation

  array[F] int<lower=1> start_idxs;        // start index per series
  array[F] int<upper=N> end_idxs;          // end index per series

  vector<lower=0>[N] area;                 // trap area per observation
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  // [3] Species-specific initial state probabilities
  array[S] simplex[2] rho;                 // was: simplex[2] rho (shared)
  vector<lower=0>[2] rho_conc;             // Dirichlet concentration hyperprior

  vector<lower=0, upper=1>[S] theta1;      // prob stay in low state per species
  vector<lower=0, upper=1>[S] theta2;      // prob stay in high state per species

  array[S] ordered[2] log_means;           // [s][1]=low mean, [s][2]=high mean

  // Stand-level random effects for both states
  vector[F] stand_effect_raw_low;
  vector[F] stand_effect_raw_high;

  real<lower=0> sigma_low;                 // SD of stand effects, low state
  real<lower=0> sigma_high;                // SD of stand effects, high state

  // [2] Partially pooled overdispersion
  // Hyperparameters
  real mu_log_phi1;                        // mean log-overdispersion, low state
  real mu_log_phi2;                        // mean log-overdispersion, high state
  real<lower=0> sigma_log_phi1;            // SD across species, low state
  real<lower=0> sigma_log_phi2;            // SD across species, high state
  // Species-level (non-centered)
  vector[S] log_phi1_raw;                  // was: vector[S] log_phi1
  vector[S] log_phi2_raw;                  // was: vector[S] log_phi2
}

transformed parameters {
  // [2] Reconstruct phi from non-centered parameterization
  vector<lower=0>[S] phi1 = exp(mu_log_phi1 + log_phi1_raw * sigma_log_phi1);
  vector<lower=0>[S] phi2 = exp(mu_log_phi2 + log_phi2_raw * sigma_log_phi2);

  // Stand-level means for both states
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;

  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    log_alpha_low[f]  = log_means[s][1] + stand_effect_raw_low[f]  * sigma_low;
    log_alpha_high[f] = log_means[s][2] + stand_effect_raw_high[f] * sigma_high;
  }

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
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_high[f] + log_area_ratio[t], phi2[s]);
    }
  }
}

model {
  // [3] Hyperprior on Dirichlet concentration, then species-specific rho
  rho_conc ~ gamma(2, 1);                  // weakly informative, allows flexibility
  for (s in 1:S)
    rho[s] ~ dirichlet([8.0, 2.0]'); 

  // [1] Widened transition priors
  theta1 ~ beta(3, 1);                     // was beta(5, 1)
  theta2 ~ beta(1, 3);                     // was beta(1, 4)

  // Species-specific mean priors (unchanged — biologically motivated)
  log_means[1][1] ~ normal(0, 0.8);
  log_means[1][2] ~ normal(5.2, 0.6);
  
  log_means[2][1] ~ normal(2.3, 0.8);
  log_means[2][2] ~ normal(5.8, 0.8);

  log_means[3][1] ~ normal(0, log(30)/2.57);
  log_means[3][2] ~ normal(5.2, 0.5);

  log_means[4][1] ~ normal(4, 0.5);
  log_means[4][2] ~ normal(7, 0.5);

  log_means[5][1] ~ normal(5.5, 0.5);
  log_means[5][2] ~ normal(7.8, 0.5);

  // [1] Widened stand effect SD priors
  sigma_low  ~ normal(0, 1.0);             // was normal(0, 0.5)
  sigma_high ~ normal(0, 0.8);             // was normal(0, 0.8)

  stand_effect_raw_low  ~ normal(0, 1);
  stand_effect_raw_high ~ normal(0, 1);

  // [2] Hyperpriors for pooled overdispersion
  mu_log_phi1    ~ normal(log(2.0), 1.0);  // was: log_phi1 ~ normal(log(2.0), 0.5)
  mu_log_phi2    ~ normal(log(6.5), 1.0);  // was: log_phi2 ~ normal(log(6.5), 0.5)
  sigma_log_phi1 ~ normal(0, 0.5);         // cross-species SD, low state
  sigma_log_phi2 ~ normal(0, 0.5);         // cross-species SD, high state

  // [2] Non-centered species-level overdispersion
  log_phi1_raw ~ normal(0, 1);
  log_phi2_raw ~ normal(0, 1);

  // HMM marginal likelihood
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    int s        = sp[start_id];
    target += hmm_marginal(log_omega[, start_id:end_id], Gamma[s], rho[s]); // [3] rho[s]
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
      log_omega[, start_id:end_id], Gamma[s_stand], rho[s_stand]  // [3] rho[s_stand]
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
