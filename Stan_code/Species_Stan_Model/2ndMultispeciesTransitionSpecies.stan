// 2-state HMM for seed production - multispecies model
// State 1 = low seed production
// State 2 = high seed production (masting)
// Ordered means to prevent label switching
// Stand-level random effects for high state
// Species-specific transitions and emission parameters
// PSME removed due to poor state separation (near-unimodal distribution)
// log_area_ratio precomputed in transformed data for efficiency

data {
  int<lower=1> N;                          // total number of observations
  int<lower=1> F;                          // number of stands
  int<lower=1> S;                          // number of species (4 without PSME)
  
  array[N] int<lower=0> y;                 // seed counts
  array[N] int<lower=1, upper=S> sp;       // species index per observation
  
  array[F] int<lower=1> start_idxs;        // start index per stand
  array[F] int<upper=N> end_idxs;          // end index per stand
  
  vector<lower=0>[N] area;                 // trap area per observation
}

transformed data {
  real baseline_area = min(area);
  
  // Precomputed once — avoids recomputing at every leapfrog step
  vector[N] log_area_ratio;
  for (t in 1:N) {
    log_area_ratio[t] = log(area[t] / baseline_area);
  }
}

parameters {
  simplex[2] rho;                          // initial state probabilities
  
  // Species-specific transition probabilities
  vector<lower=0, upper=1>[S] theta1;      // prob of staying in low state
  vector<lower=0, upper=1>[S] theta2;      // prob of staying in high state
  
  // Species-specific emission means
  // log_means[s][1] = low  state mean
  // log_means[s][2] = high state mean
  // ordered enforces [1] < [2] — prevents label switching
  array[S] ordered[2] log_means;
  
  // Stand-level random effects on high state only
  vector[F] stand_effect_raw;              // non-centered parameterisation
  real<lower=0> sigma;                     // stand-level SD
  
  // Dispersion on log scale — ensures positivity
  vector[S] log_phi1;                      // dispersion for low  state
  vector[S] log_phi2;                      // dispersion for high state
}

transformed parameters {
  // Back-transform dispersion to natural scale
  vector<lower=0>[S] phi1 = exp(log_phi1);
  vector<lower=0>[S] phi2 = exp(log_phi2);
  
  // Stand-level high state mean
  vector[F] log_alpha;
  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    log_alpha[f] = log_means[s][2] + stand_effect_raw[f] * sigma;
  }
  
  // Species-specific transition matrices
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
      log_omega[1, t] = neg_binomial_2_log_lpmf(y[t] | log_means[s][1] + log_area_ratio[t], phi1[s]);
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha[f]    + log_area_ratio[t], phi2[s]);
    }
  }
}

model {
  // --- Priors ---
  rho ~ dirichlet(rep_vector(1.0, 2));
  
  theta1 ~ beta(3, 1);
  theta2 ~ beta(2, 2);
  
  // Species-specific informed priors based on data inspection
  // ABAM (s=1): clear masting, near-zero low state
  log_means[1][1] ~ normal(0.5, 1.0);
  log_means[1][2] ~ normal(3.0, 1.0);
  
  // CANO (s=2): strongest masting signal, very low baseline
  log_means[2][1] ~ normal(0.5, 1.0);
  log_means[2][2] ~ normal(3.5, 1.0);
  
  // THPL (s=3): moderate baseline, high masting peaks
  log_means[3][1] ~ normal(2.0, 1.0);
  log_means[3][2] ~ normal(5.0, 1.0);
  
  // TSHE (s=4): high baseline, very high masting peaks
  log_means[4][1] ~ normal(4.5, 0.5);   // tighter — states close together
  log_means[4][2] ~ normal(5.5, 0.5);   // tighter — states close together
  
  // Stand random effects
  sigma            ~ normal(0, 0.5);
  stand_effect_raw ~ normal(0, 1);
  
  // Dispersion
  log_phi1 ~ normal(log(6.5), 0.5);
  log_phi2 ~ normal(log(6.5), 0.5);
  
  // --- Likelihood ---
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
        y_rep[t] = neg_binomial_2_log_rng(log_means[s][1] + log_area_ratio[t], phi1[s]);
      } else {
        y_rep[t] = neg_binomial_2_log_rng(log_alpha[f]    + log_area_ratio[t], phi2[s]);
      }
    }
  }
}
