// =============================================================================
// 2-state HMM for seed production — multispecies model v3
// =============================================================================
// CHANGES FROM v1 (original model):
//
// [1] WIDER PRIORS throughout
//     - sigma_low  : normal(0, 0.5) → normal(0, 1.0)
//     - sigma_high : normal(0, 0.8) → normal(0, 1.0)
//     - log_phi    : sd 0.5 → 1.0
//     - theta      : beta(5,1)/beta(1,4) → now hierarchical (see [2])
//
// [2] STAND-LEVEL TRANSITION PROBABILITIES, pooled within species
//     - Previously: theta1[s], theta2[s] — one per species
//     - Now:        theta1[f], theta2[f] — one per stand, pooled within species
//     - Pooling via Beta hyperprior:
//         theta1[f] ~ beta(mu_theta1[s] * kappa_theta1,
//                         (1 - mu_theta1[s]) * kappa_theta1)
//         theta2[f] ~ beta(mu_theta2[s] * kappa_theta2,
//                         (1 - mu_theta2[s]) * kappa_theta2)
//     - mu_theta1[s], mu_theta2[s]: species-level mean transition prob
//     - kappa_theta1, kappa_theta2: concentration (shared across species)
//       controls how tightly stands cluster around species mean
//     - Biologically: all ABAM stands tend to stay in low state, but
//       individual stands can deviate based on local conditions
//
// [3] PARTIALLY POOLED OVERDISPERSION across species (from v2)
//     - log_phi1[s] ~ normal(mu_log_phi1, sigma_log_phi1)
//     - log_phi2[s] ~ normal(mu_log_phi2, sigma_log_phi2)
//     - No stand-level phi — already have stand-level mean variation
//
// [4] CENTERED parameterization for stand effects on the mean
//     - Previously (non-centered): log_alpha[f] = log_means[s] + raw[f] * sigma
//     - Now (centered):            log_alpha[f] ~ normal(log_means[s], sigma)
//     - Works better when stands have enough data (~13 years each)
//
// [5] SPECIES-SPECIFIC rho with fixed Dirichlet prior
//     - rho[s]: each species has its own initial state probability
//     - Fixed concentration [8, 2] — most stands start in low state
//     - No estimated rho_conc hyperparameter (not identifiable with S=5)
//
// NOT CHANGED:
//     - log_means priors remain species-specific and biologically motivated
//       (species are not exchangeable in mean seed production)
//     - Stand-level phi not added (too many parameters for 11-15 obs/stand)
//
// TODO (lower priority):
// [6] Correlated stand effects across states (bivariate normal / Cholesky)
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
  // [5] Species-specific initial state probabilities (fixed concentration)
  array[S] simplex[2] rho;

  // [2] Stand-level transition probabilities
  vector<lower=0, upper=1>[F] theta1;      // prob stay in low state, per stand
  vector<lower=0, upper=1>[F] theta2;      // prob stay in high state, per stand

  // [2] Species-level hyperparameters for transition probabilities
  vector<lower=0, upper=1>[S] mu_theta1;   // species mean, low-state persistence
  vector<lower=0, upper=1>[S] mu_theta2;   // species mean, high-state persistence
  real<lower=0> kappa_theta1;              // concentration, low state (shared)
  real<lower=0> kappa_theta2;              // concentration, high state (shared)

  // Species-level mean seed production
  array[S] ordered[2] log_means;           // [s][1]=low mean, [s][2]=high mean

  // [4] Centered stand effects — directly the log-mean per stand per state
  vector[F] log_alpha_low;                 // centered: ~ normal(log_means[s][1], sigma_low)
  vector[F] log_alpha_high;                // centered: ~ normal(log_means[s][2], sigma_high)

  // [1] Stand effect SDs — wider priors
  real<lower=0> sigma_low;
  real<lower=0> sigma_high;

  // [3] Overdispersion — partially pooled across species
  real mu_log_phi1;                        // hyperprior mean, low state
  real mu_log_phi2;                        // hyperprior mean, high state
  real<lower=0> sigma_log_phi1;            // hyperprior SD, low state
  real<lower=0> sigma_log_phi2;            // hyperprior SD, high state
  vector[S] log_phi1_raw;                  // non-centered species deviations
  vector[S] log_phi2_raw;
}

transformed parameters {
  // [3] Reconstruct phi
  vector<lower=0>[S] phi1 = exp(mu_log_phi1 + log_phi1_raw * sigma_log_phi1);
  vector<lower=0>[S] phi2 = exp(mu_log_phi2 + log_phi2_raw * sigma_log_phi2);

  // [2] Transition matrices — now per stand f, not per species s
  array[F] matrix[2, 2] Gamma;
  for (f in 1:F) {
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
      int s = sp[t];
      log_omega[1, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_low[f]  + log_area_ratio[t], phi1[s]);
      log_omega[2, t] = neg_binomial_2_log_lpmf(y[t] | log_alpha_high[f] + log_area_ratio[t], phi2[s]);
    }
  }
}

model {
  // [5] Species-specific initial state probabilities
  for (s in 1:S)
    rho[s] ~ dirichlet([8.0, 2.0]');       // fixed: most stands start in low state

  // [2] Hyperpriors for species-level transition probability means
  mu_theta1 ~ beta(4, 1);                  // species tend to stay in low state
  mu_theta2 ~ beta(1, 4);                  // species tend to leave high state quickly
  kappa_theta1 ~ gamma(4, 0.5);            // moderate concentration — allows stand variation
  kappa_theta2 ~ gamma(4, 0.5);

  // [2] Stand-level transition probabilities pooled within species
  // Beta reparameterized via mean (mu) and concentration (kappa):
  //   alpha = mu * kappa,  beta = (1 - mu) * kappa
  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    theta1[f] ~ beta(mu_theta1[s] * kappa_theta1,
                     (1 - mu_theta1[s]) * kappa_theta1);
    theta2[f] ~ beta(mu_theta2[s] * kappa_theta2,
                     (1 - mu_theta2[s]) * kappa_theta2);
  }

  // Species-level mean priors (biologically motivated — unchanged)
  log_means[1][1] ~ normal(0, 0.8);
  log_means[1][2] ~ normal(5.2, 0.6);

  log_means[2][1] ~ normal(2.3, 0.8);
  log_means[2][2] ~ normal(5.8, 0.8);

  log_means[3][1] ~ normal(0, log(30)/2.57);
  log_means[3][2] ~ normal(5.2, 0.5);

  log_means[4][1] ~ normal(4, 0.5);
  log_means[4][2] ~ normal(7, 0.5);

  log_means[5][1] ~ normal(5.5, 0.5);     // updated for TSHE based on data
  log_means[5][2] ~ normal(7.8, 0.5);     // updated for TSHE based on data

  // [1] Wider priors on stand effect SDs
  sigma_low  ~ normal(0, 1.0);            // was normal(0, 0.5)
  sigma_high ~ normal(0, 1.0);            // was normal(0, 0.8)

  // [4] Centered parameterization for stand-level means
  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    log_alpha_low[f]  ~ normal(log_means[s][1], sigma_low);
    log_alpha_high[f] ~ normal(log_means[s][2], sigma_high);
  }

  // [3] Hyperpriors for overdispersion
  mu_log_phi1    ~ normal(log(2.0), 1.0);
  mu_log_phi2    ~ normal(log(6.5), 1.0);
  sigma_log_phi1 ~ normal(0, 0.5);
  sigma_log_phi2 ~ normal(0, 0.5);

  // [3] Non-centered species-level overdispersion deviations
  log_phi1_raw ~ normal(0, 1);
  log_phi2_raw ~ normal(0, 1);

  // HMM marginal likelihood — Gamma now indexed by f (stand), rho by s (species)
  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    int s        = sp[start_id];
    target += hmm_marginal(log_omega[, start_id:end_id], Gamma[f], rho[s]);
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
      log_omega[, start_id:end_id], Gamma[f], rho[s_stand]
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
