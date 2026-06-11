// HMM 2 state - Multi species
// Species-level baselines for emission means (key fix for state collapse)
// Transitions: additive species + stand effects
// Delta: additive species + stand effects outside log1p_exp (with strong prior on mu_log_delta)
// Overdispersion: additive per species + per state
// Partial pooling non-centered throughout

data {
  int<lower=1> N;
  int<lower=1> F;
  int<lower=1> S;
  int<lower=1> N_stands;

  array[N] int<lower=0> y;

  array[F] int<lower=1>                 start_idxs;
  array[F] int<upper=N>                 end_idxs;
  array[F] int<lower=1, upper=N_stands> stand_id;
  array[F] int<lower=1, upper=S>        species_id;  // series-level, length F

  vector<lower=0>[N] area;
  array[N] int<lower=0, upper=1> obs_missing;
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  simplex[2] rho;

  // ── Transitions ────────────────────────────────────────────────────────────
  real grand_logit_theta1;
  real grand_logit_theta2;

  vector[S]        alpha_theta1_sp_nc;
  vector[S]        alpha_theta2_sp_nc;
  real<lower=0>    sigma_theta1_sp;
  real<lower=0>    sigma_theta2_sp;

  vector[N_stands] alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand_nc;
  real<lower=0>    sigma_theta1_stand;
  real<lower=0>    sigma_theta2_stand;

  // ── Emission means ─────────────────────────────────────────────────────────
  real             grand_mu_log_low;
  real<lower=0>    sigma_mu_log_low;
  vector[S]        mu_log_low_sp_nc;

  vector[N_stands] alpha_low_stand_nc;
  real<lower=0>    sigma_low_stand;

  real             grand_mu_log_delta;
  vector[S]        log_delta_sp_nc;
  real<lower=0>    sigma_log_delta_sp;
  vector[N_stands] log_delta_stand_nc;
  real<lower=0>    sigma_log_delta_stand;

  // ── Overdispersion ─────────────────────────────────────────────────────────
  real             mu_log_phi_sp;
  real<lower=0>    sigma_log_phi_sp;
  vector[S]        log_phi_sp_nc;
  vector[2]        log_phi_state;
}

transformed parameters {
  // Species baselines
  vector[S] mu_log_low_sp = grand_mu_log_low
                           + sigma_mu_log_low * mu_log_low_sp_nc;

  // Stand effects on low state (non-centered)
  vector[N_stands] alpha_low_stand    = sigma_low_stand       * alpha_low_stand_nc;

  // Delta effects (non-centered)
  vector[S]        log_delta_sp       = sigma_log_delta_sp    * log_delta_sp_nc;
  vector[N_stands] log_delta_stand    = sigma_log_delta_stand * log_delta_stand_nc;

  // Transition effects (non-centered)
  vector[S]        alpha_theta1_sp    = sigma_theta1_sp    * alpha_theta1_sp_nc;
  vector[S]        alpha_theta2_sp    = sigma_theta2_sp    * alpha_theta2_sp_nc;
  vector[N_stands] alpha_theta1_stand = sigma_theta1_stand * alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand = sigma_theta2_stand * alpha_theta2_stand_nc;

  // Overdispersion per species (non-centered)
  vector[S] log_phi_sp = mu_log_phi_sp + sigma_log_phi_sp * log_phi_sp_nc;

  // Per-series quantities
  vector[F] log_alpha_low_f;
  vector[F] log_alpha_high_f;
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;

  for (f in 1:F) {
    int s  = species_id[f];
    int st = stand_id[f];

    log_alpha_low_f[f]  = mu_log_low_sp[s] + alpha_low_stand[st];

    log_alpha_high_f[f] = log_alpha_low_f[f]
                        + grand_mu_log_delta
                        + log_delta_sp[s]
                        + log_delta_stand[st];

    theta1[f] = inv_logit(grand_logit_theta1
                         + alpha_theta1_sp[s]
                         + alpha_theta1_stand[st]);
    theta2[f] = inv_logit(grand_logit_theta2
                         + alpha_theta2_sp[s]
                         + alpha_theta2_stand[st]);

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
      if (obs_missing[t]) {
        log_omega[1, t] = 0;
        log_omega[2, t] = 0;
      } else {
        log_omega[1, t] = neg_binomial_2_log_lpmf(
            y[t] | log_alpha_low_f[f]  + log_area_ratio[t],
            exp(log_phi_sp[s] + log_phi_state[1]));
        log_omega[2, t] = neg_binomial_2_log_lpmf(
            y[t] | log_alpha_high_f[f] + log_area_ratio[t],
            exp(log_phi_sp[s] + log_phi_state[2]));
      }
    }
  }
}

model {
  rho ~ dirichlet(rep_vector(2.0, 2));

  // Transitions
  grand_logit_theta1 ~ normal(0.5, 1);
  grand_logit_theta2 ~ normal(-0.5, 1);

  alpha_theta1_sp_nc    ~ normal(0, 1);
  alpha_theta2_sp_nc    ~ normal(0, 1);
  sigma_theta1_sp       ~ normal(0, 0.7);
  sigma_theta2_sp       ~ normal(0, 0.7);

  alpha_theta1_stand_nc ~ normal(0, 1);
  alpha_theta2_stand_nc ~ normal(0, 1);
  sigma_theta1_stand    ~ normal(0, 0.7);
  sigma_theta2_stand    ~ normal(0, 0.7);

  // Emission means
  grand_mu_log_low   ~ normal(2.8, 0.5);
  sigma_mu_log_low   ~ normal(0, 1);
  mu_log_low_sp_nc   ~ normal(0, 1);

  alpha_low_stand_nc ~ normal(0, 1);
  sigma_low_stand    ~ normal(0, 0.5);

  grand_mu_log_delta    ~ normal(1.5, 0.5);
  log_delta_sp_nc       ~ normal(0, 1);
  sigma_log_delta_sp    ~ normal(0, 0.5);
  log_delta_stand_nc    ~ normal(0, 1);
  sigma_log_delta_stand ~ normal(0, 0.5);

  // Overdispersion
  mu_log_phi_sp    ~ normal(log(4), 0.6);
  sigma_log_phi_sp ~ normal(0, 0.5);
  log_phi_sp_nc    ~ normal(0, 1);
  log_phi_state    ~ normal(0, 0.3);

  // Likelihood
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]],
        Gamma[f],
        rho);
}

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F) {
    int s        = species_id[f];
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];

    state[start_id:end_id] = hmm_latent_rng(
        log_omega[, start_id:end_id], Gamma[f], rho);

    for (t in start_id:end_id) {
      if (obs_missing[t]) {
        y_rep[t] = 0;
      } else if (state[t] == 1) {
        y_rep[t] = neg_binomial_2_log_rng(
            log_alpha_low_f[f]  + log_area_ratio[t],
            exp(log_phi_sp[s] + log_phi_state[1]));
      } else {
        y_rep[t] = neg_binomial_2_log_rng(
            log_alpha_high_f[f] + log_area_ratio[t],
            exp(log_phi_sp[s] + log_phi_state[2]));
      }
    }
  }
}




