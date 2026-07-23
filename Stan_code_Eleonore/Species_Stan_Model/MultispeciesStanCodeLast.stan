// ===============================================================
// Multispecies Hidden Markov Model (HMM) with two latent states
// ===============================================================
//
// Biological idea:
// Each species × stand combination represents a time series of
// seed production. We assume that annual seed production switches
// between two unobserved states:
//
//   State 1 = low seed production years
//   State 2 = high seed production years
//
// The observed seed counts are generated from a negative binomial
// distribution whose mean depends on the current hidden state.
//
// Model structure:
// - Two-state HMM
// - Species and stand effects are partially pooled using hierarchical(multilevel) models
// - Transition probabilities vary among species 
// - Emission means (seed production levels) vary among species and stands
// - Negative binomial emissions account for overdispersion in counts
// - State 2 is constrained to always have higher production than state 1
// ===============================================================


// ---------------------------------------------------------------
// DATA
// ---------------------------------------------------------------
data {
  int<lower=1> N; // Total number of observations across all species × stand series
  int<lower=1> F; // Number of independent time series; each series corresponds to one species × stand combination
  int<lower=1> S; // Number of species
  int<lower=1> N_stands;// number of unique stands

  array[N] int<lower=0> y;// Observed seed counts, y contains all observations concatenated across all series


  // Start and end positions of each species × stand time series, needed because all series are stored in one vector
  array[F] int<lower=1> start_idxs;
  array[F] int<upper=N> end_idxs;
  // Identity of the stand and species associated with each time series
  array[F] int<lower=1, upper=N_stands> stand_id; 
  array[F] int<lower=1, upper=S> species_id;

  //Sampling area associated with each observation
  //Used as an offset to standardize seed counts by sampling effort
  vector<lower=0>[N] area;
  
  // Indicator variable for missing observations with 0 = observed value and 1 = missing value
  array [N] int <lower= 0, upper = 1> obs_missing;
}


// ---------------------------------------------------------------
// TRANSFORMED DATA
// ---------------------------------------------------------------
// Used the smallest sampled area as reference area so that observations are expressed relative to this sampling effort
transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio; // Log-scale offset for sampling effort. This is added to the expected mean of the negative binomial model
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

// ---------------------------------------------------------------
// PARAMETERS
// ---------------------------------------------------------------
parameters {
  
// Initial state probabilities
  // Probability that a time series starts in each state
  // rho[1] = probability of starting in low state
  // rho[2] = probability of starting in high state
  simplex[2] rho;

// Transition probabilities
  // The transition matrix for each species × stand series is:
  //             next state
  //              1       2
  //
  // state 1    theta1  1-theta1
  // state 2   1-theta2 theta2
  //
  // theta1 = probability of remaining in low state
  // theta2 = probability of remaining in high state
  // Both transition probabilities are estimated hierarchically: species are partially pooled.
  real grand_logit_theta1;
  real grand_logit_theta2;

  // Species-level deviations from the global transition probability
  // Non-centered 
  vector[S]       alpha_theta1_species_nc;
  vector[S]       alpha_theta2_species_nc;
  // Amount of variation among species
  real<lower=0>   sigma_theta1_species;
  real<lower=0>   sigma_theta2_species;

// Emission means (seed production levels)
  // Each state has its own expected seed production:
  // State 1:
  // log(mean) = global low mean
  //             + species effect
  //             + stand effect
  // State 2:
  // log(mean) = low-state mean + positive difference (delta)
  // The delta formulation forces state 2 to always have higher seed production than state 1.

  real             mu_log_low; // Global mean of the low production state
  vector[S]        alpha_low_species;// Species deviations in low-state production
  real<lower=0>    sigma_low_species;// Variation among species
  matrix [S, N_stands] alpha_stand; // Stand deviations in low-state production
  real <lower=0> sigma_stand;// Variation among stands

  real             mu_log_delta;// Difference between low and high states
  vector[S]        log_delta_species_nc;// Species-specific differences non-centered
  real<lower=0>    sigma_log_delta_species;// Variation among species
  matrix [S, N_stands] log_delta_stand_nc; // Stand-specific differences non-centered
  real <lower=0> sigma_delta_stand;// Variation among stands

// Negative binomial dispersion parameters
  vector [S] log_phi_species;// Species-specific overdispersion
  vector [2] log_phi_state;// State-specific overdispersion
}


// ---------------------------------------------------------------
// TRANSFORMED PARAMETERS
// ---------------------------------------------------------------
// This block converts hierarchical parameters into the actual species- and stand-specific parameters used in the likelihood.
//
// Non-centered parameterizations are transformed here:
// actual effect = population SD × standardized deviation

transformed parameters {
  
  // Transition probability hierarchical effects
  vector[S]        alpha_theta1_species = sigma_theta1_species * alpha_theta1_species_nc;
  vector[S]        alpha_theta2_species = sigma_theta2_species * alpha_theta2_species_nc;

  // State difference (high vs low production)
  vector[S] log_delta_species    = sigma_log_delta_species * log_delta_species_nc; // Species-specific deviation in the difference between states
  matrix[S, N_stands] log_delta_stand   = sigma_delta_stand * log_delta_stand_nc;  // Stand-specific deviation in the difference between states
  
  // Species × stand specific emission means 
  // These are the expected seed production levels for each time series before accounting for sampling area.
  // The means are stored on the log scale because the negative binomial distribution uses a positive mean parameter.
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;

  for (f in 1:F) {
    int s  = species_id[f]; // Species identity of this time series
    int st = stand_id[f];// Stand identity of this time series
    // Low production state
    // log(mean seeds) =overall low-state mean+ species effect+ stand effect
    log_alpha_low[f]  = mu_log_low
                      + alpha_low_species[s]
                      + alpha_stand[s,st];
    // High production state
    // The difference between states is constrained to be positive:log(mean high) =log(mean low)+ positive delta
    // log1p_exp() transforms any value into a positive quantity.
    // This prevents label switching where the HMM swaps the meaning of "low" and "high" states.
    log_alpha_high[f] = log_alpha_low[f]
                  + log1p_exp(mu_log_delta)
                  + log1p_exp(log_delta_species[s])
                  + log1p_exp(log_delta_stand[s,st]);
  }

  // Species × stand specific transition matrices
  // Each series has its own transition matrix:
  //        Low       High
  //
  // Low    theta1    1-theta1
  // High   1-theta2  theta2
  //
  // theta values are estimated on the logit scale and transformed to probabilities using the inverse-logit function.
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;  // One transition matrix per species × stand series

  for (f in 1:F) {
    int s  = species_id[f];
    // Probability of staying in low state
    theta1[f] = inv_logit(grand_logit_theta1
                         + alpha_theta1_species[s]);
    // Probability of staying in high state
    theta2[f] = inv_logit(grand_logit_theta2
                         + alpha_theta2_species[s]);
    // Construct transition matrix
    Gamma[f][1, 1] = theta1[f];
    Gamma[f][1, 2] = 1 - theta1[f];
    Gamma[f][2, 1] = 1 - theta2[f];
    Gamma[f][2, 2] = theta2[f];
  }

  // Emission probabilities
  // log_omega contains the probability of observing each count under each hidden state.
  // Row 1 = probability under low state
  // Row 2 = probability under high state
  // These likelihoods are later used by hmm_marginal() to integrate over all possible hidden state sequences.
  matrix[2, N] log_omega;
  for (f in 1:F) {
    int s = species_id[f];  
    for (t in start_idxs[f]:end_idxs[f]) {
      if (obs_missing[t]){ // Missing observations provide no information about emissions.The HMM can still estimate the hidden state using transitions.
        log_omega[1,t]= 0;
        log_omega [2,t] =0;
      } else {
        log_omega[1, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1])); // Probability of observed count assuming low state
      log_omega[2, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2]));// Probability of observed count assuming high state
      }
    }
  }
}



// ---------------------------------------------------------------
// Model
// ---------------------------------------------------------------

model {
  // Initial state distribution
  rho ~ dirichlet(rep_vector(2.0, 2));//The HMM is allowed to start in either state, but extreme probabilities are discouraged.

  // Transition probability priors
  grand_logit_theta1 ~ normal(0.5,   1); // Global tendency to remain in each state
  grand_logit_theta2 ~ normal(-1,   0.5); // Global tendency to remain in each state

  alpha_theta1_species_nc ~ normal(0, 1); // Species-level variation
  alpha_theta2_species_nc ~ normal(0, 1);
  sigma_theta1_species    ~ normal(0, 0.7);
  sigma_theta2_species    ~ normal(0, 0.7);

  // Emission mean priors
  mu_log_low           ~ normal(2, 1.5);// Average production in low state
  alpha_low_species ~ normal(0, sigma_low_species);// Species variation in low-state production
  sigma_low_species    ~ normal(0.5, 1);
  for (s in 1:S) 
    alpha_stand [s]   ~ normal(0, sigma_stand);// Stand variation in low-state production
  sigma_stand    ~ normal(1, 0.5);

  mu_log_delta           ~ normal(1.5, 1.5); // Difference between low and high states
  log_delta_species_nc   ~ normal(0, 1);
  sigma_log_delta_species ~ normal(1.5, 0.7);
  for (s in 1:S)
    log_delta_stand_nc [s] ~ normal(0, 1);
  sigma_delta_stand ~ normal(1, 0.5); 

  
  // Negative binomial dispersion priors
  log_phi_species ~ normal(log(4), 0.6); 
  log_phi_state ~ normal(0, 0.3); 

  //HMM likelihood
  // hmm_marginal() computes the likelihood while integrating over all possible hidden state sequences.
  // This avoids having to estimate the state sequence directly.
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]],
        Gamma[f],
        rho);
}


// ---------------------------------------------------------------
// generated quantities
// ---------------------------------------------------------------

 // Posterior predictive observations Used for posterior predictive checks:
 // Does the fitted model generate data similar to the observed data?

generated quantities {
  array[N] int<lower=0>          y_rep;// Estimated latent state sequence
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];
    // Sample one possible hidden state trajectory
    // Uses the posterior transition probabilities and emission probabilities to reconstruct the most likely sequence of
    // hidden states. This allows extraction of masting events and state duration.
    state[start_id:end_id] = hmm_latent_rng(
        log_omega[, start_id:end_id],
        Gamma[f],
        rho);

    for (t in start_id:end_id) {
  int s = species_id[f];  
  
  if (obs_missing[t]) {// Missing observations remain missing
    y_rep[t] = 0; 
  } else if (state[t] == 1) {
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1]));// Generate replicated count assuming low state
  }else{
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2])); // Generate replicated count assuming high state
      }
    }
  }
}


