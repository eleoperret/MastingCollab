//HMM 2 state - Single species
//Transition matrix partially pooled by stand (centered)
//Distribution: two NB for each state partially pooled by stand
//Overdispersion per state
//Forced state 2 above state 1 with delta
//new based on victor's comment in #66
//Sigma_low_stand also to 0.5 instead of 1
// (comes from adter the code : Single specie----Victor.stan) PLus de-centering for the theta and the high distribution
//added now the de-centering on the low emission too.
//changed from the preivous model priors : SinglespeciesPartialPooling2_delta_NCbisFULL.stan


///Prior testing
// A ) based on git hub comment #67
//mu_log_low ~ normal(2.8, 1) to normal (2.8, 0.5) : seeds from 2-122 to 6-45
//now trying to reduce the theta 2 probability from (-0.5, 1) to (-0.5, 0.5): probability of still 37% with SD (18-62%); could make lower

data {
  int<lower=1> N;
  int<lower=1> F;
  int<lower=1> N_stands;

  array[N] int<lower=0> y;

  array[F] int<lower=1> start_idxs;
  array[F] int<upper=N> end_idxs;
  array[F] int<lower=1, upper=N_stands> stand_id;

  vector<lower=0>[N] area;
  
  array [N] int<lower=0, upper = 1> obs_missing;
}

transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  simplex[2] rho;

  // Transitions
  real grand_logit_theta1;
  real grand_logit_theta2;

  vector[N_stands] alpha_theta1_stand_nc;   
  vector[N_stands] alpha_theta2_stand_nc;  
  real<lower=0> sigma_theta1_stand;
  real<lower=0> sigma_theta2_stand;

  // Emission means; 
  real mu_log_low;
  real mu_log_delta;
  

  // Stand random effects — centered
  vector [N_stands] alpha_low_stand_nc;
  vector[N_stands] log_delta_tilde;
  real<lower=0> sigma_low_stand;
  real<lower=0> sigma_log_delta;

  // Dispersion per state
  vector[2] log_phi_state;
}

transformed parameters {
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;
  
  vector[N_stands] alpha_theta1_stand = sigma_theta1_stand * alpha_theta1_stand_nc;
  vector[N_stands] alpha_theta2_stand = sigma_theta2_stand * alpha_theta2_stand_nc;
  vector[N_stands] log_low_stand = mu_log_low + sigma_low_stand * alpha_low_stand_nc;
  vector[N_stands] log_1p_delta  = log1p_exp(mu_log_delta + sigma_log_delta * log_delta_tilde);
  
  for (f in 1:F) {
    int st = stand_id[f];

    log_alpha_low [f] = log_low_stand[st];
    log_alpha_high[f] = log_low_stand[st] + log_1p_delta[st];

    
    theta1[f] = inv_logit(grand_logit_theta1 + alpha_theta1_stand[st]);
    theta2[f] = inv_logit(grand_logit_theta2 + alpha_theta2_stand[st]);

    Gamma[f][1, 1] = theta1[f];
    Gamma[f][1, 2] = 1 - theta1[f];
    Gamma[f][2, 1] = 1 - theta2[f];
    Gamma[f][2, 2] = theta2[f];
  }

  matrix[2, N] log_omega;
  for (f in 1:F) {
    for (t in start_idxs[f]:end_idxs[f]) {
      if (obs_missing[t]){
        log_omega[1, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_state[1]));
        log_omega[2, t] = neg_binomial_2_log_lpmf(
          y[t] | log_alpha_high[f] + log_area_ratio[t], exp(log_phi_state[2]));
      } else {
        log_omega [1,t] = 0;
        log_omega[2,t] =0;
      }
    }
  }
}

model {
  rho ~ dirichlet(rep_vector(2.0, 2));

  // Transitions
  grand_logit_theta1 ~ normal(0.5, 1);
  grand_logit_theta2 ~ normal(-1, 0.5);

  alpha_theta1_stand_nc ~ normal(0, 1);  
  alpha_theta2_stand_nc ~ normal(0, 1);   
  sigma_theta1_stand ~ normal(1, 0.3);
  sigma_theta2_stand ~ normal(1, 0.5);

  // Emission means
  mu_log_low                ~ normal(2, 1.5);
  alpha_low_stand_nc        ~ normal(0, 1);
  sigma_low_stand           ~ normal(1, 0.5);

  mu_log_delta              ~ normal(1.5, 1.5);
  log_delta_tilde           ~ normal(0, 1); 
  sigma_log_delta           ~ normal (1.5,0.7);

  // Dispersion
  log_phi_state ~ normal(log(4), 0.6); 

  // Likelihood
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]], Gamma[f], rho);
}

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];

    state[start_id:end_id] = hmm_latent_rng(
        log_omega[, start_id:end_id], Gamma[f], rho);

    for (t in start_id:end_id) {
      if (state[t] == 1)
        y_rep[t] = neg_binomial_2_log_rng(
            log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_state[1]));
      else
        y_rep[t] = neg_binomial_2_log_rng(
            log_alpha_high[f] + log_area_ratio[t], exp(log_phi_state[2]));
    }
  }
}

