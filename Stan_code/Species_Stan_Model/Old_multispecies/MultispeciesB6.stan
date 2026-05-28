// 2-state HMM for seed production — multispecies model 

data {
  int<lower=1> N;//number of rows for whole dataset
  int<lower=1> F;// number of rows for species x stand combination
  int<lower=1> S;//number of species

  array[N] int<lower=0> y; //observation of seed counts
  array[N] int<lower=1, upper=S> sp; //unique species

  array[F] int<lower=1> start_idxs; //start for specie x stand observation
  array[F] int<upper=N> end_idxs; //end for specie x stand observation

  vector<lower=0>[N] area; // seed trap offset
}

//offset for amount of traps per years
transformed data {
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}

parameters {
  // Initial state probability
  array[S] simplex[2] rho;

  // Transtion probability 
  vector<lower=0, upper=1>[S] theta1; // prob stay in low state, per species
  vector<lower=0, upper=1>[S] theta2; // prob stay in high state, per species

  // Species mean seed production
  array[S] ordered[2] log_means; 

  //Partial pooling (for high and low state seed production) : stand-level log-means
  vector[F] log_alpha_low;              
  vector[F] log_alpha_high;               

  //Stand effect on seed production
  real<lower=0> sigma_low;
  real<lower=0> sigma_high;

  //Overdispersion — centered, partially pooled across species
  real mu_log_phi1;                        
  real mu_log_phi2;                        
  real<lower=0> sigma_log_phi1;            
  real<lower=0> sigma_log_phi2;           
  vector[S] z_phi1;                      
  vector[S] z_phi2;                      
}

transformed parameters {
  //Now my phi is non-centered after diagnostic plots.
  vector [S] log_phi1 = mu_log_phi1 + sigma_log_phi1 * z_phi1;
  vector [S] log_phi2 = mu_log_phi2 + sigma_log_phi2 * z_phi2;
  
  //Overdispersion on natural scale
  vector<lower=0>[S] phi1 = exp(log_phi1);
  vector<lower=0>[S] phi2 = exp(log_phi2);

//To change?
  //Transition matrices per species
  array[S] matrix[2, 2] Gamma;
  for (f in 1:F) {
    Gamma[s][1, 1] = theta1[s];
    Gamma[s][1, 2] = 1 - theta1[s];
    Gamma[s][2, 1] = 1 - theta2[s];
    Gamma[s][2, 2] = theta2[s];
  }

  // Emission log-likelihoods
  // log_alpha[f]  — stand-level mean
  // phi1[s] — overdispersion, low state (varies by species)
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
  // Species-specific initial state probabilities
  for (s in 1:S)
    rho[s] ~ dirichlet([8.0, 2.0]');// most stands start in low state

  //Stand-level transition probabilities — with a fixed Beta prior
  theta1 ~ beta(4, 1);// tends to stay in low state
  theta2 ~ beta(1, 4);// tends to leave high state

  // Species-level mean priors (biologically motivated)
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
  
  // Stand effect prior
  sigma_low  ~ normal(0, 1.0);
  sigma_high ~ normal(0, 1.0);

  // Stand-level means
  // log_alpha (partially pooled for stands)
  for (f in 1:F) {
    int s = sp[start_idxs[f]];
    log_alpha_low[f]  ~ normal(log_means[s][1], sigma_low);
    log_alpha_high[f] ~ normal(log_means[s][2], sigma_high);
  }

  //Overdispersion hyperpriors 
  mu_log_phi1    ~ normal(log(2.0), 1.0);
  mu_log_phi2    ~ normal(log(6.5), 1.0);
  sigma_log_phi1 ~ normal(0, 0.5);
  sigma_log_phi2 ~ normal(0, 0.5);

  //Centered species-level overdispersion
  z_phi1 ~ std_normal();
  z_phi2 ~ std_normal();
  //log_phi1 ~ normal(mu_log_phi1, sigma_log_phi1);
  //log_phi2 ~ normal(mu_log_phi2, sigma_log_phi2);

  // HMM marginal likelihood
  // Gamma[f] — transition matrix per stand
  // rho[s]   — initial state probability per species
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
