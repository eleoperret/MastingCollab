// 2-state HMM for seed production for one species, multi-stands

data {
  int<lower=1> N; //number of observations
  int<lower=1> F; //number of stands
  array [N] int<lower=0> y; // seed production per year
  array [F] int<lower=1> start_idxs; // seed production per year
  array [F] int<upper=N> end_idxs; // seed production per year
}

//ragged array, vector
parameters {
  simplex[2] rho;            // initial state probabilities
  matrix<lower=0>[2,2] Gamma_raw; // raw transition matrix (rows to be normalized)
  ordered[2] log_mu;         // log mean seed production for low/high states
  vector<lower=1e-6>[2] phi; // dispersion per state
}

transformed parameters {
  vector[2] mu = exp(log_mu);       // natural-scale means
  matrix[2,N] log_omega;            // log-likelihood for each year × state
  matrix[2,2] Gamma;                // normalized transition matrix

  // Normalize rows of Gamma_raw to sum to 1
  for (i in 1:2)
    Gamma[i] = Gamma_raw[i] / sum(Gamma_raw[i]);
    
  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    
    for (t in start_id:end_id){
        log_omega[1,t] = neg_binomial_2_lpmf(y[t] | mu[1], phi[1]);
        log_omega[2,t] = neg_binomial_2_lpmf(y[t] | mu[2], phi[2]);
    }
  }
}


model {
  // Priors
  log_mu ~ normal(4, 1);        // allows seed production roughly up to ~400
  phi ~ gamma(2, 0.1);          // moderate overdispersion
  to_vector(Gamma_raw) ~ exponential(1); // weakly informative positive prior
  //theta~beta (1,1);

  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    target += hmm_marginal(log_omega[,start_id:end_id], Gamma, rho);
  }
}

generated quantities {
  // Use the same ragged structure as in data
  array [N] int<lower=0> y_rep;  // max(T_i) because Stan needs a fixed upper bound
  array [N] int<lower=1, upper= 2> state;
  
   for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    state[start_id:end_id] = hmm_latent_rng(log_omega[,start_id:end_id], Gamma, rho);
    
    for (t in start_id:end_id){
          y_rep[t] = neg_binomial_2_rng(mu[state[1]], phi[state[1]]);
      }
  }
}



