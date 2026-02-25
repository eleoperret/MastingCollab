// 2-state HMM for seed production for one species, multi-stands with partial pooling

//before my model assumed that both means from high and low states are identidical across stands, but now I will a stand-leved random effect. 

data {
  int<lower=1> N; //number of observations
  int<lower=1> F; //number of stands
  array [N] int<lower=0> y; // seed production per year
  array [F] int<lower=1> start_idxs; // seed production per year
  array [F] int<upper=N> end_idxs; // seed production per year
}


parameters {
  simplex[2] rho;// initial state probabilities
  matrix<lower=0>[2,2] Gamma_raw; // raw transition matrix (rows to be normalized)
  ordered[2] alpha; // global state intercepts
  matrix[F,2] stand_effect_raw;// stand deviations
  vector<lower=0>[2] sigma_stand;// stand-level SD
  vector<lower=1e-6>[2] phi;// dispersion per state
}

transformed parameters {
  //vector[2] mu = exp(log_mu);modified from old code
  matrix[F,2] log_mu;
for (f in 1:F)
  for (s in 1:2)
    log_mu[f,s] = alpha[s] + stand_effect_raw[f,s] * sigma_stand[s]; //newly added
  matrix[2,N] log_omega;// log-likelihood for each year × state
  matrix[2,2] Gamma;// normalized transition matrix

  // Normalize rows of Gamma_raw to sum to 1
  for (i in 1:2)
    Gamma[i] = Gamma_raw[i] / sum(Gamma_raw[i]);
    
  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    
    for (t in start_id:end_id){
        log_omega[1,t] = neg_binomial_2_lpmf(y[t] | exp(log_mu[f,1]), phi[1]);//modfied some stuff here too
        log_omega[2,t] = neg_binomial_2_lpmf(y[t] | exp(log_mu[f,2]), phi[2]);//modfied some stuff here too
    }
  }
}

model {
  // Priors
  alpha ~ normal(4, 1); //global intercept changed from log_mu
  sigma_stand ~ exponential(1);//newly added: stand level SD
  to_vector(stand_effect_raw) ~ normal(0,1);//newly added; standardized deviations
  phi ~ gamma(2, 0.1);
  to_vector(Gamma_raw) ~ exponential(1);

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
          y_rep[t] = neg_binomial_2_rng(exp(log_mu[f,state[t]]), phi[state[t]]);
      }
  }
}




