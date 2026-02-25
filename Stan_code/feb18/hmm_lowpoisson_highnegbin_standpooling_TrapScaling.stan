// 2-state HMM for seed production for one species, multi-stands with partial pooling

//before my model assumed that both means from high and low states are identidical across stands, but now I will a stand-leved random effect. 

data {
  int<lower=1> N; //number of observations
  int<lower=1> F; //number of stands
  array [N] int<lower=0> y; // seed production per year
  array [F] int<lower=1> start_idxs; // seed production per year
  array [F] int<upper=N> end_idxs; // seed production per year
  vector<lower=0>[N] area;//or array[N] real area int<lower=0> area ?;
  
}



parameters {
  simplex[2] rho;// initial state probabilities
  matrix<lower=0>[2,2] Gamma_raw; // raw transition matrix (rows to be normalized)//changer ca
  real log_lambda; // average seed production in low years
  //real<lower=log_lambda> log_mu;// replaced because issues as parameter from parameter
  real log_mu_raw;
  
  vector[F] stand_effect_raw;// stand deviations
  real<lower=0> sigma;// stand-level SD
  real<lower=0> phi;// dispersion per state
}

transformed parameters {
  real log_mu;
  log_mu = log_lambda + exp(log_mu_raw); // ensures log_mu > log_lambda

  // stand-level random effect for mast years
  vector[F] log_alpha;
  for (f in 1:F)
    log_alpha[f] = log_mu + stand_effect_raw[f] * sigma;

  // log-likelihoods
  matrix[2,N] log_omega;
  matrix[2,2] Gamma;

  // normalize transition matrix
  for (i in 1:2)
    Gamma[i] = Gamma_raw[i] / sum(Gamma_raw[i]);

  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];

    for (t in start_id:end_id){
      log_omega[1,t] = poisson_log_lpmf(y[t] | log_lambda + log(area[t]));
      log_omega[2,t] = neg_binomial_2_log_lpmf(y[t] | log_alpha[f] + log(area[t]), phi);
    }
  }
}

model {
  // Priors
  log_lambda ~ normal(0, log(5)/2.57); // approx between 0 and 5
  log_mu ~ normal(log(200), 0.1); // approx between 150 and 250
  sigma ~ normal(0, 0.5/2.57);
  
  stand_effect_raw ~ normal(0,1);//newly added; standardized deviations
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
      if(state[t] == 1){
        y_rep[t] = poisson_log_rng(log_lambda+ log(area[t]));
      }else{
        y_rep[t] = neg_binomial_2_log_rng(log_alpha[f]+ log(area[t]), phi);
      }
          
    }
  }
}




