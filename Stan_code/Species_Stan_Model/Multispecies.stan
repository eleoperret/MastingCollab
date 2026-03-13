// 2-state HMM for seed production with stand-specific transitions
// Trap area scaling and stand-level random effects
// Two NB as distribution for the low and high seed production


data {
  int<lower=1> N;// number of observations
  int<lower=1> F;// number of stands
  int<lower=1> S; //number of species
  
  array[N] int<lower=0> y;// seed production per year
  array[N] int<lower=1, upper= S> sp; //species for each observations tells stan which species produced ibservation t
  
  array[F] int<lower=1> start_idxs; // start index per stand
  array[F] int<upper=N> end_idxs;// end index per stand
  
  vector<lower=0>[N] area;// total trap area per stand-year
}

parameters {
  simplex[2] rho;// initial state probabilities
  //NEW
  real<lower=0, upper=1> theta1[F]; // State 1
  real<lower=0, upper=1> theta2[F]; // State 2
  
  vector[S] log_lambda;// mean for low years
  vector[S] log_mu;// base for high years
  
  vector[F] stand_effect_raw;// stand deviations
  real<lower=0> sigma;// stand SD
  
  vector<lower=0>[S] phi1;// dispersion for low years
  vector<lower=0>[S] phi2;// dispersion for high yars
}

transformed parameters {
  real baseline_area = min(area);

  vector[F] log_alpha;
  for (f in 1:F){
    int start_id = start_idxs[f];
    int s= sp[start_id]; //species of this sequence
    log_alpha[f] = log_mu [s] + stand_effect_raw[f] * sigma;
  }
    

  matrix[2,2] Gamma[F]; // stand-specific transition matrices
  matrix[2,N] log_omega; // log-likelihood per year × state

  for (f in 1:F){
    Gamma[f][1,1] = theta1[f];
    Gamma[f][1,2] = 1 - theta1[f];
    Gamma[f][2,1] = 1 - theta2[f];
    Gamma[f][2,2] = theta2[f];
  }

  // compute log-likelihood for each year
  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    for (t in start_id:end_id){
      int s = sp[t];
      
      log_omega[1,t] = neg_binomial_2_log_lpmf(y[t] | log_lambda [s] + log(area[t]/ baseline_area), phi1 [s]);
      
      log_omega[2,t] = neg_binomial_2_log_lpmf(y[t] | log_alpha[f] + log(area[t]/ baseline_area), phi2 [s]);
    }
  }
}

model {
  // Priors
  rho ~ dirichlet(rep_vector(1.0,2));// initial state, equal probability of being mast or non-mast
  for (f in 1:F){
    theta1[f] ~ beta(3,1);// mean ~ 0.83 before 5,1 
    theta2[f] ~ beta(2,2); // mean ~ 0.17 before 1,5
  }

  log_lambda ~ normal(0, log(30)/2.57);
  log_mu ~ normal(log(500), 1) ; // wider as specie 5 has in the thousands seeds (thse)
  
  sigma ~ normal(0, 0.5/2.57);
  stand_effect_raw ~ normal(0,1);
  
  phi1 ~ normal(6.5, 4);
  phi2 ~ normal(6.5, 3);

  // likelihood
  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    target += hmm_marginal(log_omega[,start_id:end_id], Gamma[f], rho);
  }
}

generated quantities {
  array[N] int<lower=0> y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    
    //specie of this sequence
    int s_seq= sp[start_id];
    real log_alpha_f = log_mu[s_seq] + stand_effect_raw[f]* sigma;
    
    //draw the states for this sequence
    state[start_id:end_id] = hmm_latent_rng(log_omega[,start_id:end_id], Gamma[f], rho);
    
    for (t in start_id:end_id){
      int s = sp[t]; //specie of this observation
      
      if (state[t] == 1){
        //low state
        y_rep[t] = neg_binomial_2_log_rng(log_lambda [sp[t]] + log(area[t]/ baseline_area), phi1 [sp[t]]);
        
      } else {
        
        //high state
        y_rep[t] = neg_binomial_2_log_rng(log_alpha[f] + log(area[t]/ baseline_area), phi2 [sp[t]]);
      }
    }
  }
}

