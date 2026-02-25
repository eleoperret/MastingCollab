// 2-state HMM for seed production with stand-specific transitions
// Trap area scaling and stand-level random effects

data {
  int<lower=1> N;// number of observations
  int<lower=1> F;// number of stands
  array[N] int<lower=0> y;// seed production per year
  array[F] int<lower=1> start_idxs; // start index per stand
  array[F] int<upper=N> end_idxs;// end index per stand
  vector<lower=0>[N] area;// total trap area per stand-year
}

parameters {
  simplex[2] rho;// initial state probabilities
  //NEW
  real<lower=0, upper=1> theta1[F]; // P(mast->mast) per stand
  real<lower=0, upper=1> theta2[F]; // P(non-mast->non-mast) per stand
  
  real log_lambda;// mean for low years
  real log_mu_raw;// base for high years
  
  vector[F] stand_effect_raw;// stand deviations
  real<lower=0> sigma;// stand SD
  real<lower=0> phi;// dispersion for high years
}

transformed parameters {
  real log_mu;
  log_mu = log_lambda + exp(log_mu_raw);

  vector[F] log_alpha;
  for (f in 1:F)
    log_alpha[f] = log_mu + stand_effect_raw[f] * sigma;

  matrix[2,2] Gamma[F]; // stand-specific transition matrices
  matrix[2,N] log_omega; // log-likelihood per year × state

  for (f in 1:F){
    Gamma[f][1,1] = theta1[f];// mast -> mast
    Gamma[f][1,2] = 1 - theta1[f];// mast -> non-mast
    Gamma[f][2,1] = 1 - theta2[f];// non-mast -> mast
    Gamma[f][2,2] = theta2[f];// non-mast -> non-mast
  }

  // compute log-likelihood for each year
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
  rho ~ dirichlet(rep_vector(1.0,2));// initial state, equal probability of being mast or non-mast
  for (f in 1:F){
    theta1[f] ~ beta(1,5);// low probability mast->mast
    theta2[f] ~ beta(5,1); // high probability non-mast->non-mast
  }

  log_lambda ~ normal(0, log(5)/2.57);
  log_mu_raw ~ normal(0, 1) ; // old normal(log(200), 0.1) changed because this was causing issues when plotting the PPC. Too many NA's produced. 
  sigma ~ normal(0, 0.5/2.57);
  stand_effect_raw ~ normal(0,1);
  phi ~ gamma(2, 0.1);

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
    state[start_id:end_id] = hmm_latent_rng(log_omega[,start_id:end_id], Gamma[f], rho);
    
    for (t in start_id:end_id){
      if (state[t] == 1){
        y_rep[t] = poisson_log_rng(log_lambda + log(area[t]));
      } else {
        y_rep[t] = neg_binomial_2_log_rng(log_alpha[f] + log(area[t]), phi);
      }
    }
  }
}

