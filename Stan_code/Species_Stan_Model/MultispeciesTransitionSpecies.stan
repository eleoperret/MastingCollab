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

transformed data { //new to cause less issues
  real baseline_area = min(area);
}

parameters {
  simplex[2] rho;// initial state probabilities
  
  vector<lower=0, upper=1>[S] theta1; // one value per species
  vector<lower=0, upper=1>[S] theta2;

  array[S] ordered [2] log_means;
  
  //Species-level emission parameters
  //vector[S] log_lambda;// mean log_count for low years
  //vector[S] log_mu;// mean log_count for high years
  
  //Stand level random effects (partial pooling?)
  vector[F] stand_effect_raw;// stand deviations
  real<lower=0> sigma;// stand SD
  
  //Dispersion on a log_scale to ensure positive values
  vector[S] log_phi1;// dispersion for low years
  vector[S] log_phi2;// dispersion for high yars
}

transformed parameters {
  //Back-transform dispersion
  vector<lower=0>[S] phi1 = exp(log_phi1);
  vector<lower= 0>[S] phi2 = exp(log_phi2);
  
  //Stand-level high state mean
  vector[F] log_alpha;
  for (f in 1:F){
    int s= sp[start_idxs[f]]; //species of this sequence
    log_alpha[f] = log_means [s][2] + stand_effect_raw[f] * sigma;
  }
    
//Build stand *species transition matrix
//Gamma [f] is the transitoin matrix for stand f
//Using the species of that stand
  array[S] matrix[2, 2] Gamma;
for (s in 1:S) {
  Gamma[s][1, 1] = theta1[s];
  Gamma[s][1, 2] = 1 - theta1[s];
  Gamma[s][2, 1] = 1 - theta2[s];
  Gamma[s][2, 2] = theta2[s];
}
  
  //Emission log_likelihood
  matrix[2, N] log_omega;
    // compute log-likelihood for each year
  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    for (t in start_id:end_id){
      int s = sp[t];
      real log_area_ratio = log(area[t]/ baseline_area);
      
      log_omega[1,t] = neg_binomial_2_log_lpmf(y[t] | log_means [s][1] + log_area_ratio, phi1 [s]);
      
      log_omega[2,t] = neg_binomial_2_log_lpmf(y[t] | log_alpha[f] + log_area_ratio, phi2 [s]);
    }
  }
}

model {
  // Priors
  
  //Initial state
  rho ~ dirichlet(rep_vector(1.0,2));
  
  theta1 ~ beta(3, 1);  // mean 0.75 — species tend to stay in low state
  theta2 ~ beta(2, 2);  // mean 0.50 — less certain about high state persistence
  
  //Emission priors
  for (s in 1:S){
  log_means [s][1] ~ normal(1.5, 2.0);//exp(2)~7 for low state
  log_means [s][2] ~ normal(5, 1.5) ; 
  }
 
  
  //stand random effects
  sigma ~ normal(0, 0.5);
  stand_effect_raw ~ normal(0,1);
  
  //Dispersion priors on log scale
  log_phi1 ~ normal(log(6.5), 0.5);
  log_phi2 ~ normal(log(6.5), 0.5);

  // likelihood
  for (f in 1:F) {
  int start_id = start_idxs[f];
  int end_id   = end_idxs[f];
  int s        = sp[start_id]; // species of this stand
  target += hmm_marginal(log_omega[, start_id:end_id], Gamma[s], rho); // correct
  }
}

generated quantities {
  array[N] int<lower=0> y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F){
    int start_id = start_idxs[f];
    int end_id = end_idxs[f];
    int s_stand  = sp[start_id]; 
    
    state[start_id:end_id] = hmm_latent_rng(log_omega[, start_id:end_id], Gamma[s_stand], rho
    );
      
    for (t in start_id:end_id){
      int s = sp[t]; //specie of this observation
      real log_area_ratio = log (area[t]/baseline_area);

      if (state[t] == 1){
        //low state
        y_rep[t] = neg_binomial_2_log_rng(log_means [s][1] + log_area_ratio, phi1 [s]);
        
      } else {
        
        //high state
        y_rep[t] = neg_binomial_2_log_rng(log_alpha[f] + log_area_ratio, phi2 [s]);
      }
    }
  }
}

