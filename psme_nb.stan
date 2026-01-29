data {
  int<lower=1> N;                 
  vector[N] elev;                 
  array[N] int<lower=0> y;        
}

parameters {
  real alpha;                     
  real beta;                      
  real<lower=0> phi;              
}

model {
  // Priors
  alpha ~ normal(1.5, 0.5);       
  beta  ~ normal(-0.5, 0.5);      
  phi   ~ gamma(2, 2);            

  // Likelihood
  y ~ neg_binomial_2_log(alpha + beta * elev, phi);
}

generated quantities {
  array[N] int y_sim;
  for (n in 1:N) {
    y_sim[n] = neg_binomial_2_log_rng(
      alpha + beta * elev[n],
      phi
    );
  }
}
