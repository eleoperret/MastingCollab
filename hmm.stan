data {
  int<lower=1> N;                 // total observations
  int<lower=1> S;                 // number of stands
  int<lower=1> J;                 // number of traps
  int<lower=1> T;                 // number of years

  array[N] int<lower=1, upper=S> stand;
  array[N]int<lower=1, upper=J> trap;
  array[N] int<lower=1, upper=T> year;

  array[N] int<lower=0> y;              // seeds
  vector[N] log_offset;           // log(trap area or exposure). Is it necessary?
}

parameters {
  // Fixed effects
  real alpha;

  // Random effects
  vector[S] b_stand;
  vector[J] c_trap;
  real<lower=0> sigma_stand;
  real<lower=0> sigma_trap;

  // Latent AR(1) process or my latent continous variable
  matrix[S, T] eta;
  real<lower=-1, upper=1> rho;
  real<lower=0> sigma_eta;

  // Overdispersion parameter because of the negative binomial
  real<lower=0> phi;
}

model {
  // ---- Priors ----
  alpha ~ normal(5, 3);// for psme alpha ~ normal(0, 2)

  sigma_stand ~ normal(0, 1);
  sigma_trap  ~ normal(0, 1);
  sigma_eta   ~ normal(0, 1);

  b_stand ~ normal(0, sigma_stand);
  c_trap  ~ normal(0, sigma_trap);

  rho ~ normal(-0.5, 0.3);   // maybe could be changed? 
  phi ~ exponential(1);

  // ---- Latent process or hidden state ----
  for (s in 1:S) {
    eta[s, 1] ~ normal(0, sigma_eta);
    for (t in 2:T) {
      eta[s, t] ~ normal(rho * eta[s, t - 1], sigma_eta);
    }
  }

  // ---- Emission model ----
  for (n in 1:N) {
    real log_mu;
    log_mu =
      alpha +
      eta[stand[n], year[n]] +
      b_stand[stand[n]] +
      c_trap[trap[n]] +
      log_offset[n];

    y[n] ~ neg_binomial_2_log(log_mu, phi);
  }
}


generated quantities {
  array[N] int y_prior;// for prior predictive check
  array[N] int y_rep;// for posterior predictive check

  for (n in 1:N) {
    real log_mu = alpha + eta[stand[n], year[n]] + b_stand[stand[n]] + c_trap[trap[n]] + log_offset[n];

    y_prior[n] = neg_binomial_2_log_rng(log_mu, phi);  // prior predictive
    y_rep[n]   = neg_binomial_2_log_rng(log_mu, phi);  // posterior predictive
  }
}


