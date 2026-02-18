// 2-state HMM for seed production (single stand, T years)
data {
  int<lower=1> T;            // number of years
  array [T] int<lower=0> y;         // seed production per year
}

parameters {
  simplex[2] rho;            // initial state probabilities
  matrix<lower=0>[2,2] Gamma_raw; // raw transition matrix (rows to be normalized)
  ordered[2] log_mu;         // log mean seed production for low/high states
  vector<lower=1e-6>[2] phi; // dispersion per state
}

transformed parameters {
  vector[2] mu = exp(log_mu);       // natural-scale means
  matrix[2,T] log_omega;            // log-likelihood for each year × state
  matrix[2,2] Gamma;                // normalized transition matrix

  // Normalize rows of Gamma_raw to sum to 1
  for (i in 1:2)
    Gamma[i] = Gamma_raw[i] / sum(Gamma_raw[i]);

  // Compute log-likelihood for each state/year
  for (t in 1:T)
    for (s in 1:2)
      log_omega[s,t] = neg_binomial_2_lpmf(y[t] | mu[s], phi[s]);
}

model {
  // Priors
  log_mu ~ normal(4, 1);        // allows seed production roughly up to ~400
  phi ~ gamma(2, 0.1);          // moderate overdispersion
  to_vector(Gamma_raw) ~ exponential(1); // weakly informative positive prior

  // HMM likelihood
  target += hmm_marginal(log_omega, Gamma, rho);
}

generated quantities {
  array[T]int y_rep;       // posterior predictive seeds
  array[T]int z;           // latent state sequence

  z = hmm_latent_rng(log_omega, Gamma, rho);  // simulate latent states

  for (t in 1:T)
    y_rep[t] = neg_binomial_2_rng(mu[z[t]], phi[z[t]]);
}

