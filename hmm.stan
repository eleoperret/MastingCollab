// 2-state Hidden Markov Model for seed production

data {
  int<lower=1> S;          // number of stands
  int<lower=1> T;          // number of years
   array[S, T] int<lower=0> y;    // aggregated seeds per stand x year
}

parameters {
  ordered[2] mu_state;             // log mean seed production per state (low < high)
  array[2] simplex[2] P;           // transition matrix: rows sum to 1
  simplex[2] pi;                    // initial state probabilities
  real<lower=0> phi;                // overdispersion
}

model {
  // Priors
  mu_state ~ normal(0, 2);
  phi ~ exponential(1);
  for (i in 1:2) {
    P[i] ~ dirichlet(rep_vector(1.0, 2));
  }
  pi ~ dirichlet(rep_vector(1.0, 2));

  // Likelihood via forward algorithm
  for (s in 1:S) {
    vector[2] log_alpha;       // log forward probabilities
    // first year
    for (k in 1:2) {
      log_alpha[k] = log(pi[k]) + neg_binomial_2_log_lpmf(y[s,1] | mu_state[k], phi);
    }
    // subsequent years
    for (t in 2:T) {
      vector[2] log_alpha_new;
      for (k in 1:2) {
        vector[2] temp;
        for (j in 1:2) {
          temp[j] = log_alpha[j] + log(P[j][k]);  // note: array-of-simplex syntax
        }
        log_alpha_new[k] = log_sum_exp(temp) + neg_binomial_2_log_lpmf(y[s,t] | mu_state[k], phi);
      }
      log_alpha = log_alpha_new;
    }
    target += log_sum_exp(log_alpha);  // marginal likelihood for stand s
  }
}
