data {
  int<lower=1> N;// total observations
  int<lower=1> T;// number of years
  array[N] int<lower=0> y;//seed prodction

  int<lower=1> n_stand;//because there are changes based on the stand
  int<lower=1> n_trap;//because there are changes based on the trap

  array[N]int<lower=1, upper=T> year_id;
  array[N]int<lower=1, upper=n_stand> stand_id;
  array[N] int<lower=1, upper=n_trap> trap_id;
}

parameters {

  // HMM, creating a matrix with both states, low and high. 
  simplex[2] init_state;
  array[2] simplex[2] trans;
  
  ordered[2] alpha;//because one state is higher than the other

  vector[n_stand] beta_stand;

  vector[n_trap] u_trap;
  real<lower=0> sigma_trap;

  // Overdispersion
  real<lower=0> phi;
}


model {

  // Priors
  beta_stand ~ normal(0, 1);

  u_trap ~ normal(0, sigma_trap);
  sigma_trap ~ normal(0, 1);

  phi ~ normal(1, 1);

  // Transition matrix prior 
  for (s in 1:2)
    trans[s] ~ dirichlet(rep_vector(1.0, 2));

  init_state ~ dirichlet(rep_vector(1.0, 2));

  matrix[T, 2] log_emission;

  for (t in 1:T) {
    for (s in 1:2) {

      real sum_ll = 0;

      for (n in 1:N) {
        if (year_id[n] == t) {

          real eta =
            alpha[s]
            + beta_stand[stand_id[n]]
            + u_trap[trap_id[n]];

          real mu = exp(eta);

          sum_ll += neg_binomial_2_lpmf(y[n] | mu, phi);
        }
      }

      log_emission[t, s] = sum_ll;
    }
  }


  vector[2] log_alpha;
  vector[2] log_alpha_new;

  // Initial year
  for (s in 1:2)
    log_alpha[s] =
      log(init_state[s]) +
      log_emission[1, s];

  // Recursion
  for (t in 2:T) {
    for (s in 1:2) {

      vector[2] temp;

      for (sp in 1:2)
        temp[sp] =
          log_alpha[sp] +
          log(trans[sp][s]);

      log_alpha_new[s] =
        log_emission[t, s] +
        log_sum_exp(temp);
    }

    log_alpha = log_alpha_new;
  }

  target += log_sum_exp(log_alpha);
}





generated quantities {

  matrix[T, 2] log_emission;

  // Recompute emission log-likelihoods
  for (t in 1:T) {
    for (s in 1:2) {

      real sum_ll = 0;

      for (n in 1:N) {
        if (year_id[n] == t) {

          real eta =
            alpha[s]
            + beta_stand[stand_id[n]]
            + u_trap[trap_id[n]];

          real mu = exp(eta);

          sum_ll += neg_binomial_2_lpmf(y[n] | mu, phi);
        }
      }

      log_emission[t, s] = sum_ll;
    }
  }

  // Forward pass
  matrix[T, 2] log_alpha;

  for (s in 1:2)
    log_alpha[1, s] =
      log(init_state[s]) +
      log_emission[1, s];

  for (t in 2:T) {
    for (s in 1:2) {

      vector[2] temp;

      for (sp in 1:2)
        temp[sp] =
          log_alpha[t-1, sp] +
          log(trans[sp][s]);

      log_alpha[t, s] =
        log_emission[t, s] +
        log_sum_exp(temp);
    }
  }

  // Backward pass
  matrix[T, 2] log_beta;

  for (s in 1:2)
    log_beta[T, s] = 0;

  for (t in (T-1):1) {
    for (s in 1:2) {

      vector[2] temp;

      for (sp in 1:2)
        temp[sp] =
          log(trans[s][sp]) +
          log_emission[t+1, sp] +
          log_beta[t+1, sp];

      log_beta[t, s] =
        log_sum_exp(temp);
    }
  }

  // Posterior state probabilities
  matrix[T, 2] state_prob;

  for (t in 1:T) {

    vector[2] temp;

    for (s in 1:2)
        temp[s] = log_alpha[t, s] + log_beta[t, s];

    // Stable softmax
    real max_temp = max(temp);
    temp = temp - max_temp;
    vector[2] probs = exp(temp) / sum(exp(temp));

    state_prob[t, 1] = probs[1];
    state_prob[t, 2] = probs[2];
  }
}

