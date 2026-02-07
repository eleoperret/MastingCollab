// Hurdle Negative Binomial Model for One Species
// Includes stand, year random effects, and previous-year seed production

//DATA BLOCK
data {
  int<lower=1> N;               // number of observations
  array[N] int<lower=0> Seeds;  // seed counts
  int<lower=1> N_Stand;         // number of stands (define ranfom effects dimensions)
  array[N] int<lower=1> Stand;  // stand ID
  int<lower=1> N_Year;          // number of years (define ranfom effects dimensions)
  array[N] int<lower=1> Year;   // year ID
  array[N] real lag_Seeds;       // previous-year seed production
}

//PARAMETER BLOCK
parameters {
  // Fixed effects
  real alpha_occ;                 // baseline logit probability (hurdle); baseline probability of seed being present on the logit scale
  real beta_count;                // baseline log seed count (log scale); baseline expected log seed count if present 
  real gamma_lag;                 // effect of previous-year seeds on occurrence on current year

  // Random effects
  vector[N_Stand] alpha_stand;   // stand effect on occurrence; random intercepts for stands when absence
  vector[N_Stand] beta_stand;    // stand effect on counnt ; random intercept for stand when presence
  vector[N_Year] alpha_year;     // year effect on occurrence; random intercept for when absent
  vector[N_Year] beta_year;      // year effect on count; random intercept for when present. 

  // Standard deviations for random effects for all parameters. There is no random effect for the fixed effect as they are the same for all observations, they to do not vary across stands or years.Whereas the random effect vary by group. And each effect are assumed to come from a normal distribution with mean 0 and standard deviation sigma. 
  real<lower=0> sigma_alpha_stand;
  real<lower=0> sigma_beta_stand;
  real<lower=0> sigma_alpha_year;
  real<lower=0> sigma_beta_year;

  // Negative binomial dispersion; dispersion parameter for negative binomial that controls for overdispersion
  real<lower=0> phi;
}

model {
  // Priors
  alpha_occ ~ normal(0.3, 0.5);//centering probaility closer to 50-60% = baseline occurence from the data. 
  beta_count ~ normal(log(4), 0.5);// mean count,
  gamma_lag ~ normal(0, 0.2);//lag effect

  alpha_stand ~ normal(0, sigma_alpha_stand);
  beta_stand ~ normal(0, sigma_beta_stand);
  alpha_year ~ normal(0, sigma_alpha_year);
  beta_year ~ normal(0, sigma_beta_year);

//smaller variability than 1
  sigma_alpha_stand ~ normal(0, 0.5);
  sigma_beta_stand ~ normal(0, 0.5);
  sigma_alpha_year ~ normal(0, 0.5);
  sigma_beta_year ~ normal(0, 0.5);

  phi ~ normal(5, 1);  // dispersion for NB or try gamma (2,0.5)

  // Likelihood
  for (i in 1:N) {
    if (Seeds[i] == 0) {
      // Probability of zero
      target += bernoulli_lpmf(0 | inv_logit(
        alpha_occ
        + alpha_stand[Stand[i]]
        + alpha_year[Year[i]]
        + gamma_lag * lag_Seeds[i]
      ));
    } else {
      // Probability of >0 (hurdle) * count
      real p_occ = inv_logit(
        alpha_occ
        + alpha_stand[Stand[i]]
        + alpha_year[Year[i]]
        + gamma_lag * lag_Seeds[i]
      );
      target += bernoulli_lpmf(1 | p_occ);  // hurdle part
      target += neg_binomial_2_lpmf(Seeds[i] | exp(beta_count+ beta_stand[Stand[i]]+ beta_year[Year[i]]),phi); // count part
    }
  }
}

//goes through each observation in my dataset(trap/ year/species combo) and calculates the likelihood of what was observed (zero or counts).
//For each observations i in the dataset, is Seeds = 0 then only I care about the probability of no seeds appearing: given the stand, year and previous year seed effect, how likely is it that I observe no seeds in this trap? with alpha_occ, + random effect + lag / invit_logit converts it to a probability between 0 and 1 and bernoulli_lmpf gives the lo probability of observving 0 seeds


generated quantities {
  array[N] int Seeds_prior;
  for (i in 1:N) {
    real p_occ = inv_logit(alpha_occ
                           + alpha_stand[Stand[i]]
                           + alpha_year[Year[i]]
                           + gamma_lag * lag_Seeds[i]);
    if (bernoulli_rng(p_occ) == 0)
      Seeds_prior[i] = 0;
    else
      Seeds_prior[i] = neg_binomial_2_rng(exp(beta_count
                                             + beta_stand[Stand[i]]
                                             + beta_year[Year[i]]),
                                          phi);
  }
}


