//Hierarchical Bayesian state-space model with a continuous latent AR(1) process

//This model estimates the true, unobserved seed production of each stand in each year and where production depend on the previous year and is link to noisy trap counts. 

//y[n]= number of seeds found in trap j, in stand s, and in year t. 

//my seed production is noisy and overdisped, affected by trap size and stand specific differences. so y[n]!= seed measured by us. 

//the true seed production level is eta[s,t] where I will have one number per stand and if the number is high, good production; low = poor production

//Seed prduction is not random from year to year (there is depletion after very good year of resources) which is coded in my eta (eta[s,t-1])
//rho = how strong and in which direction that dependence is : is it masting specie, not etc...
//Sigma_eta is how unpredictibale are my year to year changes. 

//I consider here trap as a unit of sampling not a biological unit but could be discussed as trap layout changed and could be adding some complexity. 

//Structure of the code 1) data, 2) parameters, 3) model and 4) generated quantiles for prior and posterior predictive check 

//Data: what I observe (all my data and how they look like)
data {
  int<lower=1> N;// number of rows in my dataset
  int<lower=1> S;// number of stands
  int<lower=1> J;// number of traps
  int<lower=1> T;// number of years

  array[N] int<lower=1, upper=S> stand;
  array[N]int<lower=1, upper=J> trap;
  array[N] int<lower=1, upper=T> year;

  array[N] int<lower=0> y;// seeds, my actual response variable and it can be 0 or higher. 
  vector[N] log_offset; // log(trap area or exposure). Is it necessary? maybe maybe not? 
}

//Unknowns Stan is going to estimate
parameters {
  // Fixed effects
  real alpha; //overall average log seed production

  // Random effects
  vector[S] b_stand;//stand difference
  vector[J] c_trap;//some traps can catch more seeds
  real<lower=0> sigma_stand;// how different are stands from eachother
  real<lower=0> sigma_trap;//same for trap

  // Latent AR(1) process or my latent continous variable and creates one hidden value for each stand x year
  matrix[S, T] eta;
  real<lower=-1, upper=1> rho;// limited to -1,1 so that the process stays stable
  real<lower=0> sigma_eta;//how noisy are years to year changes

  // Overdispersion parameter 
  real<lower=0> phi;
}

//Model: assumption on how my data behave: here a negative binomial with overdispersion. 

model {
  // ---- Priors ----
  alpha ~ normal(5, 3);// for psme alpha ~ normal(0, 2) but how can I do for all species?
  sigma_stand ~ normal(0, 1);// core bayesian concept where Sd is unknown and I need to get a probability distribution-- partial pooling?
  sigma_trap  ~ normal(0, 1);
  sigma_eta   ~ normal(0, 1);

  b_stand ~ normal(0, sigma_stand);
  c_trap  ~ normal(0, sigma_trap);

  rho ~ normal(-0.5, 0.3);// maybe could be changed? for th emoment I have a negative autocorrealtion (masting ) with uncertainity. For example I could use a more standard one : normal (0, 0.5) or uniform (-1,1): which means I have no idea. 
  phi ~ exponential(1);//some overdispersion

  // ---- Latent process or hidden state ----
  for (s in 1:S) {
    eta[s, 1] ~ normal(0, sigma_eta);
    for (t in 2:T) {
      eta[s, t] ~ normal(rho * eta[s, t - 1], sigma_eta);
    }
  }

  // ---- Emission model ----
  for (n in 1:N) {
    real log_mu;//temporary varaible for expected log mean
    log_mu =
      alpha +
      eta[stand[n], year[n]] +
      b_stand[stand[n]] +
      c_trap[trap[n]] +
      log_offset[n];

    y[n] ~ neg_binomial_2_log(log_mu, phi);
  }
}

//simulate data for checking prior and model
generated quantities {
  array[N] int y_prior;// for prior predictive check
  array[N] int y_rep;// for posterior predictive check

  for (n in 1:N) {
    real log_mu = alpha + eta[stand[n], year[n]] + b_stand[stand[n]] + c_trap[trap[n]] + log_offset[n];

    y_prior[n] = neg_binomial_2_log_rng(log_mu, phi);  // prior predictive
    y_rep[n]   = neg_binomial_2_log_rng(log_mu, phi);  // posterior predictive
  }
}


