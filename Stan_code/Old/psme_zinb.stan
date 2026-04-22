//Stan script for my zero-inflated negative binomial model.

data {
  int<lower=1> N;// number of observations
  vector[N] elev;// predictor :scaled elaation
  array[N] int<lower=0> y;// response: total seeds
}
//In this block I define the data that will be fed to my model
//int<lower=1> N. N is the number of observations, here it is the number of stands I have (n=18)
//vector[N] elev. It is a numeric vetor of length N and the predictor (elevation) which was scaled in the R script. and it means that each stand has one elevation value.
//int<lower=0> y[N]. Y is the vector with a lenght of N and each value is an interger above or equal to 0. Which is my response variable total seeds and because seeds are count the value cannot be below 0.

parameters {
  real alpha;// intercept for log(count) part
  real beta;// slope of elevation
  real<lower=0> phi;// NB dispersion
 // real<lower=0, upper=1> theta; // probability of structural zero
 real gamma0; //intercept 
 real gamma1; //slope
}

transformed parameters {
  vector[N] theta_n;
  for (n in 1:N)
  theta_n[n]= inv_logit(gamma0 + gamma1 *elev[n]);//Structural zeros depending on elevation
}
//This block is what STAN estimates
//real is a real number


//This is the model block. THe assumption and the learning from my data.
//It defines my bayesian model 

//liklihood is basically my model 
// I loop over each observation in the dataset. and then I need to consider if the 0 is a structural 0 (an observation that is always 0) or a sampling zero (a zero that come from the negative binomial distribution by chance)
//theta is the mixing properties for structural zeros. and log mix comupter log(p*exp(log_f1) + (1-p)*exp(log_f2)) to see if the zero is a structural or if it comes from the negative binomial


model {
  // Priors
  alpha ~ normal(1.5, 0.5); 
  beta ~ normal(-0.5, 0.5); 
  //phi ~ gamma(2,2);
  phi~lognormal(-1, 0.5);
  gamma0 ~ normal(-2,1);
  gamma1 ~normal(2,1);

  // Likelihood
  for (n in 1:N) {
   if (y[n] == 0) {
     target += log_mix(theta_n[n],0,//log_prob for structural zero
     neg_binomial_2_log_lpmf(y[n] | alpha + beta * elev[n], phi)
     );
   } else {
     y[n] ~ neg_binomial_2_log(alpha + beta * elev[n], phi);
  }
}
  
}






//This is to simulated new data from the model and I will use that to do my prior predictive check (before I feed my real y) and the posterior predictive check (after fitting the model)
generated quantities {
  array[N] int y_sim;//a new simulated observation
  for (n in 1:N) {
    //draw from Bernoulli with observation_specfis theta
    if (bernoulli_rng(theta_n[n]) == 1)
      y_sim[n] = 0;//structural 0
    else
    //draws from negative binomial now
      y_sim[n] = neg_binomial_2_log_rng(alpha + beta * elev[n], phi);
  } // for loop
} // generated quantities
