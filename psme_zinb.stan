//this is a comment on STAN
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
  real<lower=0, upper=1> theta; // probability of structural zero
}

//This block is what STAN estimates
//real is a real number
//This could also be an interger, a numeric vector or a matrix and I can also constrain it.
//Alpha. There is a unknown quantity or quantities that I want to estimate.It is the expected log seed count when elevation is 0. The baseline seed production 
//Why log? Because the values cannot be negative with a poisson. 
//Beta is the slope. How stronlgy seed production changes with elevation if beta is >0 then the seed production increase otherwise descrease. It is the strength of my effect
//What is my effect? It is the strength at with my response changes when a predictor changes by one unit. 

model {
  // Priors
  alpha ~ normal(1.5, 1); 
  beta ~ normal(0, 0.5); 
  phi ~ gamma(2,0.5);
  theta ~ beta(2, 2);

  // Likelihood
  for (n in 1:N) {
   if (y[n] == 0) {
     target += log_mix(theta,
                       bernoulli_lpmf(1 | 1),
                       neg_binomial_2_log_lpmf(y[n] | alpha + beta * elev[n], phi)
     );
   } else {
     y[n] ~ neg_binomial_2_log(alpha + beta * elev[n], phi);
  }
}
  
}

//last changes 28.01.26 : tried the gamma but looks worse than the exponential and changed the alpha to 1.5 so maybe try again next time. 

//What to try : Alpha 1.4,1/Beta 0,0.5 or (0,1)/theta beta(5,5)or (10,10) or (20,1) or (2,2)/ phi exponential (1) or lognnormal(-1,0.7) or gamma(2,1)

//Justification on how I chose my priors.
//model : log(expected number of seeds)= alpha + beta * elevation
//alpha. As my mean seed is around 4 so if I do log(4) which is 1,4. Then if you look at the density plot my values are mostly around 0-17 seeds with some case of extreme seed production (which is handled by my phi-- overdispersion), so then I should select for a SD that would be within this range of most datapoint which one would be good so I would have alpha ranging between 0,4 and 2,4 which would meanbetween 1,5 and 11 seeds (????maybe I could use another one like 1.1 or 1.3 ????). The SD here has to reflect the plausible variation on the log scale not the extreme outliers that will be handled by phi. 
//theta. Because I know that I have 50% of my data that is 0, this means that I have a and b that are symmetric so using (5,5) means that I have 0,5 for my theta prior which reflects my data. 


//This is the model block. THe assumption and the learning from my data.
//It defines my bayesian model 
//The priors. It is that bfore seing my data, I assume that those are the distribution of my parameters. It is what I believe about a parameter befor seing the data.
//alpha= is on a log scale. So log(4)~1,4 from the mean of the data. It is the set baseline when elev = 0. First the mean of the normal distribution and then the standard deviation(spread of plausible alpha values)
//beta is the slope of elevation and we assume always at the beginning that there is no effect of elevation.(Mean = 0 so no effect, and SD, plausible range)
//phi is the dispersion parameter that allows my variance to be larger than my mean. in a NB, Var(y)= mean (y) + mean(y)^2/phi. MEean = overdispersion around 2 and SD allows phi to vary. and then thrucated at 0 because phi must be positiv
//theta is the probability of a structural zero. It accounts for the excta 0's. Theta = 0 no structural 0, 1 all 0's are structural. Beta is bounded between 0 an 1. Shape of the parameters is 2,2 so roughly symmetric. THis means that structural zeros are equally likely as sampling zeros. 

//liklihood is basically my model 
// I loop over each observation in the dataset. and then I need to consider if the 0 is a structural 0 (an observation that is always 0) or a sampling zero (a zero that come from the negative binomial distribution by chance)
//theta is the mixing properties for structural zeros. and log mix comupter log(p*exp(log_f1) + (1-p)*exp(log_f2)) to see if the zero is a structural or if it comes from the negative binomial

generated quantities {
  array[N] int y_sim;
  for (n in 1:N) {
    if (bernoulli_rng(theta) == 1)
      y_sim[n] = 0;
    else
      y_sim[n] = neg_binomial_2_log_rng(alpha + beta * elev[n], phi);
  } // for loop
} // generated quantities
