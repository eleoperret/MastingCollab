#Getting to learn STAN

##Issues with the code. 
##I can't get the model to not overcompute 0's. There are a lot of them and even if the model is merging and there is technically an output I'm not sure it is working. 
##I'm not sure what would be the next good approach
##Should I change my priors?
#I tried changing my theta distribution, I tried changing alpha, I also change the likelihood function but nothing worked. So I will go back to the first one I did because I'm lost on what to do.
#28.01.2026
#Do I need to simulate differently? Like simulated with completely other data and priors matching my expectations without seing the data? Or what should I do ?
#Clean coomments and explanation for .R and .stan 
#Figure out next steps. 

#Analysis on seed production for one specie over elevation (stand used as proxy)
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork) #To lay out two plots together

getwd()
list.files()


seed_data<-read.csv("SeedData_all.csv")


# ABAM --------------------------------------------------------------------
#selecting one species 
Abam_data<-seed_data%>%
  filter(spp=="ABAM")

total_stand<-Abam_data%>%
  group_by(stand)%>%
  summarise(total_seeds=sum(total_viable_sds),na.rm=TRUE)

ggplot(total_stand, aes(x = stand, y = total_seeds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Stand",
       y = "Total viable seeds",
       title = "Total viable seeds per stand") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


#Site elevation based on the handbook
stand_elevations <- c(
  "TO11" = 600,
  "TO04" = 700,
  "TA01" = 700,
  "TB13" = 850,
  "AV02" = 850,
  "AO03" = 900,
  "AG05" = 950,
  "AV06" = 1060,
  "AX15" = 1090,
  "AB08" = 1100,
  "AV14" = 1150,
  "PP17" = 1150,
  "AM16" = 1200,
  "AE10" = 1450,
  "AR07" = 1450,
  "PARA" = 1600,
  "SPRY" = 1700,
  "SUNR" = 1800
)

#Adding the elevation of the stands to the results
#Creating a dataframe with the stand data
elev_df <- data.frame(
  stand = names(stand_elevations),
  elevation = as.numeric(stand_elevations)
)
#Adding it 
total_stand <- merge(total_stand, elev_df, by = "stand")

ggplot(total_stand,
       aes(x = reorder(stand, elevation), y = total_seeds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Stand (ordered by elevation)",
       y = "Total viable seeds") +
  theme(axis.text.x = element_text(angle = 90))
#Really interesting as there are less year of data for SUNR than for other stands. But I have to admit that there is also the possibility that the ABAM seeds where misidentified as ABLA seeds maybe I should try with another specie? 

# PSME --------------------------------------------------------------------
#selecting one species 
psme_data<-seed_data%>%
  filter(spp=="PSME")

total_stand_psme<-psme_data%>%
  group_by(stand)%>%
  summarise(total_seeds=sum(total_viable_sds),na.rm=TRUE)

ggplot(total_stand_psme, aes(x = stand, y = total_seeds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Stand",
       y = "Total viable seeds",
       title = "Total viable seeds per stand") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#Adding elevation
total_stand_psme <- merge(total_stand_psme, elev_df, by = "stand")

ggplot(total_stand_psme,
       aes(x = reorder(stand, elevation), y = total_seeds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Stand (ordered by elevation)",
       y = "Total viable seeds") +
  theme(axis.text.x = element_text(angle = 90))
#much more what I would expect. I will use this specie for the moment. 



# Modelling production/elevation with Stan --------------------------------
str(psme_data)

#Key points: 
#Everything that has to do with definin the distribution, the prior, the parameter, the likelihood is in the .stan file
#Anything with data prep, scaling, plotting, extracting results is in R. 

###STEP 1 : Understand the data
#Attaching my elevation data to my dataset.
psme_data <- merge(psme_data, elev_df, by = "stand")

#Checking the results with a scatterplot to see my distribution 
ggplot(psme_data, aes(x = elevation, y = total_viable_sds)) +
  geom_jitter(width = 20, height = 0, alpha = 0.3) +
  theme_bw() +
  labs(x = "Elevation (m)",
       y = "Total viable seeds (PSME)")
#Looks like the seed production is mostly around 1000m. 

#Checking my data 
mean(psme_data$total_viable_sds)   # 3.998735
sd(psme_data$total_viable_sds)     # 9.753935
range(psme_data$total_viable_sds)  # 0-87
#The variance is much larger than the mean which indicates that the data is overdispersed. 

#Checking for the amount of 0's in my data
# Count number of zeros
num_zeros <- sum(psme_data$total_viable_sds == 0)
# Total number of observations
total_obs <- nrow(psme_data)
# Fraction of zeros
fraction_zeros <- num_zeros / total_obs
fraction_zeros #0.5243517
#52% of my observations are zero seeds (which is a lot of 0'S)


#I need to scale my elevation for STAN stability? why? because it keeps the predictor around 0 and the SD= 1. Large numbers can cause instability. 
#Is instability bad?
psme_data$elev_sc <- scale(psme_data$elevation)[,1]
#and build a stan list where I define my N, y etc...
stan_data <- list(
  N = nrow(psme_data),
  y = psme_data$total_viable_sds,
  elev = psme_data$elev_sc
)

#Density plot of my data
#This shows me where my data is concentrated. I can see that most of it is close to 0. 
#X-axis : Each point is one observation (year/stand/trap)
#Y-axis: How concentrated is the data around a specific x -value. So here I can see that my data is highly concentrated around 0. 
ggplot(psme_data, aes(x = total_viable_sds)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(
    title = "Density of total viable seeds",
    x = "Total viable seeds",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(0, 50))

#Another good approach acn be to compare density plots and histogram as with count data histograms can be more informative. 
p1 <- ggplot(psme_data, aes(x = total_viable_sds)) +
  geom_histogram(bins = 30, fill = "grey70") +
  ggtitle("Histogram")
p2 <- ggplot(psme_data, aes(x = total_viable_sds)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  ggtitle("Density")
p1 + p2
#So seed count is relatively low in general. 
#Interpretation: Seeds are usually not high but occasionally some high years or high sites which are shown with the end of the curve. So this could be masting or some sites with higher production?

#Another way of visualizing my data to see the distribution of the seed productio nbased on each stand
ggplot(psme_data, aes(x = total_viable_sds)) +
  geom_density() +
  facet_wrap(~ stand)

summary(psme_data$elev_sc)
range(psme_data$elev_sc)
#This can help me define the range for my beta prior




###STEP 2: Choosing the likelihood distribution (the model I will use)
#My options
# a) I thoutgh that I could try first with a *poisson distribution* as it is count data but as the poisson distribution assumes that the mean = var (see line 126-127) I think it is not appropriate. 
# b) I could try now a *negative binomial* (adding a dispersion parameter so that my variance > mean which is useful for overdispersion of the seed counts.) But I have an excess of 0'S (see line 138) and this my underestimate the probabiliy of 0's
# c) I can also try a *zero-inflated Negative Binomial* if I have a lot of 0's which is the case. I think it is the best option as it can handle overdispersion and many 0's 

#Model justification : why log? Seed counts are non-negative integers and doing exp() ensure that alpha + beta * elevation is always =or > than 0. but the multiplicative effect (from beta) becomes additive so it is easier to interpret and alows for a stable model. It stabilizes the variance for overdispesed count data


###In .stan file now

###STEP 3: Specify model structure and choosing priors.  Creating my STAN file that I will name (psme_poisson.stan) where I choose my model structure and my priors (the set of values that are possible) and where I can write my prior predictive check (~are the model assumptions correct~shapiro test)


###STEP 4 : Prior predictive check. 
#If simulated count are way to high or all zeros than my priors need adjusting
#Once my prior predictive check looks good I can fit my model with the real data 

stan_data_ppc <- list(
  N = nrow(psme_data),
  y = rep(0, nrow(psme_data)),  # dummy values instead of the real ones, creates a vecto of 0 to satisfy STAN data requirement. 
  elev = psme_data$elev_sc
)

#This is getting what I wrote in the .rstan file
# Translating my .stan code into C++ code and then compiling the C++ into a standalone executable and then saving it so it can run the MCMC sampling. 
mod <- cmdstan_model("psme_zinb.stan")

fit_ppc <- mod$sample(
  data = stan_data_ppc,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
)

print(fit_ppc, variables = c("alpha", "beta", "phi", "theta"))

#How to interpet this : 
#mad is the median absolute deviation of posterior draws. It measrures spread and sometimes reported as an alternative to SD. Is it important? 
#q5 and q95 are the 5th and 95th percentiles of the posterior draws they form the 90% interval where most of the parameter values lie given the data. 
#rhat is the gelman-rubin convergence diagnostic and compares within chain vs between chain variance (1= chains converge well) 
#Ess-bulk is the effective sample size and measures how many independent draws your posterior is equivalent to for the main bulk of the distribution (not sure I understand). High ESS is reliable estimates, low ess is high autocorrelation. 
#ESS-tail: effective sample for the tail of the posterior, important for credible intervals (especially if the distribution is skewdwd). Low ess tail is uncertainity in the extremes. 



# Extract prior predictive draws
y_sim <- fit_ppc$draws("y_sim")
#Plot prior predictive check
# Extract y_sim as 2D because the ysim is actually 3D I need to reduce: iterations x observations (combine all chains automatically)
yrep <- as_draws_matrix(fit_ppc$draws("y_sim"))
# 50 draws for plotting
yrep <- yrep[1:200, ]
# Plot density overlay
ppc_dens_overlay(y = psme_data$total_viable_sds, yrep = yrep) +
  ggtitle("Prior Predictive Check: PSME seed counts") +
  xlab("Total viable seeds") +
  ylab("Density") +
  coord_cartesian(xlim = c(0, 50))   # optional: zoom on realistic seed range


###STEP 5: Fitting the model to the real data
#Stan will estimate a posterior distribution
#Stan uses both priors an dlikelihood from real data. 
#the y_sim in generated quantiles (.stan file) is now a posterior predictive draw

fit <- mod$sample(
  data = stan_data,   # this includes real y
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)


###STEP 6 : Posterior predictive check. 
# ~ similar to residual diagnostic pots 
y_post <- fit$draws("y_sim")  # posterior predictive

ypost <- as_draws_matrix(fit$draws("y_sim"))#same as before with reducing the matrix size

ypost_plot <- ypost[1:200, ] #reducing the layout 

ppc_dens_overlay(
  y = psme_data$total_viable_sds,  # observed seed counts
  yrep = ypost_plot
) +
  ggtitle("Posterior Predictive Check: Density Overlay") +
  xlab("Total viable seeds") +
  ylab("Density") +
  coord_cartesian(xlim = c(0, 100))   
#It looks like my model is creating way to many 0's

#I want to see how much 0's is my model computing
frac_zero <- function(x) mean(x == 0) #creating a function
zero_fracs <- apply(ypost, 1, frac_zero)
summary(zero_fracs)

hist(zero_fracs, breaks = 20, xlim= c(0,1))
abline(v = mean(psme_data$total_viable_sds == 0), col = "red", lwd = 2)

#Model output
fit$summary()
fit$cmdstan_diagnose()





#Old Stan likelihood function 

# for (n in 1:N) {
#   if (y[n] == 0) {
#     target += log_mix(theta,
#                       bernoulli_lpmf(1 | 1),
#                       neg_binomial_2_log_lpmf(y[n] | alpha + beta * elev[n], phi)
#     );
#   } else {
#     y[n] ~ neg_binomial_2_log(alpha + beta * elev[n], phi);
#   }
# }


#To figure out what distribution I can use for my priors 

#Type of distribution for the theta parameter (probability of 0's)
curve(dbeta(x,1,1), 0,1)   # uniform, mean 0.5, max variance
curve(dbeta(x,2,2), 0,1)   # symmetric, mean 0.5, less spread
curve(dbeta(x,5,5), 0,1)   # symmetric, mean 0.5, tighter
curve(dbeta(x,10,5), 0,1)  # skewed, mean 0.667

#Type of distribution for the phi parameter (dispersion)
curve(dexp(x,1),0,5)        # exponential
curve(dlnorm(x,-1,0.7),0,5) # lognormal
curve(dgamma(x,2,1),0,5)    # gamma


#A prior's check

alpha_draws <- rnorm(1000, mean = 1.4, sd = 1)  # log mean
mu_draws <- exp(alpha_draws)
summary(mu_draws)  # check range

theta_draws <- rbeta(1000, 2, 2)
hist(theta_draws)  # mostly between 0.1 and 0.9

alpha <- rnorm(1000, 1.4, 1)
beta  <- rnorm(1000, 0, 0.5)
phi   <- rgamma(1000, 2, 0.5)
theta <- rbeta(1000, 2, 2)

summary(phi)   # gamma(2,0.5) → mean 4
summary(beta)  # mean 0, SD 0.5


#Simulationg priors
set.seed(123)

n_draws <- 1000      # number of prior draws
n_obs   <- 100       # number of observations to simulate
elev    <- runif(n_obs, -2, 2)  # scaled elevation

# Draw parameters from priors
alpha <- rnorm(n_draws, 1.5, 1)
beta  <- rnorm(n_draws, 0, 0.5)
phi   <- rgamma(n_draws, 2, 0.5)   # NB overdispersion
theta <- rbeta(n_draws, 2, 2)      # probability of structural zero

#Simulating from the negative binomial ditstribution
y_sim <- matrix(NA, nrow = n_draws, ncol = n_obs)

for(i in 1:n_draws){
  mu <- exp(alpha[i] + beta[i] * elev)  # NB mean for each obs
  for(j in 1:n_obs){
    if(runif(1) < theta[i]){
      y_sim[i,j] <- 0                 # structural zero
    } else {
      # NB: rnbinom parametrized with size=phi, mean=mu
      y_sim[i,j] <- rnbinom(1, size = phi[i], mu = mu[j])
    }
  }
}

#Checking for the 0's
frac_zero <- apply(y_sim, 1, function(x) mean(x==0))
summary(frac_zero)
hist(frac_zero, breaks = 20, col = "grey",
     main = "Fraction of zeros from prior predictive draws",
     xlab = "Fraction of zeros")

# Add observed fraction of zeros for reference
observed_frac_zero <- 0.524  # from your data
abline(v = observed_frac_zero, col = "red", lwd = 2)

#Checking distribution of simulated seeds
y_all <- as.vector(y_sim)

hist(y_all, breaks = 50, col = "lightblue",
     main = "Prior predictive distribution of seeds",
     xlab = "Simulated seed counts")

