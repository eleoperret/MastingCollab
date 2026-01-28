#Getting to learn STAN

##Issues with the code. 
##I can't get the model to not overcompute 0's. There are a lot of them and even if the model is merging and there is technically an output I'm not sure it is working. 

#Analysis on seed production for one specie over elevation (stand used as proxy)
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(posterior)

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

###STEP 2: Choosing the likelihood distribution (the model I will use)
#My options
# a) I thoutgh that I could try first with a *poisson distribution* as it is count data but as the poisson distribution assumes that the mean = var (see line 126-127) I think it is not appropriate. 
# b) I could try now a *negative binomial* (adding a dispersion parameter so that my variance > mean which is useful for overdispersion of the seed counts.) But I have an excess of 0'S (see line 138) and this my underestimate the probabiliy of 0's
# c) I can also try a *zero-inflated Negative Binomial* if I have a lot of 0's which is the case. I think it is the best option as it can handle overdispersion and many 0's 


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
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9766  0.9962  0.9981  0.9975  0.9994  1.0000
hist(zero_fracs, breaks = 20, xlim= c(0,1))
abline(v = mean(psme_data$total_viable_sds == 0), col = "red", lwd = 2)

#Model output
print(fit, variables = c("alpha", "beta"))

fit$summary()
fit$cmdstan_diagnose()




















library(bayesplot)
library(ggplot2)

#--------------------------
# Simulated example
#--------------------------

# Observed data: y (10 observations)
y <- c(2, 3, 0, 5, 1, 4, 2, 0, 3, 6)

# Posterior predictive draws: yrep (5 draws, same 10 observations)
set.seed(123)
yrep <- matrix(
  rpois(50, lambda = 3),   # 50 simulated values = 5 draws * 10 observations
  nrow = 5,
  ncol = 10
)

#--------------------------
# Posterior predictive density overlay
#--------------------------
ppc_dens_overlay(y = y, yrep = yrep) +
  ggtitle("Illustration: Observed y (black) vs Multiple yrep draws (colored)") +
  xlab("Counts") +
  ylab("Density")





# 
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