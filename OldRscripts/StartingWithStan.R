#Getting to learn STAN
#Analysis on seed production for one specie over elevation (stand used as proxy)
#Eléonore Perret - February 2026

##Issues with the code. 
#What to work on and questions for myself: 
#Worked on the HMM model and it looked really good for PSME but not for TSHE. Also do I change my priors for each species? How would that look like in a multi-species model?
#For HMM I did not use a transition matrix but a latent continuous variable as my hidden state. Is that correct? I have 
#Is my approach correct to check priors? 
#What is a bin probability? and why plotting bin probability vs bin width? Should I use this? 
#Think about my parameters. Are they correlated, are they constrained? 

#Re do with same distribution and then check with a partial pooling for the difference between stands. Check the parameters differnece between model. Create git issue and with the summary adn the description of the model (likelihood) adn thoughts also the posterior check visuals (using mickeals diagnostics, or density together with histogramms). What I see and what to do next. WHcih parameters to change based on the fact that wach stand has its own distribution.

library(dplyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(patchwork) #To lay out two plots together
library(tidyr)

getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")
list.files()


seed_data<-read.csv("SeedData_all.csv")


# SiteInfos ---------------------------------------------------------------


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
  filter(spp=="PSME")%>%
  #Added this constraint
  filter(!stand %in% c("AE10", "AR07", "PARA", "SUNR", "SPRY", "AV14", "AM16"))

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

summary(psme_data$total_viable_sds)
mean(psme_data$total_viable_sds) #6.48657
var(psme_data$total_viable_sds) #139.4104
#This is overdispersed Var>> mean

#I can see here that some stand have 0 seed production in 15 years (very unlikely that it will change) except for SPRY and AR07 (but to be honest I'm a bit suprised; in total it is 2 seeds found for the AR07 and 3 for the SPRY over the last years; what to do with that. I honestly cannot refute the misidentification because there are no adult tree anywhere close to this area so something to think about)

#Changing my dataset to PSME where only stands with seeds -- removing AE10, AR07, PARA, SUNR, SPRY, AV14, AM16


# TSHE --------------------------------------------------------------------
#selecting one species 
tshe_data<-seed_data%>%
  filter(spp=="TSHE")%>%
  #Added this constraint
  filter(!stand %in% c("AE10", "AR07", "PARA", "SUNR", "SPRY"))

total_stand_tshe<-tshe_data%>%
  group_by(stand)%>%
  summarise(total_seeds=sum(total_viable_sds),na.rm=TRUE)

# Modelling production/elevation with Stan  --------------------------------
str(psme_data)

#Key points for myself: 
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
#52% of my observations are zero seeds (which is a lot of 0'S), probably coming a lot from those high elevations that have 0 seed prduction

##Thinking about the elevation now
#I need to scale my elevation for STAN stability? why? because it keeps the predictor around 0 and the SD= 1. Large numbers can cause instability. 
#Is instability bad?
psme_data$elev_sc <- scale(psme_data$elevation)[,1]

length(unique(psme_data$stand))
length(unique(psme_data$elevation))
#Some stands share the same elevation
#Not a problem for this example but I think that it could be a problem for the bigger model as there are some stands with similar elevation yes but with diffirent species composition and slope and orientation (N/S) so I think if we want to use the stand as a proxy for elevation or climatic/environemental conditions it would have to be taken into account

#Thinking more about my data
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
#Interpretation: Seeds are usually not high but occasionally some high years or because of the fact that some sites never produce seeds which are shown with the end of the curve. So this could be masting or an artefact of the discrpency between sites

#Another way of visualizing my data to see the distribution of the seed productio nbased on each stand *Similar as above but here I can see the density distribution for each stand
ggplot(psme_data, aes(x = total_viable_sds)) +
  geom_density() +
  facet_wrap(~ stand)

#And year
ggplot(psme_data, aes(x = total_viable_sds)) +
  geom_density() +
  facet_wrap(~ year)
#Interesting as we can see that some years have very low seed production and others quite high, which correspons to the literature. I think that it is really from what I read where PSME have ~7 years cycle so there is natural variation maybe pushed by other species but the very high years are around this 7 year cycle.
ggplot(psme_data,
       aes(x = year, y = total_viable_sds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Year",
       y = "Total viable seeds") +
  theme(axis.text.x = element_text(angle = 90))


##*Before going to the next step there a few things to consider. Some years have low and high seed production I think this should not be an issue but should I but this as a prior?
##*There are stands where there is no seed production so this is quite tricky because it makes it harder to model or to have a proper model; at least with what I tried




###STEP 2: Choosing the likelihood distribution (the model I will use)
#My options
# a) I thoutgh that I could try first with a *poisson distribution* as it is count data but as the poisson distribution assumes that the mean = var (see line 137) I think it is not appropriate because Var>>Mean. 
# b) I could try a *negative binomial* (adding a dispersion parameter so that my variance > mean which is useful for overdispersion of the seed counts.) But I have an excess of 0'S and this my underestimate the probabiliy of 0's
# c) I think the best option is to try a *zero-inflated Negative Binomial* if I have a lot of 0's which is the case. I think it is the best option as it can handle overdispersion and many 0's 

#Model justification : why log? Seed counts are non-negative integers and doing exp() ensure that alpha + beta * elevation is always =or > than 0. but the multiplicative effect (from beta) becomes additive so it is easier to interpret and alows for a stable model. It stabilizes the variance for overdispesed count data


###STEP 3: Specify model structure and choosing priors.  Creating my STAN file that I will name (psme_poisson.stan) where I choose my model structure and my priors (the set of values that are possible) and where I can write my prior predictive check (~are the model assumptions correct~shapiro test)

# Justification on how I chose my priors.
# model : log(expected number of seeds)= alpha + beta * elevation

#Alpha. (the intercept). 
#It is the expected log seed count when elevation is 0. The baseline seed production 
# alpha~normal(a,b). First the mean of the normal distribution and then the standard deviation(spread of plausible alpha values)
#As my mean seed is around 4 and alpha is on a log scale. So log(4)~1,4. Then if you look at the density plot my values are mostly around 0-17 seeds with some case of extreme seed production (which is handled by my phi-- overdispersion), so then I should select for a SD that would be within this range of most datapoint which one would be good so I would have alpha ranging between 0,4 and 2,4 which would meanbetween 1,5 and 11 seeds (????maybe I could use another one like 1.1 or 1.3 ????). The SD here has to reflect the plausible variation on the log scale not the extreme outliers that will be handled by phi. 

#Beta. (the slope)
#beta is the slope of elevation and we assume always at the beginning that there is no effect of elevation.(Mean = 0 so no effect, and SD, plausible range)
# It is how stronlgy seed production changes with elevation if beta is >0 then the seed production increase otherwise descrease. It is the strength of my effect
#What is my effect? It is the strength at with my response changes when a predictor changes by one unit. 
summary(psme_data$elev_sc)
range(psme_data$elev_sc) #-1.573070  2.021847
#This can help me define the range for my beta prior
#As I know that my seed produciton decrease with elevation 
lm_fit <- lm(log(total_viable_sds + 1) ~ elev_sc, data = psme_data)
summary(lm_fit)
#As I did not check for the validity of my lm but just used it to have a sort of reference about the decrease of seed production, I will use a negative slope of 0,5 and a SD of 0,5 (maybe too large?).Shoudl I define this more clearly? Because it doesn't seem linear so using an lm is for sure wrong. I don't know.. 


#Phi. (the overdispersion). 
#phi is the dispersion parameter that allows my variance to be larger than my mean. in a NB, Var(y)= mean (y) + mean(y)^2/phi. 
#Type of distribution for the phi parameter (dispersion)
p2
curve(dexp(x,1),0,5)        # exponential
curve(dlnorm(x,-1,0.7),0,5) # lognormal
curve(dgamma(x,2,1),0,5)    # gamma
#Based on the curves, it would make the most sense to use a exp. or maybe lognorm 
#exponential use a rate ^ and I can use a moderate one (1)or a slighty wider prior but than the rate is lower (0,5). This would account for the higly variable counts. 
#However, I tried and it doesn't work well so I will  phi ~ gamma(2, 2)
#It still doesn't work so I could potentially remove the 0 inflation. --> Did not work
#I tried with lognormal(-1, 0.5)

#Theta. (the zero inflated)
#theta is the probability of a structural zero. It accounts for the excta 0's. and is bounded between 0 an 1. 
#Because I know that I have 50% of my data that is 0, this means that I have a and b that are symmetric so using (5,5) means that I have 0,5 for my theta prior which reflects my data.
#Type of distribution for the theta parameter (probability of 0's)
curve(dbeta(x,1,1), 0,1)   # uniform, mean 0.5, max variance
curve(dbeta(x,2,2), 0,1)   # symmetric, mean 0.5, less spread
curve(dbeta(x,5,5), 0,1)   # symmetric, mean 0.5, tighter
curve(dbeta(x,10,5), 0,1)  # skewed, mean 0.667
#I tried that but doesn't work because my 0's are not distributed equally. I have the high sites with 0 seed prodcution so I decided to create a parameter that depends on the elevation where I have gamma0 (probability of 0) and gamma1 (effect of elevation) that control the mean and the slope and because I know that there are more 0's with elevations. Gamma0 is the baselline  and the gamma 1 is the slope 
#theta is now : inv_log(gamma0 + gamma1 *elev), where the inv_logit turns numbers nto a probability between 0 and 1 
#So gamma0 ~normal (-2,1) where inv_logit(-2)~0,12 so low probability of 0's at low elevation and SD 1 so allows uncertainity
# Gamma 1~normal(2,1) which 2 indicating a steep increase in 0's (maybe this can be changes too, into 1 see plot afterwards) and again 1 for allowing uncertainity 
#Checking it visually
gamma0 <- -2
gamma1 <- 2
elev <- seq(-2, 2, length.out=100)
theta <- plogis(gamma0 + gamma1*elev)
plot(elev, theta, type='l', ylab='θ', xlab='Elevation', main='Prior θ curve')
#looking at the 0'S from my data
zero_curve <- psme_data %>%
  group_by(elev_sc) %>%
  summarise(
    n = n(),
    n_zero = sum(total_viable_sds == 0),
    prop_zero = mean(total_viable_sds == 0)
  )
ggplot(zero_curve, aes(x = elev_sc, y = prop_zero)) +
  geom_point(size=3, color="red") +      
  geom_line(color="blue", size=1) + 
  geom_smooth(method="glm", method.args = list(family = "binomial"),
              se=FALSE, color="darkgreen", linetype="dashed") +
  ylim(0,1) +
  xlab("Scaled Elevation") +
  ylab("Proportion of zeros") +
  ggtitle("Observed proportion of zero seed counts by elevation") +
  theme_minimal()
#Looks quite similar maybe?


###STEP 4 : Prior predictive check. 
#If simulated count are way to high or all zeros than my priors need adjusting
#Once my prior predictive check looks good I can fit my model with the real data 

#First I need to create a list for the priors that I call ppc
stan_data_ppc <- list(
  N = nrow(psme_data),
  y = rep(0, nrow(psme_data)),  # dummy values instead of the real ones, creates a vecto of 0 to satisfy STAN data requirement. 
  elev = psme_data$elev_sc
)

#This is getting what I wrote in the .rstan file
# Translating my .stan code into C++ code and then compiling the C++ into a standalone executable and then saving it so it can run the MCMC sampling. 
mod <- cmdstan_model("psme_zinb.stan")
#alternative without zero inlfation
#mod <- cmdstan_model("psme_nb.stan")

fit_ppc <- mod$sample(
  data = stan_data_ppc,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500, #can be modified depending on how many times I want to do
  iter_sampling = 500
)

print(fit_ppc, variables = c("alpha", "beta", "phi", "gamma0", "gamma1"))
#print(fit_ppc, variables = c("alpha", "beta", "phi"))


#How to interpet this : 
#mad is the median absolute deviation of posterior draws. It measrures spread and sometimes reported as an alternative to SD. Is it important? 
#q5 and q95 are the 5th and 95th percentiles of the posterior draws they form the 90% interval where most of the parameter values lie given the data. 
#rhat is the gelman-rubin convergence diagnostic and compares within chain vs between chain variance (1= chains converge well) 
#Ess-bulk is the effective sample size and measures how many independent draws your posterior is equivalent to for the main bulk of the distribution (not sure I understand). High ESS is reliable estimates, low ess is high autocorrelation. 
#ESS-tail: effective sample for the tail of the posterior, important for credible intervals (especially if the distribution is skewdwd). Low ess tail is uncertainity in the extremes. 


#Prior predictive check.
#I have to first extract y_sim (from the simulated data) and then convert it to a 2D matrix (iterations*observations) because Stan give (chains*draws*N). I then select 200 draws for plotting and then plot to see it the priors predicts reasonabily 
# Extract prior predictive draws
y_sim <- fit_ppc$draws("y_sim")
# Extract y_sim as 2D 
yrep <- as_draws_matrix(fit_ppc$draws("y_sim"))
yrep <- yrep[1:200, ] # 200 draws for plotting
# Plot density overlay
ppc_dens_overlay(y = psme_data$total_viable_sds, yrep = yrep) +
  ggtitle("Prior Predictive Check: PSME seed counts") +
  xlab("Total viable seeds") +
  ylab("Density") +
  coord_cartesian(xlim = c(0, 50))   # optional: zoom on realistic seed range


###STEP 5: Fitting the model to the real data
#Stan will estimate a posterior distribution
#Stan uses both priors an dlikelihood from real data. 

#list for posterior draws
stan_data <- list(
  N = nrow(psme_data),
  y = psme_data$total_viable_sds,
  elev = psme_data$elev_sc
)

fit <- mod$sample(
  data = stan_data,   # this includes real y
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

#Model output
fit$summary()
fit$cmdstan_diagnose()
#Not good... still too many 0's

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


# Hurdle-model ------------------------------------------------------------
#As the zinb did not really work and I'm not sure if it comes from me or from the fact that this is not a good model, I will try the hurdle model now. 
#What is a hurdle model?
#the hurdle model is constructing two likelihood. One for the probability that it is 0 and the other one for when the probability is not 0 what is the probability?
#STEP1:
#First I would like then to see if there is a stand and year effect on seed production
#How much of seed production is explained by year and stand. 
# Potential step afterwards, 
#Adding time without hard constraints,adding a lag effect, was last year a mast or not?
#Is there evidence that last year reproduction affected this years reproduction?
#is this strong, weak or zero?


#Year: I will treat years as random effect, why ? because if I treat each year as a fixed effect, I'm giving each year its own parameter. and then the model treats each year as completely seperated from others. And then if one year has too few sample, the estimate for this year can be extreme. So I will have year as a random effect (partial pooling). What does this do? Random effect assume that each year comes from a common distribution, all year share information and extreme years are shrunk towards the overall mean. This prevents overfitting 

#Defining the model components: 
#Part 1 : Occurence (0 vs >0)
#Occurence (yi)~Bernouilli yi
#logit (yi)= alpha + alpha stand [i] + delta year [i]+ gamma*Seeds (i,t-1)
#Alpha = overal baseline probability
#Alpha stand [i]= random effect for the stand 
#delta year= random effect for the year
# gamma *Seeds i,t-1 )= effect of previous year seed production
#Part 2: Count part (how many seeds |>0 ; when production is not 0). 
#yi count ~negbinomial (mean i and standard deviation and dispersion parameter)
#log (meani)= beta + beta stand [i]+ n year[i]
#beta = baseline mean
#beta stand [i]= stand random effect
#beta_year[i]= year random effect


psme_data <- psme_data %>%
  arrange(stand, year) %>%
  group_by(stand) %>%
  mutate(lag_seeds = lag(total_viable_sds, 1)) %>%
  ungroup()


# prior predictive checks
stan_data_ppc <- list(
  N = nrow(psme_data),                    # number of observations
  Seeds = rep(0, nrow(psme_data)),        # dummy response for Stan
  N_Stand = length(unique(psme_data$stand)),
  Stand = as.integer(factor(psme_data$stand)),  
  N_Year  = length(unique(psme_data$year)),
  Year = as.integer(factor(psme_data$year)),    
  lag_Seeds = ifelse(is.na(psme_data$lag_seeds), 0, psme_data$lag_seeds)  
)

mod <- cmdstan_model("hurdle.stan")

fit_ppc <- mod$sample(
  data = stan_data_ppc,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500, #can be modified depending on how many times I want to do
  iter_sampling = 500
)
#Model converging issues-- need to work on that. I think changing my priors and maybe thinking about what could potentially affect my model. Because overproduction of 0's but better than the previous one I tried :) 0.99 before 1 :p

fit_ppc$summary()

#Checking if my priors are good. 

#1
prior_pred <- fit_ppc$draws("Seeds_prior")          # draws x obs x chains
prior_pred_mat <- as_draws_matrix(prior_pred)       # combine chains

#2
# Observed zeros (comparing my data with what the priors are producing)
prop_zero_obs <- mean(psme_data$total_viable_sds == 0)
prop_zero_prior <- mean(prior_pred_mat == 0)

cat("Observed proportion of zeros:", round(prop_zero_obs, 2), "\n")
cat("Prior predictive proportion of zeros:", round(prop_zero_prior, 2), "\n")

#3
# Only positive counts
positive_obs <- psme_data$total_viable_sds[psme_data$total_viable_sds > 0]
positive_seeds <- prior_pred_mat[prior_pred_mat > 0]

summary(positive_obs)
summary(positive_seeds)

#4
# define inv_logit
inv_logit <- function(x) 1 / (1 + exp(-x))

# example prior draws
alpha_occ_draws <- rnorm(1000, 0, 5)   # prior for intercept
gamma_lag_draws <- rnorm(1000, 0, 2)   # prior for lag effect

# probability of reproducing with lag = 0 and lag = 10
lag0 <- inv_logit(alpha_occ_draws + gamma_lag_draws * 0)
lag10 <- inv_logit(alpha_occ_draws + gamma_lag_draws * 10)

summary(lag0)
summary(lag10)

#5
#Checking the observed vs prior results. 
# Combine data
df_all <- data.frame(
  Seeds = c(psme_data$total_viable_sds, prior_pred_mat),
  Type = c(rep("Observed", nrow(psme_data)),
           rep("Prior", length(prior_pred_mat)))
)

# Plot using density
ggplot(df_all, aes(x = Seeds + 1, fill = Type, color = Type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
  scale_x_continuous(trans = "log1p") +  # log(1 + x) to show zeros
  labs(title = "Observed vs Prior Predictive Seeds (density, log scale)",
       x = "log1p(Seeds)", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Observed" = "tomato", "Prior" = "skyblue")) +
  scale_color_manual(values = c("Observed" = "tomato", "Prior" = "skyblue"))

# ##Posterior predictive check
# stan_data <- list(
#   N = nrow(species_data),
#   Seeds = species_data$Seeds,          # can be dummy for prior predictive
#   N_Stand = length(unique(species_data$Stand)),
#   Stand = species_data$Stand,
#   N_Year  = length(unique(species_data$Year)),
#   Year = species_data$Year,
#   lag_Seeds = species_data$lag_Seeds
# )
# 




# HMM Model AR(1) ---------------------------------------------------------------

#A hidden markov model is a model for situations where the system moves through a sequence of states over time, those states are not directly obsevable (hidden) and they generate observable data that I can measure. It is about observing the outcome not the true underlying process. 
#there is a process happening in the background (hidden states), and this process follows rules: the next state depends only on the current state (markov part) and each hidden state produces observations with certain probabilities. 

#With the data I have I have a few issues. 1) a lot of zeros and 2) occasional huge mast years 3) strong temporal autocorrelation because of the trees capacity to produce seeds so in that context my hidden state would be the tree's physiological/reproductive state?

#There could be three hidden states: low, or no repoduction (recovery)/ moderate reproduction/ mast year. But those states are latent and I can't observe then directly. 

#Difference with a hurdle model: a hurdle model explains how zeros vs non-zeros happen, why HMM explains why years cluster into low and  mast phase over time. The hurdle model splits the process into two parts (whether any seeds are produced; usually a Bernoulli(logistic model) and how many seeds are produced given production is >0; poisson or negative binomial). 

#In the hurdle : 0 = no reproduction this year and positives = reproduction occured but it doesn't encode memory, a mast year has no effect except if you force it. (I can use a lagged seed production)

#In the HMM, it assumes that the tree/stands switch among latent reproductive states and the states persists over time. 

#How to compare both of my models. THe Hurdle is going to try to explain the 0's by dispersal limitation, trap geometry, species absence and local stochasticity-- no discrete reproductive phases whereas the HMM the 0's happen because the system is sometimes in a low reproductive state and sampling (the system switches between discrete reproductive phases). 
#How to compare both models: look at the posteriror predictive check. Are both models: able to reproduce the lengths of runs of zeros, the amplitude of high years, the timing of peaks? I can also use the predictive performance (LOO/WAIC): large improvement-->structure matters/ small or no improvement : simpler model wins. 





stan_data <- list(
  N = nrow(psme_data),
  S = length(unique(psme_data$stand)),
  J = length(unique(psme_data$trapno)),
  T = length(unique(psme_data$year)),
  stand = as.integer(as.factor(psme_data$stand)),
  trap  = as.integer(as.factor(psme_data$trapno)),
  year  = as.integer(as.factor(psme_data$year)),
  y     = as.integer(psme_data$total_viable_sds),
  log_offset = log(psme_data$size)
)

stan_data_tshe <- list(
  N = nrow(tshe_data),
  S = length(unique(tshe_data$stand)),
  J = length(unique(tshe_data$trapno)),
  T = length(unique(tshe_data$year)),
  stand = as.integer(as.factor(tshe_data$stand)),
  trap  = as.integer(as.factor(tshe_data$trapno)),
  year  = as.integer(as.factor(tshe_data$year)),
  y     = as.integer(tshe_data$total_viable_sds),
  log_offset = log(tshe_data$size)
)


mod <- cmdstan_model("ar1.stan")

#I think I might need to change the prior acordingly for each species--how can I do that for a multipspecie model then. 

fit_prior <- mod$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 0,
  iter_sampling = 1000,
  fixed_param = TRUE  # VERY IMPORTANT: sample from priors only
)

y_prior_draws <- fit_prior$draws("y_prior")
y_prior_mean <- apply(as_draws_matrix(y_prior_draws), 2, mean)  # mean per trap/year

hist(y_prior_mean, breaks = 20,
     main = "Prior predictive distribution of seed counts",
     xlab = "Simulated seed counts",
     xlim=c(0,100))
abline(v=c(0,6), col="red", lty=2)  # reference range from the summary

summary(psme_data$total_viable_sds)
summary(tshe_data$total_viable_sds)


fit <- mod$sample(
  data = stan_data_tshe,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 12
)

fit$summary(c("rho", "sigma_eta", "sigma_stand", "sigma_trap"))
fit$summary(c("eta"))
fit$summary()


# Extract y_rep draws
y_rep_draws <- fit$draws("y_rep")  

y_rep_matrix <- as_draws_matrix(y_rep_draws)  # necessary for plotting
dim(y_rep_matrix) 

y_obs <- stan_data$y
y_rep_flat <- as.numeric(y_rep_matrix)

ggplot() +
  geom_histogram(aes(x = y_rep_flat, y = after_stat(density)),
                 bins = 50, fill = "skyblue", alpha = 0.5) +
  geom_histogram(aes(x = y_obs, y = after_stat(density)),
                 bins = 50, fill = "red", alpha = 0.5) +
  labs(x = "Seed counts", y = "Density",
       title = "Posterior predictive check") +
  theme_minimal()

ggplot() +
  geom_histogram(aes(x = y_rep_flat, y = ..density..),
                 bins = 50, fill = "skyblue", alpha = 0.5) +
  geom_histogram(aes(x = y_obs, y = ..density..),
                 bins = 50, fill = NA, color = "red", size = 1.2) +
  labs(x = "Seed counts", y = "Density",
       title = "Posterior predictive check") +
  theme_minimal()

eta_draws <- fit$draws("eta")  # draws × stands × years

length(unique(psme_data$stand))
length(unique(psme_data$year))


# Flatten chains and iterations into one "draws" dimension, necessary for plotting
eta_array <- as_draws_matrix(eta_draws)  

# Number of stands and years
S <- 11
T <- 16


# changing the columns names
eta_cols <- grep("^eta\\[", colnames(eta_array))
eta_mat <- eta_array[, eta_cols]  

# Convert to tidy long format for plotting
eta_long <- as.data.frame(eta_mat) %>%
  mutate(draw = 1:n()) %>%  
  pivot_longer(cols = -draw, names_to = "eta_var", values_to = "eta") %>%
  separate(eta_var, into = c("tmp", "stand", "year"), sep = "\\[|,|\\]", convert = TRUE) %>%
  select(-tmp)

# Compute posterior mean per stand/year
eta_mean <- eta_long %>%
  group_by(stand, year) %>%
  summarise(eta = mean(eta), .groups = "drop")

ggplot(eta_mean, aes(x = year, y = eta, color = factor(stand), group = stand)) +
  geom_line(linewidth =1) +
  theme_minimal() +
  labs(y = "Latent reproductive effort (eta)", x = "Year", color = "Stand")





# HMM model Matrix --------------------------------------------------------




# Aggregate seeds per stand per year
psme_sy <- psme_data %>%
  group_by(stand, year) %>%
  summarise(y = sum(total_viable_sds), .groups = "drop") %>%
  mutate(
    stand_id = as.integer(as.factor(stand)),
    year_id  = as.integer(as.factor(year))
  )

S <- length(unique(psme_sy$stand_id))
T <- length(unique(psme_sy$year_id))

# Create matrix: rows=stands, columns=years
y_mat <- matrix(0, nrow = S, ncol = T)
for(i in 1:nrow(psme_sy)) {
  y_mat[psme_sy$stand_id[i], psme_sy$year_id[i]] <- psme_sy$y[i]
}

stan_data <- list(
  S = S,
  T = T,
  y = y_mat
)


mod <- cmdstan_model("hmm.stan")

fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

fit$summary()



stand_ids <- sort(unique(psme_data$stand))
year_ids  <- sort(unique(psme_data$year))

S <- length(stand_ids)
T <- length(year_ids)

# Create y[S,T] matrix
y_obs <- matrix(NA, nrow = S, ncol = T)
for(s in 1:S){
  for(t in 1:T){
    y_obs[s,t] <- sum(psme_data$total_viable_sds[
      psme_data$stand == stand_ids[s] & psme_data$year == year_ids[t]
    ])
  }
}


n_samples <- 200

mu_state_draws <- as_draws_matrix(fit$draws(variables = "mu_state"))
phi_draws      <- as_draws_matrix(fit$draws(variables = "phi"))
P_draws        <- as_draws_matrix(fit$draws(variables = "P"))
pi_draws       <- as_draws_matrix(fit$draws(variables = "pi"))


set.seed(123)
y_rep <- array(NA, dim = c(n_samples, S, T))

for(draw in 1:n_samples){
  mu  <- mu_state_draws[draw, ]       # log mean per state
  phi <- phi_draws[draw, 1]           # overdispersion
  P   <- matrix(P_draws[draw, ], 2, 2, byrow = TRUE)  # transition matrix
  pi  <- pi_draws[draw, ]             # initial state probabilities
  
  for(s in 1:S){
    # Simulate latent states
    z <- integer(T)
    z[1] <- sample(1:2, 1, prob = pi)
    if(T > 1){
      for(t in 2:T){
        z[t] <- sample(1:2, 1, prob = P[z[t-1], ])
      }
    }
    
    # Simulate observed counts
    for(t in 1:T){
      y_rep[draw, s, t] <- rnbinom(1, mu = exp(mu[z[t]]), size = phi)
    }
  }
}

y_obs_vec <- as.vector(y_obs)
y_rep_vec <- as.vector(y_rep)

df <- data.frame(
  count = c(y_obs_vec, y_rep_vec),
  type = rep(c("Observed", "Simulated"), c(length(y_obs_vec), length(y_rep_vec)))
)


ggplot(df, aes(x = count + 1, fill = type)) +  # +1 to avoid log(0)
  geom_density(alpha = 0.5) +
  scale_x_log10() +  # optional for heavy-tailed counts
  theme_minimal() +
  labs(
    x = "Seed count (log scale)",
    y = "Density",
    title = "Posterior Predictive Check for HMM"
  )

fit$summary("mu_state")



# hmm_zeroinflated --------------------------------------------------------

# # Encode categorical variables as integers
# df <- psme_data %>%
#   mutate(
#     year_id  = as.integer(factor(year)),
#     stand_id = as.integer(factor(stand)),
#     trap_id  = as.integer(factor(trapno))
#   )
# 
# # Total counts
# N       <- nrow(df)
# T_years <- length(unique(df$year_id))
# n_stand <- length(unique(df$stand_id))
# n_trap  <- length(unique(df$trap_id))
# 
# # Build data list for Stan
# stan_data <- list(
#   N       = N,
#   T       = T_years,
#   y       = df$total_viable_sds,
#   n_stand = n_stand,
#   n_trap  = n_trap,
#   year_id = df$year_id,
#   stand_id = df$stand_id,
#   trap_id  = df$trap_id
# )
# 
# mod <- cmdstan_model("hmm_zero.stan")
# 
# fit <- mod$sample(
#   data = stan_data,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 2000,
#   seed = 123
# )
# 
# fit$summary()
# 
# library(cmdstanr)
# library(dplyr)
# library(ggplot2)
# 


# params <- as_draws_df(fit$draws())
# posterior_draws <- params[1:200,]  # use first 200 draws for speed
# 
# T_years <- length(unique(df$year_id))
# N <- nrow(df)
# 
# # Pre-allocate arrays
# all_state_probs <- array(0, dim = c(T_years, 2, nrow(posterior_draws)))
# posterior_y_pred <- matrix(0, nrow = N, ncol = nrow(posterior_draws))
# 

# compute_state_probs <- function(y, year_id, stand_id, trap_id,
#                                 alpha, beta_stand, u_trap,
#                                 init_state, trans, phi, T_years) {
#   
#   log_emission <- matrix(0, nrow = T_years, ncol = 2)
#   for (t in 1:T_years) {
#     for (s in 1:2) {
#       inds <- which(year_id == t)
#       eta <- alpha[s] + beta_stand[stand_id[inds]] + u_trap[trap_id[inds]]
#       mu <- exp(eta)
#       log_emission[t,s] <- sum(dnbinom(y[inds], size=phi, mu=mu, log=TRUE))
#     }
#   }
#   
#   # Forward pass
#   log_alpha <- matrix(0, nrow = T_years, ncol = 2)
#   log_alpha[1,] <- log(init_state) + log_emission[1,]
#   for (t in 2:T_years) {
#     for (s in 1:2) {
#       temp <- log_alpha[t-1,] + log(trans[,s])
#       log_alpha[t,s] <- log_emission[t,s] + log(sum(exp(temp - max(temp))) + 1e-12) + max(temp)
#     }
#   }
#   
#   # Backward pass
#   log_beta <- matrix(0, nrow = T_years, ncol = 2)
#   for (t in (T_years-1):1) {
#     for (s in 1:2) {
#       temp <- log(trans[s,]) + log_emission[t+1,] + log_beta[t+1,]
#       log_beta[t,s] <- log(sum(exp(temp - max(temp))) + 1e-12) + max(temp)
#     }
#   }
#   
#   # Posterior probabilities
#   state_prob <- matrix(0, nrow = T_years, ncol = 2)
#   for (t in 1:T_years) {
#     temp <- log_alpha[t,] + log_beta[t,]
#     temp <- temp - max(temp)
#     probs <- exp(temp) / sum(exp(temp))
#     state_prob[t,] <- probs
#   }
#   
#   return(state_prob)
# }
# 

# for (i in 1:nrow(posterior_draws)) {
#   
#   # Extract parameters
#   alpha <- as.numeric(posterior_draws[i, c("alpha[1]", "alpha[2]")])
#   beta_stand <- as.numeric(posterior_draws[i, grep("^beta_stand", colnames(posterior_draws))])
#   u_trap <- as.numeric(posterior_draws[i, grep("^u_trap", colnames(posterior_draws))])
#   init_state <- as.numeric(posterior_draws[i, c("init_state[1]", "init_state[2]")])
#   trans <- matrix(as.numeric(posterior_draws[i, c("trans[1,1]","trans[1,2]","trans[2,1]","trans[2,2]")]), nrow=2, byrow=TRUE)
#   phi <- as.numeric(posterior_draws[i, "phi"])
#   
#   # Compute posterior state probabilities
#   all_state_probs[,,i] <- compute_state_probs(
#     y = df$total_viable_sds,
#     year_id = df$year_id,
#     stand_id = df$stand_id,
#     trap_id = df$trap_id,
#     alpha = alpha,
#     beta_stand = beta_stand,
#     u_trap = u_trap,
#     init_state = init_state,
#     trans = trans,
#     phi = phi,
#     T_years = T_years
#   )
#   
#   # Posterior predictive simulation
#   for (n in 1:N) {
#     t <- df$year_id[n]
#     state_probs <- all_state_probs[t,,i]
#     z <- sample(1:2, size=1, prob = state_probs)
#     eta <- alpha[z] + beta_stand[df$stand_id[n]] + u_trap[df$trap_id[n]]
#     mu <- exp(eta)
#     posterior_y_pred[n,i] <- rnbinom(1, size=phi, mu=mu)
#   }
# }
# 

# summary_probs <- data.frame(
#   year = 1:T_years,
#   mean = apply(all_state_probs[,2,], 1, mean),
#   lower = apply(all_state_probs[,2,], 1, quantile, 0.025),
#   upper = apply(all_state_probs[,2,], 1, quantile, 0.975)
# )
# 
# ggplot(summary_probs, aes(x = year, y = mean)) +
#   geom_line(color = "darkgreen", size=1) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill="lightgreen") +
#   ylab("Posterior P(High Production)") +
#   xlab("Year") +
#   theme_minimal()
# 

# ppc_summary <- data.frame(
#   observed = df$total_viable_sds,
#   mean_pred = rowMeans(posterior_y_pred),
#   lower = apply(posterior_y_pred, 1, quantile, 0.025),
#   upper = apply(posterior_y_pred, 1, quantile, 0.975)
# )
# 
# ggplot(ppc_summary, aes(x = observed, y = mean_pred)) +
#   geom_point(alpha = 0.6) +
#   geom_errorbar(aes(ymin=lower, ymax=upper), alpha=0.3) +
#   geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
#   xlab("Observed seed counts") +
#   ylab("Posterior predictive mean ± 95% CI") +
#   theme_minimal()
# 
# ggplot(ppc_summary, aes(x = observed, y = mean_pred)) +
#   geom_point(color = "red", alpha = 0.6) +
#   geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
#   xlab("Observed seed counts") +
#   ylab("Predicted seed counts") +
#   theme_minimal()


# Newmodel_HMM ------------------------------------------------------------

psme_data_simple2<-psme_data%>%
  group_by(year,stand) %>%
  summarise(y = sum(total_viable_sds), .groups = "drop")

head(psme_data_simple2)

# Pivot to a stand x year matrix
y_mat <- psme_data_simple2 %>%
  pivot_wider(
    names_from = year,
    values_from = y,
    values_fill = 0   # fill missing years with 0 seeds
  ) %>%
  arrange(stand) %>%     # optional, sort stands
  select(-stand) %>%     # remove stand column for matrix
  as.matrix()            # convert to numeric matrix




stand_ranges <- psme_data_simple2 %>%
  group_by(stand) %>%
  summarise(
    start_year = min(year),
    end_year   = max(year),
    .groups = "drop"
  )

stand_ranges <- stand_ranges %>%
  mutate(
    year_seq = purrr::map2(start_year, end_year, seq)
  )


missing_years <- psme_data_simple2 %>%
  group_by(stand) %>%
  summarise(
    start_year = min(year),
    end_year   = max(year),
    observed_years = list(sort(unique(year))),
    full_seq       = list(seq(start_year, end_year)),
    missing_years  = list(setdiff(full_seq[[1]], observed_years[[1]])),
    .groups = "drop"
  )

str(missing_years)

# Sort data properly
data_sorted <- psme_data_simple2 %>%
  arrange(stand, year)

# Create stacked observation vector
y <- data_sorted$y   # <-- define y here

# Count number of observations per stand
T_i <- data_sorted %>%
  count(stand) %>%
  pull(n)            # <-- define T_i here

# Define dimensions
F <- length(T_i)     # number of stands
N <- length(y)       # total number of observations


start_idxs <- c()
end_idxs <- c()
id <- 1
for(s in 1:F){
  start_idxs <- c(start_idxs, id)
  id <- id + T_i[s]-1
  end_idxs <- c(end_idxs, id)
  id <- id + 1
}


data_list <- list(
  F   = F,
  N   = N,
  start_idxs = start_idxs, 
  end_idxs = end_idxs, 
  y   = y
)


mod <- cmdstan_model("hmm_stand.stan")

fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4
)

fit$summary()
summary(fit$draws(c("theta","log_mu")))


yrep <- fit$draws("y_rep")

yrep_mat <- as_draws_matrix(yrep)

summary(as.vector(yrep_mat))

hist(as.vector(yrep_mat), breaks = 50)

plot(density(as.vector(yrep_mat)))

boxplot(yrep_mat[,1:20], outline = FALSE)

df <- data.frame(y_pred = as.vector(yrep_mat))

ggplot(df, aes(x = y_pred)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  theme_minimal()


# 1. Extract posterior predictive draws
yrep_mat <- as_draws_matrix(fit$draws("y_rep"))

# 2. Create stand index using start/end indices
stand_id <- rep(NA, data_list$N)

for (f in 1:data_list$F) {
  stand_id[data_list$start_idxs[f]:data_list$end_idxs[f]] <- f
}

lookup <- data.frame(
  idx = 1:data_list$N,
  stand = stand_id
)

# 3. Convert to long format
yrep_long <- as.data.frame(yrep_mat)
colnames(yrep_long) <- 1:data_list$N

yrep_long <- yrep_long %>%
  mutate(draw = 1:n()) %>%
  pivot_longer(
    cols = -draw,
    names_to = "idx",
    values_to = "y_pred"
  ) %>%
  mutate(idx = as.integer(idx)) %>%
  left_join(lookup, by = "idx")

# 4. Add observed data with same stand indexing
obs_df <- data.frame(
  idx = 1:data_list$N,
  y = data_list$y,
  stand = stand_id
)

# 5. Plot posterior predictive densities per stand
ggplot() +
  geom_density(
    data = yrep_long,
    aes(x = y_pred),
    fill = "lightblue",
    alpha = 0.5
  ) +
  geom_density(
    data = obs_df,
    aes(x = y),
    color = "black",
    size = 1
  ) +
  facet_wrap(~stand, scales = "free") +
  theme_minimal() +
  labs(
    x = "Seed count",
    y = "Density",
    title = "Posterior Predictive Densities per Stand"
  )








# # Convert to a matrix: (iterations * chains) × T
# yrep_mat <- as_draws_matrix(yrep)  # posterior package handles flattening
# 
# yrep_quants <- apply(yrep_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))

# years <- psme_data_simple2$year
# y_obs <- psme_data_simple$y
# 
# df_plot <- data.frame(
#   year = years,
#   y_obs = y_obs,
#   lower = yrep_quants[1, ],
#   median = yrep_quants[2, ],
#   upper = yrep_quants[3, ]
# )
# 
# 
# ggplot(df_plot, aes(x = year)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
#   geom_line(aes(y = median), color = "blue", size = 1) +
#   geom_point(aes(y = y_obs), color = "red", size = 2) +
#   labs(
#     y = "Seed counts",
#     x = "Year",
#     title = "Posterior Predictive Check"
#   ) +
#   theme_minimal()
# 
# #changer avec deux distributions plot : avec les données. 
# 
# # Suppose your fitted model is 'fit'
# posterior <- fit$draws(variables = c("mu", "phi"), format = "df")
# 
# # Check the first few rows
# head(posterior)
# 
# 
# 
# set.seed(123)
# 
# # Take a subset of posterior draws to reduce plotting size
# posterior_sub <- posterior[sample(nrow(posterior), 500), ]
# 
# # Simulate distributions
# sim_data <- data.frame()
# 
# for (i in 1:nrow(posterior_sub)) {
#   # State 1 (low)
#   sim_low <- rnbinom(100, size = posterior_sub$`phi[1]`[i],
#                      mu = posterior_sub$`mu[1]`[i])
#   sim_data <- rbind(sim_data,
#                     data.frame(seed = sim_low, state = "low"))
#   
#   # State 2 (high)
#   sim_high <- rnbinom(100, size = posterior_sub$`phi[2]`[i],
#                       mu = posterior_sub$`mu[2]`[i])
#   sim_data <- rbind(sim_data,
#                     data.frame(seed = sim_high, state = "high"))
# }
# 
# 
# ggplot(sim_data, aes(x = seed, fill = state)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(name = "Seed Production") +
#   scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
#   theme_minimal() +
#   ggtitle("Posterior Predictive Distributions of Seed Production per State")
# 
# 
# mu_means <- colMeans(posterior_sub[, c("mu[1]", "mu[2]")])
# 
# ggplot(sim_data, aes(x = seed, fill = state)) +
#   geom_density(alpha = 0.5) +
#   geom_vline(xintercept = mu_means[1], color = "blue", linetype = "dashed") +
#   geom_vline(xintercept = mu_means[2], color = "red", linetype = "dashed") +
#   scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
#   theme_minimal() +
#   xlim(c(0,300))
# 
# ggplot(psme_data_simple, aes(x = y)) +
#   geom_density(fill = "skyblue", alpha = 0.5) +
#   ggtitle("Density")
# 
# 
# ggplot() +
#   # Posterior predictive densities (by state)
#   geom_density(data = sim_data, aes(x = seed, fill = state), alpha = 0.5) +
#   geom_vline(xintercept = mu_means[1], color = "blue", linetype = "dashed") +
#   geom_vline(xintercept = mu_means[2], color = "red", linetype = "dashed") +
#   scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
#   # Actual observed data density
#   geom_density(data = psme_data_simple, aes(x = y),
#                fill = "black", alpha = 0.3) +
#   theme_minimal() +
#   xlim(c(0, 300)) +
#   ggtitle("Posterior Predictive Densities (by state) with Observed Data")
