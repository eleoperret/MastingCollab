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

#I can see here that some stand have 0 seed production in 15 years (very unlikely that it will change) except for SPRY and AR07 (but to be honest I'm a bit suprised; in total it is 2 seeds found for the AR07 and 3 for the SPRY over the last years; what to do with that. I honestly cannot refute the misidentification because there are no adult tree anywhere close to this area so something to think about)


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

#Hurdle model (zero vs non-zero). Does the site produce any seeds at all? if it produces seeds how many. 
#Difference with zinb, is that in the zinb the zeros can come from two different sources whereas the hurdle only from the hurdle bever from the model. 


y <- psme_data$total_viable_sds
elev <- psme_data$elev_sc

# Indicator: 1 = produces seeds, 0 = zero
is_pos <- as.integer(y > 0)

# Positive counts only
y_pos <- y[y > 0]
elev_pos <- elev[y > 0]

stan_data <- list(
  N = length(y),
  N_pos = length(y_pos),
  y = y,
  y_pos = y_pos,
  elev = elev,
  elev_pos = elev_pos
)








