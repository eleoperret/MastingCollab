#HMM model for PSME
#Analysis on seed production for one specie over elevation (stand used as proxy)
#Eléonore Perret - February 2026


#TO DO: 
#Description of the script.  
#Add titles to each plots to make sure I know what I'm seing and making sure all labels are correct
#Adding the plots Lizzie wants from what Christophe has been doing

library(dplyr)
library(ggplot2)
library(rstan)
library(posterior)
library(bayesplot)
library(tidyr)
library(reshape2)

getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

seed_data<-read.csv("SeedData_all.csv")


#Make this in an .R file so I can save it. And add this to a basic script
stand_elev <- data.frame(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13",
            "AO03","AG05","AV06","AX15","AB08","PP17",
            "AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600,668,700,850,1450,850,
                900,950,1060,1090,1100,1150,
                1150,1200,1450,1600,1700,1800)
)
saveRDS(stand_elev, "stand_elevation_table.rds")

#Figure out what to do with that
stand_order <- stand_elev %>%
  arrange(elevation) %>%
  pull(stand)

# PSME --------------------------------------------------------------------
#Changing my dataset to PSME where only stands with seeds -- removing AE10, AR07, PARA, SUNR, SPRY, AV14, AM16
#Creating something like this for each species with a loop. or storing it somewhere as a .R ?
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

summary(psme_data$total_viable_sds)
mean(psme_data$total_viable_sds) #6.48657
var(psme_data$total_viable_sds) #139.4104
#This is overdispersed Var>> mean


min(psme_data$size); max(psme_data$size);var(psme_data$size)
min(psme_data$total_viable_sds); max(psme_data$total_viable_sds);var(psme_data$total_viable_sds)


# HMM: One specie x Stands (Partial Pooled)  x Traps (as offset) x same start year and also with new transition matrix per stand ---------------------------------------

#Creating the correct dataset for my analysis.
stand_year_df <- psme_data %>%
  group_by(stand, year) %>%
  filter(year!="2009") %>% #All starting in 2010 
  summarise(
    y = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE), #Summarizing the total trap aerea
    .groups = "drop"
  ) %>%
  arrange(stand, year)


# Compute seed density (seeds per m²)
stand_year_df <- stand_year_df %>%
  mutate(seed_density = y / area)
# Seed density over time
ggplot(stand_year_df, aes(x = year, y = seed_density, color = stand, group = stand)) +
  geom_line(size = 1) +
  geom_point() +
  labs(
    title = "Seed Density (seeds/m²) per Stand over Time",
    x = "Year",
    y = "Seed Density"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

#For stan code
years_per_stand <- stand_year_df %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")

F <- nrow(years_per_stand)
T_i <- years_per_stand$T_i

start_idxs <- c()
end_idxs <- c()
id <- 1

for (s in 1:F) {
  start_idxs <- c(start_idxs, id)
  id <- id + T_i[s] - 1
  end_idxs <- c(end_idxs, id)
  id <- id + 1
}

y <- stand_year_df$y
area <- stand_year_df$area
N <- length(y)

stan_data <- list(
  F = F,              
  N = N,              
  T_i = T_i,          
  start_idxs = start_idxs,
  end_idxs = end_idxs,
  y = y,              
  area = area        
)



mod1 <- stan_model("Stan_code/feb18/hmm_2highnegbin_standpooling_TrapScaling_TransitionStand.stan")
mod2 <- stan_model("Stan_code/feb18/hmm_2highnegbin_standpooling_TrapScaling_TransitionStand2.stan")#different priors but same model (based on older models)

fit1 <- sampling(
  mod1,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)
fit2 <- sampling(
  mod2,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)

####
##Plots
#Based on Mike's code 
util<- new.env()
#creating a source file to make all the plots. 
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit1)
diagnostics <- util$extract_hmc_diagnostics(fit2)
util$check_all_hmc_diagnostics(diagnostics)

# extraire les posterior values
samples <- util$extract_expectand_vals(fit1)
samples <- util$extract_expectand_vals(fit2)

# diagnostics parametre par parametre
base_samples <- util$filter_expectands(samples,
                                       c("rho", "theta1","theta2",
                                         "log_lambda", "log_mu",
                                         "stand_effect_raw", "phi1", "phi2",                                         "sigma"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

#Ordering the stand based on elevation for PSME
original_order <- unique(stand_year_df$stand)

# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data$y)

par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:stan_data$F){
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  subsamples <-  util$filter_expectands(samples,
                                        paste0("y_rep[", idxs, "]"))
  util$plot_hist_quantiles(subsamples, "y_rep",
                           bin_min = 0, bin_max = 500, bin_delta = 10,
                           baseline_values = stan_data$y[idxs])
}


#For plotting
#Adding the stand elevation and names
perm <- match(stand_order, original_order)


#Probability of being in mast or not
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples,names,
    baseline_values = stan_data$y[idxs],
    main = paste0(stand_order[i], " (",
                  stand_elev$elevation[stand_elev$stand == stand_order[i]],
                  " m)")
  )
}

#Scaled
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples,names,                                    baseline_values = stan_data$y[idxs],display_ylim = c(0,400),main = paste0(stand_order[i], " (",       stand_elev$elevation[stand_elev$stand == stand_order[i]]," m)")
  )
}


#Parameter : **rho** (Probability in starting in a specific state)
par(mfrow = c(1,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["rho[1]"]], 50, flim=c(0,1),display_name = bquote(rho * "(Low)" ))
util$plot_expectand_pushforward(samples[["rho[2]"]], 50, flim=c(0,1),display_name = bquote(rho * "(High)" ))


#Parameter : **theta** (probability of changing state)
#General Probability of switching (transition matrix per stand)
par(mfrow = c(1,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["theta1[1]"]], 50, flim = c(0,1),display_name = bquote(theta * "(Staying in Low)" ))
util$plot_expectand_pushforward(samples[["theta2[1]"]], 50, flim = c(0,1),display_name = bquote(theta * "(Staying in Mast)" ))


#Parameter : **log_lambda** (seed production for non-mast)
par(mfrow = c(1,1), mar = c(5,5,1,1))
util$plot_expectand_pushforward(exp(samples[["log_lambda"]]), 20,
                                display_name = bquote(lambda * " (low)"))


#Parameter : **log_mu** (mean seed production for mast)
util$plot_expectand_pushforward(exp(samples[["log_mu"]]), 50, display_name = bquote(mu * " (high all stands)"))


#Parameter : log_alpha (grand mean seed production for mast including partial pooling of my stands)
#per stand
par(mfrow=c(4,3))
for(f in 1:F){
  util$plot_expectand_pushforward(
    exp(samples[[paste0("log_alpha[", f, "]")]]), 50,display_name = bquote(alpha * "(high)"))#remove the exp if wanting the log scale
}


#Parameter: **phi** (overdispersion)
par(mfrow=c(2,1))
util$plot_expectand_pushforward(samples[["phi1"]], 50, flim=c(0,1),display_name = bquote(phi * " (NB for nonmasting)"))
util$plot_expectand_pushforward(samples[["phi2"]], 50, flim=c(0,1),display_name = bquote(phi * " (NB for masting)"))



#Other things
str(Cleaned_mapping_2017_Rainier)


treestand<-Cleaned_mapping_2017_Rainier%>%
  group_by(stand_id)%>%
  summarise(species_list = list(unique(species)))
treestandpsme <- Cleaned_mapping_2017_Rainier %>%
  filter(species == "PSME") %>%
  pull(stand_id) %>%
  unique()

#Plots for Lizzie
#Run the model without the likelihood
n_prior <- 1e4

#model based on Victor's prior
prior_df <- data.frame(
  log_lambda = rnorm(n_prior, 0, log(30)/2.57),
  log_mu = rnorm(n_prior, log(200), 0.25),
  sigma = rnorm(n_prior, 0, 0.5/2.57),
  phi1 = rnorm(n_prior, 6.5, 1.79),
  phi2 = rnorm(n_prior, 6.5, 1.79)
)

#Model based on old priors
prior_df <- data.frame(
  log_lambda = rnorm(n_prior, 0, log(5)/2.57),
  log_mu = rnorm(n_prior, log(200), 0.1),
  sigma = rnorm(n_prior, 0, 0.5/2.57),
  phi1 = rnorm(n_prior, 2, 0.1),
  phi2 = rnorm(n_prior, 2, 0.1)
)

Stand <- stan_data$F
for(f in 1:Stand){
  prior_df[[paste0("stand_effect_", f)]] <- rnorm(n_prior, 0, 1)
}

post_samples <- rstan::extract(fit2)

posterior_df <- data.frame(
  log_lambda = post_samples$log_lambda,
  log_mu = post_samples$log_mu,
  sigma = post_samples$sigma,
  phi1 = post_samples$phi1,
  phi2 = post_samples$phi2
)

# add stand effects
for(f in 1:F){
  posterior_df[[paste0("stand_effect_", f)]] <- post_samples$stand_effect_raw[,f] * post_samples$sigma
}

prior_long <- reshape(
  prior_df,
  direction = "long",
  varying = list(names(prior_df)),
  v.names = "value",
  timevar = "parameter",
  times = names(prior_df),
  idvar = "draw"
)

posterior_long <- reshape(
  posterior_df,
  direction = "long",
  varying = list(names(posterior_df)),
  v.names = "value",
  timevar = "parameter",
  times = names(posterior_df),
  idvar = "draw"
)

ggplot() +
  geom_density(data = prior_long,
               aes(x = value, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = posterior_long,
               aes(x = value, colour = "Posterior", group = parameter),
               linewidth = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  labs(title = "Prior vs Posterior", x = "Parameter value", y = "Density", color = "Curve") +
  theme_minimal()
