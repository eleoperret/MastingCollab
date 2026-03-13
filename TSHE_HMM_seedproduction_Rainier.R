#HMM model for TSHE
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

# TSHE --------------------------------------------------------------------
##Should run based on a species name
#Data cleaning. 
#Add a cleaning code that decides based on the species presence or absence from the dataset here : 
# str(Cleaned_mapping_2017_Rainier)
# treestand<-Cleaned_mapping_2017_Rainier%>%
#   group_by(stand_id)%>%
#   summarise(species_list = list(unique(species)))
# treestandpsme <- Cleaned_mapping_2017_Rainier %>%
#   filter(species == "PSME") %>%
#   pull(stand_id) %>%
#   unique()

tshe_data<-seed_data%>%
  filter(spp=="TSHE")%>%
  #Added this constraint
  filter(!stand %in% c("AE10", "AR07", "PARA", "PP17","SUNR", "SPRY"))


total_stand_tshe <- tshe_data %>%
  group_by(stand) %>%
  summarise(total_seeds = sum(total_viable_sds, na.rm = TRUE))

total_stand_tshe <- total_stand_tshe %>%
  left_join(stand_elev, by = "stand") %>%
  arrange(elevation)

ggplot(total_stand_tshe, aes(x = reorder(stand,elevation), y = total_seeds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Stand",
       y = "Total viable seeds",
       title = "Total viable seeds per stand") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

summary(tshe_data$total_viable_sds)
mean(tshe_data$total_viable_sds) #6.48657
var(tshe_data$total_viable_sds) #139.4104
#This is overdispersed Var>> mean


min(tshe_data$size); max(tshe_data$size);var(tshe_data$size)
min(tshe_data$total_viable_sds); max(tshe_data$total_viable_sds);var(tshe_data$total_viable_sds)


# HMM: One specie x Stands (Partial Pooled)  x Traps (as offset) x same start year and also with new transition matrix per stand ---------------------------------------

##Should run based on a species name

#Creating the correct dataset for my analysis.
stand_year_df <- tshe_data %>%
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



mod <- stan_model("Stan_code/feb18/tshe.stan")

fit <- sampling(
  mod,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)


####
##Should run based on a species name
#Plots
#Based on Mike's code 
util<- new.env()
#creating a source file to make all the plots. 
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

# extraire les posterior values
samples <- util$extract_expectand_vals(fit)


# diagnostics parametre par parametre
base_samples <- util$filter_expectands(samples,
                                       c("rho", "theta1","theta2",
                                         "log_lambda", "log_mu",
                                         "stand_effect_raw", "phi1", "phi2",                                         "sigma"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)


# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 3000, bin_delta = 50,
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
                           bin_min = 0, bin_max = 2000, bin_delta = 50,
                           baseline_values = stan_data$y[idxs])
}


#For plotting

#Ordering the stand based on elevation for PSME
original_order <- unique(stand_year_df$stand)
stand_ordertshe <- stand_elev %>%
  filter(stand %in% original_order) %>%  # keep only stands present in original_order
  arrange(elevation) %>%
  pull(stand)

#Adding the stand elevation and names
perm <- match(stand_ordertshe, original_order)


#Probability of being in mast or not
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples,names,
                                       baseline_values = stan_data$y[idxs],
                                       main = paste0(stand_order[i], " (", stand_elev$elevation[stand_elev$stand == stand_order[i]],
                                                     " m)")
  )
}

#Scaled
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples,names, baseline_values = stan_data$y[idxs],display_ylim = c(0,3000),main = paste0(stand_order[i], " (",       stand_elev$elevation[stand_elev$stand == stand_order[i]]," m)")
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
util$plot_expectand_pushforward(samples[["phi1"]], 50,display_name = bquote(phi * " (NB for nonmasting)"))
util$plot_expectand_pushforward(samples[["phi2"]], 50,display_name = bquote(phi * " (NB for masting)"))

##Add saving the resulting plots for all posteriror check somewhere 


#Other things
str(Cleaned_mapping_2017_Rainier)


treestand<-Cleaned_mapping_2017_Rainier%>%
  group_by(stand_id)%>%
  summarise(species_list = list(unique(species)))
treestandpsme <- Cleaned_mapping_2017_Rainier %>%
  filter(species == "PSME") %>%
  pull(stand_id) %>%
  unique()

#Plots for Lizzie of the priors! 
#Run the model without the likelihood
n_prior <- 1e4

##to change based on the priors used. 
prior_df <- data.frame(
  log_lambda = rnorm(n_prior, 4,0.5),
  log_mu = rnorm(n_prior, 6,0.5),
  sigma = rnorm(n_prior, 0, 0.5/2.57),
  phi1 = rnorm(n_prior, 2, 0.1),
  phi2 = rnorm(n_prior, 2, 0.1)
)

Stand <- stan_data$F
for(f in 1:Stand){
  prior_df[[paste0("stand_effect_", f)]] <- rnorm(n_prior, 0, 1)
}

post_samples <- rstan::extract(fit)

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

##Add saving the resulting plots somewhere 