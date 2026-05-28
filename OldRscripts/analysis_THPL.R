library(dplyr)
library(ggplot2)
library(rstan)



getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

#Defining the specie
species <- "THPL"

seed_data<-read.csv("SeedData_all.csv")


source("functions.R")


# 1) Filter for the species in general
species_data <- filter_species_data(seed_data, species)

# 2) Remove invalid stands for that species

stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R

species_data <- remove_invalid_stands(
  data = species_data,
  species_name = species,
  stands_per_species = stands_per_species
)

#only stands where the specie is present. 
unique(species_data$stand)


# Load stand elevations from saved RDS
stand_elev <- readRDS("data/stand_elevation_table.rds")#from General.Data.R
species_data <- left_join(species_data, stand_elev, by = "stand")

plot_total_seeds_per_stand(species_data)


stand_year_df <- create_stand_year_dataset(species_data) #Based on 2) using my function defined in functions.R 

plot_seed_density_timeseries(stand_year_df) #Based on 3) using my function defined in functions.R 

stan_data <- prepare_stan_data(stand_year_df) #Based on 4) using my function defined in functions.R 

#figuring out the priors
min(species_data$size); max(species_data$size);var(species_data$size)
min(species_data$total_viable_sds); max(species_data$total_viable_sds);var(species_data$total_viable_sds)

fit <- fit_hmm_model("Stan_code/Species_Stan_Model/thpl.stan", stan_data) #Based on 5) using my function defined in functions.R 


#Posterior predictive check
source("PPCfunctions.R")

F <- length(unique(stand_year_df$stand))
run_hmm_diagnostics(
  fit,
  stan_data,
  stand_year_df,
  stand_elev,
  F,
  species = species
)


#Prior predictive check
source("PriorVsPosteriorFunction.R")
n_prior <- 1e4

#TO change base on priors
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

plot_prior_posterior(
  prior_df = prior_df,
  fit = fit,
  stan_data = stan_data
)



#For non-working plots: 
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


#For plotting
# --- Prepare stand ordering ---
original_order <- unique(stand_year_df$stand)

# stand names sorted by elevation
stand_order <- stand_elev %>%
  filter(stand %in% original_order) %>%
  arrange(elevation) %>%
  pull(stand)

# perm[i] gives the index in original_order for stand_order[i]
perm <- match(stand_order, original_order)

# --- PPC per stand ---
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]   # correct index in stan_data
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  subsamples <- util$filter_expectands(samples, paste0("y_rep[", idxs, "]"))
  
  util$plot_hist_quantiles(
    subsamples, "y_rep",
    bin_min = 0,
    bin_max = 2000,
    bin_delta = 50,
    baseline_values = stan_data$y[idxs],
    main = paste0(
      stand_order[i], " (",
      stand_elev$elevation[stand_elev$stand == stand_order[i]], " m)"
    )
  )
}

# --- Probability of being in mast or not (unscaled) ---
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]  # index in stan_data
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names_yrep <- paste0("y_rep[", idxs, "]")
  
  util$plot_disc_pushforward_quantiles(
    samples, names_yrep,
    baseline_values = stan_data$y[idxs],
    main = paste0(
      stand_order[i], " (",
      stand_elev$elevation[stand_elev$stand == stand_order[i]], " m)"
    )
  )
}

# --- Probability of being in mast or not (scaled) ---
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data$F){
  s <- perm[i]  # index in stan_data
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names_yrep <- paste0("y_rep[", idxs, "]")
  
  util$plot_disc_pushforward_quantiles(
    samples, names_yrep,
    baseline_values = stan_data$y[idxs],
    display_ylim = c(0,3000),
    main = paste0(
      stand_order[i], " (",
      stand_elev$elevation[stand_elev$stand == stand_order[i]], " m)"
    )
  )
}

# --- log_alpha per stand ---
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(f in 1:stan_data$F){
  s <- perm[f]  # correct index in original data
  
  util$plot_expectand_pushforward(
    exp(samples[[paste0("log_alpha[", f, "]")]]),
    50,
    display_name = bquote(alpha * "(high)"),
    main = paste0(
      stand_order[f], " (",
      stand_elev$elevation[stand_elev$stand == stand_order[f]], " m)"
    )
  )
}
