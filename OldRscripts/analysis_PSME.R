library(dplyr)
library(ggplot2)
library(rstan)



getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

#Defining the specie
species <- "PSME"

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
min(species_data$total_viable_sds); max(species_data$total_viable_sds);var(species_data$total_viable_sds);median(species_data$total_viable_sds)
min(stand_year_df$y); max(stand_year_df$y);var(stand_year_df$y);median(stand_year_df$y)


fit05phi4 <- fit_hmm_model("Stan_code/Species_Stan_Model/PSME.stan", stan_data) #Based on 5) using my function defined in functions.R 
saveRDS(fit05phi4, "fit05phi4_PSME.rds")

#Changing the log_mu prior
fit025 <- readRDS("fit025_PSME.rds") # only 1 divergent transition
fit05 <- readRDS("fit05_PSME.rds") # only 1 divergent transition
fit07 <- readRDS("fit07_PSME.rds") #no divergence
fit09 <- readRDS("fit09_PSME.rds") #no divergence

fit05phi2 <- readRDS("fit05phi2_PSME.rds") #phi is same for both and normal(6.5, 3)
fit05phi3 <- readRDS("fit05phi2_PSME.rds") #phi is different normal(6.5, 1 (low)/3 (high))
fit05phi4 <- readRDS("fit05phi2_PSME.rds") #phi is different normal(6.5, 2(low)/3 (high))

#Posterior predictive check
util <- new.env()
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)


diagnostics <- util$extract_hmc_diagnostics(fit05phi4)
util$check_all_hmc_diagnostics(diagnostics)


samples <- util$extract_expectand_vals(fit05phi4)

base_samples <- util$filter_expectands(samples,
                                       c("rho","theta1","theta2","log_lambda","log_mu","stand_effect_raw","phi1","phi2","sigma"),check_arrays = TRUE)

util$check_all_expectand_diagnostics(base_samples)

all_plots <- list()

names_yrep <- paste0("y_rep[", 1:stan_data$N, "]")




# -------------------------
# 1) Simulate priors
# -------------------------
n_prior <- 1e4

prior_df <- data.frame(
  log_lambda = rnorm(n_prior, 0, log(30)),
  log_mu     = rnorm(n_prior, log(200), 0.5),
  sigma      = rnorm(n_prior, 0, 0.5/2.57),
  phi1       = rnorm(n_prior, 6.5, 4),
  phi2       = rnorm(n_prior, 6.5, 3)
)


prior_df$draw <- 1:n_prior

# -------------------------
# 2) Extract posterior
# -------------------------
post_samples <- rstan::extract(fit05phi4)

posterior_df <- data.frame(
  log_lambda = post_samples$log_lambda,
  log_mu    = post_samples$log_mu,
  sigma      = post_samples$sigma,
  phi1       = post_samples$phi1,
  phi2       = post_samples$phi2
)


posterior_df$draw <- 1:nrow(posterior_df)

# -------------------------
# 3) Convert to long format
# -------------------------
prior_long <- reshape(
  prior_df,
  direction = "long",
  varying = list(names(prior_df)[!names(prior_df) %in% "draw"]),
  v.names = "value",
  timevar = "parameter",
  times = names(prior_df)[!names(prior_df) %in% "draw"],
  idvar = "draw"
)

posterior_long <- reshape(
  posterior_df,
  direction = "long",
  varying = list(names(posterior_df)[!names(posterior_df) %in% "draw"]),
  v.names = "value",
  timevar = "parameter",
  times = names(posterior_df)[!names(posterior_df) %in% "draw"],
  idvar = "draw"
)

# -------------------------
# 4) Plot
# -------------------------
ggplot() +
  geom_density(data = prior_long, aes(x = value, colour = "Prior"), linewidth = 0.8) +
  geom_density(data = posterior_long, aes(x = value, colour = "Posterior"), linewidth = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  labs(
    title = paste0("Prior vs Posterior (0.5/4:3) - ", species),
    x = "Parameter value",
    y = "Density",
    color = "Curve"
  ) +
  theme_minimal()



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
    display_ylim = c(0,700),
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


str(seed_data)


#1) sort the data

