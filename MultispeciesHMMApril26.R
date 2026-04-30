##Code for multispecies HMM model 
# Stan code: 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

####Doesn't include PARA and SUNR as I do not have a map of the species. 

##Partially pooling also species
##Theta partially pooled by species and stand
##Seed production on the low and high state partially pooled also 


###STOPPED at partial pooling non-centered on both the species and stand for overdispersion
##To do : Put the results plot on github and then try the partial pooling on only the species
## Then see if better or not.
##then try the PPC plot again and maybe do a posterior vs prior check and then do a widening of the priors
##Thing about discussion with Ken 


##Next things to do : 

#1 : Add the 2009 year, handle the SUNR missing year (add PARA and SUNR)
#2 : Trap-level
#3 : Extract values for the synchrony.
#4 : Clean my code (structure of the code and clean the top part)


#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)
library(tidyverse)
library(posterior)

options(mc.cores = parallel::detectCores())


#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

seed_data<-read.csv("SeedData_all.csv")

# 1) Keep species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","ABLA","CANO","PSME","TSHE","THPL"))
unique(seed_data$stand)
#I include ABLA now

# 2) Removing invalid stands 
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)

# Count number of species per stand per year
species_count_check <- stand_year_all %>%
  group_by(stand, year) %>%
  summarise(n_species = n_distinct(spp), .groups = "drop")

# Look at raw data for one species x stand
stand_year_all %>%
  filter(spp == "ABLA", stand == "AE10") %>%
  ggplot(aes(x = year, y = y / area)) +
  geom_line() + geom_point() +
  labs(title = "Raw seed density - ABLA AE10",
       y = "Seeds per m²")

stand_year_all %>%
  mutate(density = y / area) %>%
  group_by(spp) %>%
  summarise(
    low_approx  = median(density),           # rough low state
    high_approx = quantile(density, 0.9),    # rough high state
    log_low     = log(median(density) + 0.5),  # +0.5 to avoid log(0)
    log_high    = log(quantile(density, 0.9))
  )


# Data prep ---------------------------------------------------------------

stand_year_all <- seed_filtered %>%
  filter(year != 2009) %>% #because like this they start all the same year
  group_by(spp, stand, year) %>%
  summarise(
    y    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(spp, stand, year)

str(stand_year_all)
unique(stand_year_all$area)

#Defining my stan list
# Species index
stand_year_all$species_id <- as.numeric(as.factor(stand_year_all$spp))
S <- length(unique(stand_year_all$spp)) #which is 5 at the moment

# Checking my species order before running
print(levels(as.factor(stand_year_all$spp)))

# Creating my stan data list
years_per_series <- stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands <- length(unique(years_per_series$stand_id))  # NEW

G          <- nrow(years_per_series) #More rows now because I added ABLA
T_i        <- years_per_series$T_i #different for every stands and species (56 rows)
start_idxs <- cumsum(c(1, T_i[-G])) #cumulative sum: ragged vector 
end_idxs   <- cumsum(T_i)

# Checking everything
cat("N =", nrow(stand_year_all), "\n")
cat("F =", G, "\n")
cat("S =", S, "\n")
cat("Indices correct:", all(end_idxs - start_idxs + 1 == T_i), "\n")
cat("Final index matches N:", tail(end_idxs, 1) == nrow(stand_year_all), "\n")

#Final Stan data list #added a stand parameter 
stan_data_all <- list(
  N          = nrow(stand_year_all),
  F          = G,
  S          = S,
  N_stands   = N_stands, #Nouveau
  stand_id   = years_per_series$stand_id,  # Nouveau
  sp         = stand_year_all$species_id,
  start_idxs = start_idxs,
  end_idxs   = end_idxs,
  y          = stand_year_all$y,
  area       = stand_year_all$area
)

unique(years_per_series$spp)

# Fitting Model -----------------------------------------------------------

fit_allPhi2 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_deltaPhi5.stan",
  data    = stan_data_all,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_allPhi3 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_deltaPhi6.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_allPhi4 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_deltaPhi7.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# Diagnostic plots --------------------------------------------------------


#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_allPhi4, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_allPhi4, 
      pars = c("sigma_theta1_species", "sigma_theta2_species",
               "sigma_theta1_stand", "sigma_theta2_stand",
               "sigma_low_species", "sigma_low_stand",
               "sigma_log_delta_high_species", "sigma_log_delta_high_stand"),
      condition = "accept_stat__")

#Model : Multispecies Paritial Pooling Delta Phi 6

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_allPhi4)
util$check_all_hmc_diagnostics(diagnostics)
#there is 1% of divergence

samples <- util$extract_expectand_vals(fit_allPhi3)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_species[', 1:6, ']'),
                                         paste0('alpha_theta2_species[', 1:6, ']'),
                                         paste0('alpha_theta1_stand[',   1:15, ']'),
                                         paste0('alpha_theta2_stand[',   1:15, ']'),
                                         'sigma_theta1_species', 'sigma_theta2_species',
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_species',
                                          'sigma_log_delta_high_stand',
                                          paste0('log_delta_high_species[', 1:6, ']'),
                                          
                                          # Phi by species and state (new)
                                          paste0('log_phi_species_state[', 1:6, ',1]'),
                                          paste0('log_phi_species_state[', 1:6, ',2]'),
                                          'mu_log_phi',
                                          'sigma_log_phi',
                                          
                                          # Low state means
                                          'grand_mean_low',
                                          paste0('alpha_low_species[', 1:6, ']'),
                                          'sigma_low_species',
                                          'sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)

#Extracting some infos
draws       <- rstan::extract(fit_allPhi4)
state_draws <- draws$state


# Partial pooling effect --------------------------------------------------
# Check the actual mapping
levels(as.factor(stand_year_all$spp))

# Define sp_names first
sp_names <- c("ABAM", "ABLA", "CANO", "PSME", "THPL", "TSHE")

shrinkage_df <- map_dfr(seq_len(G), function(f) {
  idx     <- start_idxs[f]:end_idxs[f]
  s_id    <- stand_year_all$species_id[idx[1]]
  s_name  <- sp_names[s_id]
  st_name <- years_per_series$stand[f]
  
  tibble(
    series      = f,
    sp_name     = s_name,        
    stand       = st_name,
    obs_mean    = mean(stan_data_all$y[idx]),
    obs_max     = max(stan_data_all$y[idx]),
    p_mast_mean = mean(colMeans(state_draws[, idx, drop = FALSE] == 2)),
    n_years     = length(idx)
  )
})

ggplot(shrinkage_df, aes(obs_mean, p_mast_mean, label = stand, colour = sp_name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_colour_brewer(palette = "Dark2") +
  labs(title    = "Partial pooling influence: observed mean vs P(mast)",
       subtitle = "High P(mast) with low counts = over-pooled by partial pooling",
       x = "Observed mean seed count",
       y = "Mean P(mast state)",
       colour = "Species") +
  theme_minimal()

# PPC ---------------------------------------------------------------------

#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_all$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_all$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G) {
  idx    <- start_idxs[f]:end_idxs[f]
  years_f <- stand_year_all$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_all$y[idx], pch = 16, cex = 0.8)
}

#Problematic things
par(mfrow = c(1, 3))
for (f in c(4, 26, 32)) {
  idx     <- start_idxs[f]:end_idxs[f]
  years_f <- stand_year_all$year[idx]
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples, pred_names_f, years_f,
    xlab = paste("Series", f), ylab = "y"
  )
  points(years_f, stan_data_all$y[idx], pch = 16, cex = 0.8)
}

f <- 26
idx     <- start_idxs[f]:end_idxs[f]
years_f <- stand_year_all$year[idx]
y_obs_f <- stan_data_all$y[idx]
pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))

util$plot_conn_pushforward_quantiles(
  samples, pred_names_f, years_f,
  xlab = paste("Series", f), ylab = "y"
)
# Add points using original coordinates, clip manually
points(years_f, pmin(y_obs_f, 400), pch = 16, cex = 0.8)



# State identification ----------------------------------------------------


# Prior vs posterior ------------------------------------------------------


