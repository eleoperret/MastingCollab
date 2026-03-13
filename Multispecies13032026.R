##Code for multispecies HMM model 
# Stan code: 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

####Doesn't include PARA and SUNR as I do not have a map of the species. 
####Issues with ABLA so I will remove it. 

#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)

options(mc.cores = parallel::detectCores())


#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

seed_data<-read.csv("SeedData_all.csv")

# 1) Keep species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","CANO","PSME","TSHE","THPL"))
#I did not include Abla because it is present only at two stands and the other ones are the PARA and SUNR which I still do not have the data for


# 2) Removing invalid stands 
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)

total_stand<-seed_filtered%>%
  group_by(spp,year)%>%
  summarise(total_seeds=sum(total_viable_sds),na.rm=TRUE)

# -------------------------------------------------------
# 1. Data preparation — all species including PSME
# -------------------------------------------------------

stand_year_all <- seed_filtered %>%
  filter(year != 2009) %>%
  group_by(spp, stand, year) %>%
  summarise(
    y    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(spp, stand, year)

#Checking the distribution of my species
stand_year_all %>%
  ggplot(aes(x = log1p(y))) +
  geom_histogram(bins = 30) +
  facet_wrap(~spp, scales = "free") +
  labs(x = "log(seeds + 1)", title = "Are distributions bimodal?")
#I have some species that are really heavy seed producer compared to others. 

#Checking more things to define my priors
stand_year_all %>%
  group_by(spp) %>%
  summarise(
    n_stands    = n_distinct(stand),
    year_min    = min(year),
    year_max    = max(year),
    n_obs       = n(),
    n_zeros     = sum(y == 0),
    pct_zeros   = round(100 * mean(y == 0), 1),
    mean_y_nonzero = round(mean(y[y > 0]), 1)
  )

# Look at quantiles per species to understand state separation better
stand_year_all %>%
  filter(y > 0) %>%
  group_by(spp) %>%
  summarise(
    q10  = quantile(y, 0.10),
    q25  = quantile(y, 0.25),
    q50  = quantile(y, 0.50),
    q75  = quantile(y, 0.75),
    q90  = quantile(y, 0.90),
    max  = max(y)
  )

# Species index
stand_year_all$species_id <- as.numeric(as.factor(stand_year_all$spp))
S <- length(unique(stand_year_all$spp))

# Checking my species order before running
cat("Species order — check this matches your Stan priors:\n")
print(levels(as.factor(stand_year_all$spp)))

# Creating my stan data list
years_per_series <- stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(T_i = n(), .groups = "drop")

G          <- nrow(years_per_series)
T_i        <- years_per_series$T_i
start_idxs <- cumsum(c(1, T_i[-G]))
end_idxs   <- cumsum(T_i)

# Checking everything
cat("N =", nrow(stand_year_all), "\n")
cat("F =", G, "\n")
cat("S =", S, "\n")
cat("Indices correct:", all(end_idxs - start_idxs + 1 == T_i), "\n")
cat("Final index matches N:", tail(end_idxs, 1) == nrow(stand_year_all), "\n")

#Final Stan data list
stan_data_all <- list(
  N          = nrow(stand_year_all),
  F          = G,
  S          = S,
  sp         = stand_year_all$species_id,
  start_idxs = start_idxs,
  end_idxs   = end_idxs,
  y          = stand_year_all$y,
  area       = stand_year_all$area
)

# -------------------------------------------------------
# 2. Initial values
# -------------------------------------------------------
# Adjust low/high state starting values per species
# based on your histogram inspection:
# ABAM(1): low~0.5,  high~3.0
# CANO(2): low~0.5,  high~3.5
# PSME(3): low~3.0,  high~4.5  <- higher both states, close together
# THPL(4): low~2.0,  high~5.0
# TSHE(5): low~4.5,  high~5.5

init_fn_all <- function() {       # returns a LIST — Stan expects this format
  list(
    rho = c(0.8, 0.2),            # transition probs: start chains assuming
    # mostly low state (70%) vs high (30%)
    
    theta1 = rep(0.80, stan_data_all$S),        # initial state probs for state 1, one per species
    theta2 = rep(0.20, stan_data_all$S),        # initial state probs for state 2, one per species
    
    log_means = matrix(           # S×2 matrix of starting log-scale means
      c(0.5, 1.5, 0.5, 3.0, 4.0, #   col 1: low  state (ABAM, CANO, PSME, THPL, TSHE)
        5.2, 5.4, 5.2, 6.6, 7.4),#   col 2: high state
      nrow = stan_data_all$S, ncol = 2
    ),                            # note: these match your q90 comments exactly
    
    log_phi1 = rep(log(2), stan_data_all$S),  # log-scale dispersion for state 1, all species
    log_phi2 = rep(log(6.5), stan_data_all$S),  # log-scale dispersion for state 2, all species
    
    sigma = 0.5,                  # starting value for some SD parameter
    
    stand_effect_raw = rep(0, stan_data_all$F)  # random effects for F stands, start at zero
    # (neutral / no stand effect)
  )
}

# -------------------------------------------------------
# 3. Fit
# -------------------------------------------------------

fit_all <- stan(
  file    = "Stan_code/Species_Stan_Model/Multispecies4.stan",
  data    = stan_data_all,
  iter    = 4000,
  warmup  = 2000,
  chains  = 4,
  seed    = 123,
  init    = init_fn_all,
  control = list(
    adapt_delta   = 0.95,
    max_treedepth = 12
  )
)


# 4. Diagnostics --------------------------------------------------------

#New style

print(fit_all, pars = c("theta1", "theta2", "log_means",
                        "log_phi1", "log_phi2", "sigma", "rho"))

rhat_vals <- rhat(fit_all)
cat("\nRhat > 1.05 (excluding state and log_omega):\n")
bad_rhat <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.05]
print(bad_rhat[!grepl("state|log_omega|y_rep", names(bad_rhat))])

cat("\nESS ratio < 0.1 (excluding state and log_omega):\n")
neff_vals <- neff_ratio(fit_all)
bad_neff  <- neff_vals[!is.na(neff_vals) & neff_vals < 0.1]
print(bad_neff[!grepl("state|log_omega|y_rep", names(bad_neff))])

# Trace plots for emission means — most important to check
posterior_all <- as.array(fit_all)
mcmc_trace(posterior_all,
           pars = c("log_means[1,1]", "log_means[1,2]",
                    "log_means[2,1]", "log_means[2,2]",
                    "log_means[3,1]", "log_means[3,2]"))
mcmc_trace(posterior_all,
           pars = c("log_means[4,1]", "log_means[4,2]",
                    "log_means[5,1]", "log_means[5,2]",
                    "sigma"))



#Old style Based on Mike's code 

util<- new.env()
#creating a source file to make all the plots. 
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit_all)
diagnostics <- util$extract_hmc_diagnostics(fit_all)
util$check_all_hmc_diagnostics(diagnostics)

# extraire les posterior values
samples <- util$extract_expectand_vals(fit_all)

# diagnostics parametre par parametre
base_samples <- util$filter_expectands(samples,
                                       c("rho", "theta1","theta2",
                                         "log_lambda", "log_mu",
                                         "stand_effect_raw", "phi1", "phi2","sigma"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

#Ordering the stand based on elevation for PSME
original_order <- unique(stand_year_df$stand)

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

par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:stan_data_all$F){
  idxs <- stan_data_all$start_idxs[s]:stan_data_all$end_idxs[s]
  subsamples <-  util$filter_expectands(samples,
                                        paste0("y_rep[", idxs, "]"))
  util$plot_hist_quantiles(subsamples, "y_rep",
                           bin_min = 0, bin_max = 500, bin_delta = 10,
                           baseline_values = stan_data_all$y[idxs])
}


#For plotting
#Adding the stand elevation and names
perm <- match(stand_order, original_order)

elevation<-readRDS(stand_elevation_table.rds)

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




# New plotting ------------------------------------------------------------

# -------------------------------------------------------
# 0. Extract everything once
# -------------------------------------------------------

posterior_draws  <- rstan::extract(fit_all)
y_rep            <- posterior_draws$y_rep          # [draws × N]
state_draws      <- posterior_draws$state          # [draws × N]
log_means_draws  <- posterior_draws$log_means      # [draws × S × 2]
theta1_draws     <- posterior_draws$theta1         # [draws × S]
theta2_draws     <- posterior_draws$theta2         # [draws × S]
sigma_draws      <- posterior_draws$sigma          # [draws]
rho_draws        <- posterior_draws$rho            # [draws × 2]
log_phi1_draws   <- posterior_draws$log_phi1       # [draws × S]
log_phi2_draws   <- posterior_draws$log_phi2       # [draws × S]
log_alpha_draws  <- posterior_draws$log_alpha      # [draws × F]

species_names <- levels(as.factor(stand_year_all$spp))
# Should be: ABAM CANO PSME THPL TSHE

# -------------------------------------------------------
# 1. Posterior predictive checks
# -------------------------------------------------------

# Overall fit
ppc_dens_overlay(
  y     = stan_data_all$y,
  yrep  = y_rep[1:100, ],
  trans = "log1p"
) +
  ggtitle("Posterior predictive check — all species")

# Per species
ppc_dens_overlay_grouped(
  y     = stan_data_all$y,
  yrep  = y_rep[1:100, ],
  group = stand_year_all$spp,
  trans = "log1p"
) +
  ggtitle("Posterior predictive check — per species")

# Zero rate per species
ppc_stat_grouped(
  y     = stan_data_all$y,
  yrep  = y_rep,
  group = stand_year_all$spp,
  stat  = function(y) mean(y == 0)
) +
  ggtitle("Proportion of zeros — observed vs predicted")

# Mean per species
ppc_stat_grouped(
  y     = stan_data_all$y,
  yrep  = y_rep,
  group = stand_year_all$spp,
  stat  = "mean"
) +
  ggtitle("Mean seed count — observed vs predicted")

# -------------------------------------------------------
# 2. State assignment summary
# -------------------------------------------------------

# Posterior mode state per observation
state_mode <- apply(state_draws, 2,
                    function(x) as.integer(names(which.max(table(x)))))

stand_year_all$state <- state_mode

# Summary per species
state_summary <- stand_year_all %>%
  group_by(spp) %>%
  summarise(
    n_obs       = n(),
    n_low       = sum(state == 1),
    n_high      = sum(state == 2),
    pct_masting = round(100 * mean(state == 2), 1),
    .groups     = "drop"
  )
print(state_summary)

# Per stand
masting_by_stand <- stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(
    n_years     = n(),
    n_masting   = sum(state == 2),
    pct_masting = round(100 * mean(state == 2), 1),
    .groups     = "drop"
  )
print(masting_by_stand)

# -------------------------------------------------------
# 3. Masting year timeline
# -------------------------------------------------------

stand_year_all %>%
  group_by(spp, year) %>%
  summarise(
    mean_seeds  = mean(y),
    pct_masting = round(100 * mean(state == 2), 1),
    .groups     = "drop"
  ) %>%
  ggplot(aes(x = year, y = pct_masting, colour = spp)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~spp, ncol = 1) +
  labs(
    title = "Percentage of stands in masting state per year",
    y     = "% stands in high state",
    x     = "Year"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# -------------------------------------------------------
# 4. Posterior distributions — transitions
# -------------------------------------------------------

posterior_array <- as.array(fit_all)

mcmc_areas(
  posterior_array,
  pars = c("theta1[1]", "theta1[2]", "theta1[3]",
           "theta1[4]", "theta1[5]"),
  prob = 0.90
) +
  ggtitle("theta1: prob of staying in low state per species")

mcmc_areas(
  posterior_array,
  pars = c("theta2[1]", "theta2[2]", "theta2[3]",
           "theta2[4]", "theta2[5]"),
  prob = 0.90
) +
  ggtitle("theta2: prob of staying in high state per species")

# -------------------------------------------------------
# 5. Emission means on natural scale
# -------------------------------------------------------

natural_scale <- data.frame(
  species    = rep(species_names, each = 2),
  state      = rep(c("low", "high"), times = length(species_names)),
  mean_seeds = NA,
  lo90       = NA,
  hi90       = NA
)

for (s in seq_along(species_names)) {
  for (k in 1:2) {
    vals <- exp(log_means_draws[, s, k])
    natural_scale$mean_seeds[(s - 1) * 2 + k] <- median(vals)
    natural_scale$lo90[(s - 1) * 2 + k]       <- quantile(vals, 0.05)
    natural_scale$hi90[(s - 1) * 2 + k]       <- quantile(vals, 0.95)
  }
}

print(natural_scale)

ggplot(natural_scale,
       aes(x = species, y = mean_seeds,
           colour = state, shape = state)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(
    aes(ymin = lo90, ymax = hi90),
    width    = 0.15,
    position = position_dodge(width = 0.3)
  ) +
  scale_y_log10() +
  labs(
    title = "Posterior emission means — low vs high state",
    y     = "Expected seed count (log scale)",
    x     = "Species"
  ) +
  theme_bw()

# -------------------------------------------------------
# 6. Pairs plot for PSME — diagnose bimodal issue
# -------------------------------------------------------

pairs(fit_all,
      pars = c("log_means[3,1]", "log_means[3,2]",
               "log_phi1[3]",    "log_phi2[3]",
               "theta1[3]",      "theta2[3]"))

# -------------------------------------------------------
# 7. Stan's util plots if you have Mike's tools loaded
# -------------------------------------------------------

if (exists("util")) {
  
  samples <- util$extract_expectand_vals(fit_all)
  
  # HMC diagnostics
  diagnostics <- util$extract_hmc_diagnostics(fit_all)
  util$check_all_hmc_diagnostics(diagnostics)
  
  # Parameter diagnostics — using correct parameter names
  base_samples <- util$filter_expectands(
    samples,
    c("rho", "theta1", "theta2",
      "log_means", "log_phi1", "log_phi2",
      "sigma", "stand_effect_raw"),
    check_arrays = TRUE
  )
  util$check_all_expectand_diagnostics(base_samples)
  
  # PPC histogram — overall
  par(mfrow = c(1, 1))
  util$plot_hist_quantiles(
    samples, "y_rep",
    bin_min        = 0,
    bin_max        = 600,
    bin_delta      = 20,
    baseline_values = stan_data_all$y
  )
  
  # PPC per stand
  par(mfrow = c(4, 3), mar = c(5, 5, 2, 1))
  for (f in 1:stan_data_all$F) {
    idxs <- stan_data_all$start_idxs[f]:stan_data_all$end_idxs[f]
    subsamples <- util$filter_expectands(
      samples,
      paste0("y_rep[", idxs, "]")
    )
    spp_name   <- stand_year_all$spp[idxs[1]]
    stand_name <- stand_year_all$stand[idxs[1]]
    util$plot_hist_quantiles(
      subsamples, "y_rep",
      bin_min         = 0,
      bin_max         = 500,
      bin_delta       = 10,
      baseline_values = stan_data_all$y[idxs],
      main            = paste0(stand_name, " (", spp_name, ")")
    )
  }
  
  # rho
  par(mfrow = c(1, 2), mar = c(5, 5, 1, 1))
  util$plot_expectand_pushforward(
    samples[["rho[1]"]], 50, flim = c(0, 1),
    display_name = bquote(rho * " (Low)")
  )
  util$plot_expectand_pushforward(
    samples[["rho[2]"]], 50, flim = c(0, 1),
    display_name = bquote(rho * " (High)")
  )
  
  # theta per species
  par(mfrow = c(2, 5), mar = c(5, 5, 2, 1))
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      samples[[paste0("theta1[", s, "]")]], 50, flim = c(0, 1),
      display_name = paste0("theta1 ", species_names[s])
    )
  }
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      samples[[paste0("theta2[", s, "]")]], 50, flim = c(0, 1),
      display_name = paste0("theta2 ", species_names[s])
    )
  }
  
  # log_means per species on natural scale
  par(mfrow = c(2, 5), mar = c(5, 5, 2, 1))
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      exp(samples[[paste0("log_means[", s, ",1]")]]), 50,
      display_name = paste0("lambda ", species_names[s])
    )
  }
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      exp(samples[[paste0("log_means[", s, ",2]")]]), 50,
      display_name = paste0("mu ", species_names[s])
    )
  }
  
  # sigma
  par(mfrow = c(1, 1))
  util$plot_expectand_pushforward(
    samples[["sigma"]], 50,
    display_name = expression(sigma)
  )
  
  # log_alpha per stand
  par(mfrow = c(4, 3), mar = c(5, 5, 2, 1))
  for (f in 1:stan_data_all$F) {
    spp_name   <- stand_year_all$spp[stan_data_all$start_idxs[f]]
    stand_name <- stand_year_all$stand[stan_data_all$start_idxs[f]]
    util$plot_expectand_pushforward(
      exp(samples[[paste0("log_alpha[", f, "]")]]), 50,
      display_name = paste0(stand_name, " (", spp_name, ")")
    )
  }
}

