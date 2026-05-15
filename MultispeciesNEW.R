##Code for multispecies HMM model 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

##Next things to do : 

#1 : Add the 2009 year, handle the SUNR missing year 
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
  filter(spp %in% c("ABAM","ABLA","CANO","PSME","TSHE","TSME","THPL"))
unique(seed_data$stand)
#I include ABLA now

# Extra stand
# Add SPRY, SUNR, PARA for CANO, ABAM, ABLA
extra_stands <- seed_data %>%
  filter(spp %in% c("CANO", "ABAM", "ABLA"),
         stand %in% c("SPRY", "SUNR", "PARA")) %>%
  filter(year != 2009) %>%
  group_by(spp, stand, year) %>%
  summarise(
    y    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(spp, stand, year)
# Checking
extra_stands %>%
  group_by(spp, stand) %>%
  summarise(n_years = n(), mean_density = mean(y/area), .groups = "drop")


# 2) Removing invalid stands 
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)



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

stand_year_all <- bind_rows(stand_year_all, extra_stands) %>%
  arrange(spp, stand, year)

stand_year_all <- stand_year_all %>%
  filter(!(stand == "SUNR" & year == 2024))



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

print (stan_data_all)

# Fitting Model -----------------------------------------------------------

fit_all <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesNEWMike.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_all2 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesNEWNCThetaAllHigh.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

# Diagnostic plots Centered --------------------------------------------------------


#If there is a divergence : 

##Add the divergence color 

pairs(
  fit_all,
  pars = c(
    "log_phi_state[1]",
    "log_phi_state[2]"
  ),
  condition = "divergent__"
)

pairs(
  fit_all,
  pars = c(
    "mu_log_low",
    "mu_log_delta",
    "sigma_low_species",
    "sigma_low_stand",
    "sigma_log_delta_species",
    "sigma_log_delta_stand"
  ),
  condition = "divergent__"
)

pairs(
  fit_all,
  pars = c(
    "alpha_low_species [1]",
    "alpha_low_stand [1]",
    "log_delta_species [1]",
    "log_delta_stand [1]",
    "sigma_low_species",
    "sigma_low_stand",
    "sigma_log_delta_species",
    "sigma_log_delta_stand"
  ),
  condition = "divergent__"
)

pairs(
  fit_all,
  pars = c(
    "alpha_theta1_species[1]",
    "sigma_theta1_species",
    "alpha_theta2_species[1]",
    "sigma_theta2_species",
    "alpha_theta1_stand[1]",
    "sigma_theta1_stand",
    "alpha_theta2_stand[1]",
    "sigma_theta2_stand"
  ),
  condition = "divergent__"
)

# Diagnostic plots Non-centered  --------------------------------------------------------

#Non-centering on the thetas and the high seed production

pairs(
  fit_all,
  pars = c(
    "log_phi_state[1]",
    "log_phi_state[2]"
  ),
  condition = "divergent__"
)


pairs(
  fit_all,
  pars = c(
    "alpha_low_species [1]",
    "alpha_low_stand [1]",
    "log_delta_species_nc [1]",
    "log_delta_stand_nc [1]",
    "sigma_low_species",
    "sigma_low_stand",
    "sigma_log_delta_species",
    "sigma_log_delta_stand"
  ),
  condition = "divergent__"
)

pairs(
  fit_all,
  pars = c(
    "alpha_theta1_species_nc[1]",
    "sigma_theta1_species",
    "alpha_theta2_species_nc[1]",
    "sigma_theta2_species",
    "alpha_theta1_stand_nc[1]",
    "sigma_theta1_stand",
    "alpha_theta2_stand_nc[1]",
    "sigma_theta2_stand"
  ),
  condition = "divergent__"
)




# Results -----------------------------------------------------------------

summary(fit_all2)

s <- summary(fit_all2)$summary

species_list <- c("ABAM", "ABLA", "CANO", "PSME", "THPL", "TSHE", "TSME")
results <- list()

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  
  low_mean    <- exp(s["mu_log_low", "50%"] + s[paste0("alpha_low_species[", i, "]"), "50%"])
  delta_sp    <- s["mu_log_delta", "50%"] + s["sigma_log_delta_species", "50%"] * s[paste0("log_delta_species_nc[", i, "]"), "50%"]
  fold_change <- 1 + exp(delta_sp)
  high_mean   <- low_mean * fold_change
  
  p_stay_low  <- plogis(s["grand_logit_theta1", "50%"] + s["sigma_theta1_species", "50%"] * s[paste0("alpha_theta1_species_nc[", i, "]"), "50%"])
  p_stay_high <- plogis(s["grand_logit_theta2", "50%"] + s["sigma_theta2_species", "50%"] * s[paste0("alpha_theta2_species_nc[", i, "]"), "50%"])
  
  results[[sp]] <- data.frame(
    species                = sp,
    low_state_mean         = round(low_mean, 3),
    high_state_fold_change = round(fold_change, 3),
    high_state_mean        = round(high_mean, 3),
    p_stay_low             = round(p_stay_low, 3),
    p_stay_high            = round(p_stay_high, 3),
    sigma_low_stand        = round(s["sigma_low_stand",         "50%"], 3),
    sigma_log_delta_sp     = round(s["sigma_log_delta_species", "50%"], 3),
    sigma_log_delta_stand  = round(s["sigma_log_delta_stand",   "50%"], 3),
    phi_low                = round(exp(s["log_phi_state[1]",    "50%"]), 3),
    phi_high               = round(exp(s["log_phi_state[2]",    "50%"]), 3)
  )
}

summary_table <- do.call(rbind, results)
rownames(summary_table) <- NULL
print(summary_table, digits = 3)


# Mike's plot -------------------------------------------------------------



util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_all2)
util$check_all_hmc_diagnostics(diagnostics)
#there is 1% of divergence

samples <- util$extract_expectand_vals(fit_all2)

base_samples <- util$filter_expectands(samples,
                                       c(# Transitions - species NC (what's in the model)
                                         paste0('alpha_theta1_species_nc[', 1:7, ']'),
                                         paste0('alpha_theta2_species_nc[', 1:7, ']'),
                                         paste0('alpha_theta1_stand_nc[',   1:18, ']'),
                                         paste0('alpha_theta2_stand_nc[',   1:18, ']'),
                                         'sigma_theta1_species', 'sigma_theta2_species',
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Emission means
                                          'mu_log_low',
                                          paste0('alpha_low_species[', 1:7, ']'),
                                          'sigma_low_species',
                                          paste0('alpha_low_stand[',   1:18, ']'),
                                          'sigma_low_stand',
                                          
                                          # Delta parameters
                                          'mu_log_delta',
                                          'sigma_log_delta_species',
                                          'sigma_log_delta_stand',
                                          paste0('log_delta_species_nc[', 1:7, ']'),
                                          paste0('log_delta_stand_nc[',   1:18, ']'),
                                          
                                          # Overdispersion (per state, not per species)
                                          'log_phi_state[1]',
                                          'log_phi_state[2]'
                                        ))
util$check_all_expectand_diagnostics(base_samples2)



# PPC ---------------------------------------------------------------------

#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")

samples <- util$extract_expectand_vals(fit_all2)

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
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
samples <- util$extract_expectand_vals(fit_all2)
par(mfrow = c(4,4))

for (f in 1:length(stan_data_all$start_idxs)) {
  
  start_id <- start_idxs[f]
  end_id   <- end_idxs[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_all$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(stan_data_all$start_idxs), " -- Press [Enter] for next..."))
}

# Prior vs posterior ------------------------------------------------------



samples <- rstan::extract(fit_all2)

# --- Simulate priors (same sample size as posterior) ---
n <- length(samples$mu_log_low)

priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_species    = abs(rnorm(n, 0, 0.5)),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta_sp   = abs(rnorm(n, 0, 1.0)),
  sigma_log_delta_stand= abs(rnorm(n, 0, 1.0)),
  sigma_theta1_species = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_species = abs(rnorm(n, 0, 0.7)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
  phi_low              = exp(rnorm(n, log(4), 0.6)),
  phi_high             = exp(rnorm(n, log(4), 0.6))
)

posteriors <- data.frame(
  mu_log_low           = samples$mu_log_low,
  mu_log_delta         = samples$mu_log_delta,
  grand_logit_theta1   = samples$grand_logit_theta1,
  grand_logit_theta2   = samples$grand_logit_theta2,
  sigma_low_species    = samples$sigma_low_species,
  sigma_low_stand      = samples$sigma_low_stand,
  sigma_log_delta_sp   = samples$sigma_log_delta_species,
  sigma_log_delta_stand= samples$sigma_log_delta_stand,
  sigma_theta1_species = samples$sigma_theta1_species,
  sigma_theta2_species = samples$sigma_theta2_species,
  sigma_theta1_stand   = samples$sigma_theta1_stand,
  sigma_theta2_stand   = samples$sigma_theta2_stand,
  phi_low              = exp(samples$log_phi_state[, 1]),
  phi_high             = exp(samples$log_phi_state[, 2])
)

# --- Nice labels ---
param_labels <- c(
  mu_log_low            = "µ log low (grand mean seeds, low state)",
  mu_log_delta          = "µ log delta (grand mean fold-change)",
  grand_logit_theta1    = "logit θ₁ (P stay low)",
  grand_logit_theta2    = "logit θ₂ (P stay high)",
  sigma_low_species     = "σ low species",
  sigma_low_stand       = "σ low stand",
  sigma_log_delta_sp    = "σ log delta species",
  sigma_log_delta_stand = "σ log delta stand",
  sigma_theta1_species  = "σ θ₁ species",
  sigma_theta2_species  = "σ θ₂ species",
  sigma_theta1_stand    = "σ θ₁ stand",
  sigma_theta2_stand    = "σ θ₂ stand",
  phi_low               = "φ low state (overdispersion)",
  phi_high              = "φ high state (overdispersion)"
)

# --- Reshape to long ---
prior_long <- priors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(type = "Prior")

post_long <- posteriors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(type = "Posterior")

df <- bind_rows(prior_long, post_long) %>%
  mutate(
    label = param_labels[param],
    type  = factor(type, levels = c("Prior", "Posterior"))
  )

# --- Plot ---
ggplot(df, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_fill_manual(values  = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
  scale_color_manual(values = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
  labs(
    title = "Prior vs posterior distributions",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
