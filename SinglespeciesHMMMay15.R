##Code for single species HMM model (for all species)
# Stan code: 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

##Next things to do : 

#1 : Add the 2009 year, handle the SUNR missing year
#2 : Trap-level
#3 : Extract values for the synchrony.
#4 : Clean my code (structure of the code and clean the top part)

###DATA:
#6 species 
#Starting in 2010 ( as not all the stands start in the same year)
#With SPRY, SUNR and PARA manually added for ABAM; ABLA and CANO
#For SUNR data stops at 2022 (as 2023 is a missing year and do not know how to deal with this yet)
#At the stand NOT trap level
#Combination of stand x specie (aggregated) seed production per year

#Changed the priors based on : SinglespeciesPartialPooling2_delta_NCbisFULL

#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)
library(tidyverse)
library(posterior)
library(boot)

options(mc.cores = parallel::detectCores())

#To run to use Mike's package
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")


# Data import and cleaning-------------------------------------------------------------


seed_data<-read.csv("SeedData_all.csv")

# 1) Keep species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","ABLA","CANO","PSME","TSHE","THPL", "TSME"))
unique(seed_data$stand)
#I include ABLA now

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
  filter(year != 2009) %>% 
  group_by(spp, stand, year) %>%
  summarise(
    y    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(spp, stand, year)

str(stand_year_all)

stand_year_all <- bind_rows(stand_year_all, extra_stands) %>%
  arrange(spp, stand, year)

# stand_year_all <- stand_year_all %>%
#   filter(!(stand == "SUNR" & year == 2024))

missing_rows <- tibble(
  spp   = c("ABAM", "ABLA", "CANO"),
  stand = "SUNR",
  year  = 2023,
  y     = NA_integer_,
  area  = NA_real_
)

stand_year_all <- bind_rows(stand_year_all, missing_rows) %>%
  group_by(spp, stand) %>%
  mutate(area = ifelse(is.na(area), mean(area, na.rm = TRUE), area)) %>%
  ungroup() %>%
  arrange(spp, stand, year)

obs_missing <- as.integer(!is.na(stand_year_all$y))

stand_year_all <- stand_year_all %>%
  mutate(y = ifelse(is.na(y), 0L, as.integer(y)))

# ABAM --------------------------------------------------------------------
#Per species
stand_year_abam <- stand_year_all %>%
  filter(spp== "ABAM")

# Creating my stan data list
years_per_series_abam <- stand_year_abam %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_abam <- length(unique(years_per_series_abam$stand_id))  # NEW

G_abam          <- nrow(years_per_series_abam) 
T_i_abam        <- years_per_series_abam$T_i 
start_idxs_abam <- cumsum(c(1, T_i_abam[-G_abam])) 
end_idxs_abam    <- cumsum(T_i_abam)

stan_data_abam <- list(
  N          = nrow(stand_year_abam),
  F          = G_abam,
  N_stands   = N_stands_abam, 
  stand_id   = years_per_series_abam$stand_id,
  start_idxs = start_idxs_abam,
  end_idxs   = end_idxs_abam,
  y          = stand_year_abam$y,
  area       = stand_year_abam$area, 
  obs_missing = obs_missing[stand_year_all$spp=="ABAM"]
)

#print(stan_data_abam)

# ABLA  --------------------------------------------------------------------
#Per species
stand_year_abla <- stand_year_all %>%
  filter(spp== "ABLA")

# Creating my stan data list
years_per_series_abla <- stand_year_abla %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_abla <- length(unique(years_per_series_abla$stand_id))  # NEW

G_abla          <- nrow(years_per_series_abla) 
T_i_abla        <- years_per_series_abla$T_i 
start_idxs_abla <- cumsum(c(1, T_i_abla[-G_abla])) 
end_idxs_abla   <- cumsum(T_i_abla)

stan_data_abla <- list(
  N          = nrow(stand_year_abla),
  F          = G_abla,
  N_stands   = N_stands_abla, 
  stand_id   = years_per_series_abla$stand_id,
  start_idxs = start_idxs_abla,
  end_idxs   = end_idxs_abla,
  y          = stand_year_abla$y,
  area       = stand_year_abla$area,
  obs_missing = obs_missing[stand_year_all$spp=="ABLA"]
)

print(stan_data_abla)

# CANO  --------------------------------------------------------------------
#Per species
stand_year_cano <- stand_year_all %>%
  filter(spp== "CANO")

# Creating my stan data list
years_per_series_cano <- stand_year_cano %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_cano <- length(unique(years_per_series_cano$stand_id))  # NEW

G_cano          <- nrow(years_per_series_cano) 
T_i_cano       <- years_per_series_cano$T_i 
start_idxs_cano <- cumsum(c(1, T_i_cano[-G_cano])) 
end_idxs_cano   <- cumsum(T_i_cano)

stan_data_cano <- list(
  N          = nrow(stand_year_cano),
  F          = G_cano,
  N_stands   = N_stands_cano, 
  stand_id   = years_per_series_cano$stand_id,
  start_idxs = start_idxs_cano,
  end_idxs   = end_idxs_cano,
  y          = stand_year_cano$y,
  area       = stand_year_cano$area,
  obs_missing = obs_missing[stand_year_all$spp=="CANO"]
)

print(stan_data_cano)

#works
# PSME  --------------------------------------------------------------------
#Per species
stand_year_psme <- stand_year_all %>%
  filter(spp== "PSME")

# Creating my stan data list
years_per_series_psme <- stand_year_psme %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_psme <- length(unique(years_per_series_psme$stand_id))  # NEW

G_psme          <- nrow(years_per_series_psme) 
T_i_psme      <- years_per_series_psme$T_i 
start_idxs_psme <- cumsum(c(1, T_i_psme[-G_psme])) 
end_idxs_psme   <- cumsum(T_i_psme)

stan_data_psme <- list(
  N          = nrow(stand_year_psme),
  F          = G_psme,
  N_stands   = N_stands_psme, 
  stand_id   = years_per_series_psme$stand_id,
  start_idxs = start_idxs_psme,
  end_idxs   = end_idxs_psme,
  y          = stand_year_psme$y,
  area       = stand_year_psme$area,
  obs_missing = rep(1L, nrow(stand_year_psme)) 
)

print(stan_data_psme)


# TSHE  --------------------------------------------------------------------
#Per species
stand_year_tshe <- stand_year_all %>%
  filter(spp== "TSHE")

# Creating my stan data list
years_per_series_tshe <- stand_year_tshe %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_tshe <- length(unique(years_per_series_tshe$stand_id))  # NEW

G_tshe         <- nrow(years_per_series_tshe) 
T_i_tshe     <- years_per_series_tshe$T_i 
start_idxs_tshe <- cumsum(c(1, T_i_tshe[-G_tshe])) 
end_idxs_tshe   <- cumsum(T_i_tshe)

stan_data_tshe <- list(
  N          = nrow(stand_year_tshe),
  F          = G_tshe,
  N_stands   = N_stands_tshe, 
  stand_id   = years_per_series_tshe$stand_id,
  start_idxs = start_idxs_tshe,
  end_idxs   = end_idxs_tshe,
  y          = stand_year_tshe$y,
  area       = stand_year_tshe$area,
  obs_missing = rep(1L, nrow(stand_year_tshe))
)

print(stan_data_tshe)




# TSME  --------------------------------------------------------------------
#Per species
stand_year_tsme <- stand_year_all %>%
  filter(spp== "TSME")

# Creating my stan data list
years_per_series_tsme <- stand_year_tsme %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_tsme <- length(unique(years_per_series_tsme$stand_id))  # NEW

G_tsme         <- nrow(years_per_series_tsme) 
T_i_tsme     <- years_per_series_tsme$T_i 
start_idxs_tsme <- cumsum(c(1, T_i_tsme[-G_tsme])) 
end_idxs_tsme   <- cumsum(T_i_tsme)

stan_data_tsme <- list(
  N          = nrow(stand_year_tsme),
  F          = G_tsme,
  N_stands   = N_stands_tsme, 
  stand_id   = years_per_series_tsme$stand_id,
  start_idxs = start_idxs_tsme,
  end_idxs   = end_idxs_tsme,
  y          = stand_year_tsme$y,
  area       = stand_year_tsme$area,
  obs_missing = rep(1L, nrow(stand_year_tsme))
)

print(stan_data_tsme)





# THPL  --------------------------------------------------------------------
#Per species
stand_year_thpl <- stand_year_all %>%
  filter(spp== "THPL")

# Creating my stan data list
years_per_series_thpl <- stand_year_thpl %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_thpl <- length(unique(years_per_series_thpl$stand_id))  # NEW

G_thpl         <- nrow(years_per_series_thpl) 
T_i_thpl    <- years_per_series_thpl$T_i 
start_idxs_thpl <- cumsum(c(1, T_i_thpl[-G_thpl])) 
end_idxs_thpl   <- cumsum(T_i_thpl)

stan_data_thpl <- list(
  N          = nrow(stand_year_thpl),
  F          = G_thpl,
  N_stands   = N_stands_thpl, 
  stand_id   = years_per_series_thpl$stand_id,
  start_idxs = start_idxs_thpl,
  end_idxs   = end_idxs_thpl,
  y          = stand_year_thpl$y,
  area       = stand_year_thpl$area,
  obs_missing = rep(1L, nrow(stand_year_thpl))
)

print(stan_data_thpl)




# Plots of raw seed production --------------------------------------------
ggplot(stand_year_abam, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_abam$year), 
                                  max(stand_year_abam$year), by = 2)) +
  labs(
    title = "ABam seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )

ggplot(stand_year_abla, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_abla$year), 
                                  max(stand_year_abla$year), by = 2)) +
  labs(
    title = "ABLA seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )

ggplot(stand_year_cano, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_cano$year), 
                                  max(stand_year_cano$year), by = 2)) +
  labs(
    title = "CANO seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )

ggplot(stand_year_psme, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_psme$year), 
                                  max(stand_year_psme$year), by = 2)) +
  labs(
    title = "psme seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )


ggplot(stand_year_thpl, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_thpl$year), 
                                  max(stand_year_thpl$year), by = 2)) +
  labs(
    title = "THPL seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )


ggplot(stand_year_tsme, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_tsme$year), 
                                  max(stand_year_tsme$year), by = 2)) +
  labs(
    title = "TSME seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )

ggplot(stand_year_tshe, aes(x = year, y = y)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_point(size = 1.5, color = "black") +
  facet_wrap(~ stand, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(stand_year_tshe$year), 
                                  max(stand_year_tshe$year), by = 2)) +
  labs(
    title = "TSHE seed production per stand and year",
    x     = "Year",
    y     = "Seed count (y)"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank()
  )
# Fitting Model -----------------------------------------------------------

fit_ABAM <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_abam,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_ABLA <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_abla,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_CANO <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_cano,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_PSME <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_psme,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_TSHE <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_tshe,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_TSME <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_tsme,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_THPL <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesFULL_09062026.stan",
  data    = stan_data_thpl,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# Summary per species --------------------------------------------


species_list <- c("ABAM", "ABLA", "CANO", "PSME", "TSHE","TSME","THPL")

fits_list <- list(
  ABAM = fit_ABAM,
  ABLA = fit_ABLA,
  CANO = fit_CANO,
  PSME = fit_PSME,
  TSHE = fit_TSHE,
  TSME = fit_TSME,
  THPL = fit_THPL
)

results <- list()

for (sp in species_list) {
  
  cat("\n=============================\n")
  cat("Species:", sp, "\n")
  cat("=============================\n")
  
  samples <- rstan::extract(fits_list[[sp]])
  
  # new: low state mean is exp(mu_log_low)
  # new: fold change is exp(log1p_exp(mu_log_delta)) = 1 + exp(mu_log_delta)
  # new: high state mean is low * fold_change
  low_mean   <- exp(samples$mu_log_low)
  fold_change <- 1 + exp(samples$mu_log_delta)
  high_mean  <- low_mean * fold_change
  
  cat("=== Emission means ===\n")
  cat("Low state mean (seeds):  ", round(median(low_mean), 2), "\n")
  cat("High state fold-change:  ", round(median(fold_change), 2), "x\n")
  cat("High state mean (seeds): ", round(median(high_mean), 2), "\n")
  
  cat("\n=== Transition probabilities ===\n")
  cat("P(stay low  | low state):  ", round(median(plogis(samples$grand_logit_theta1)), 3), "\n")
  cat("P(stay high | high state): ", round(median(plogis(samples$grand_logit_theta2)), 3), "\n")
  
  cat("\n=== Standard deviations ===\n")
  cat("sigma_low_stand:   ", round(median(samples$sigma_low_stand), 3), "\n")
  cat("sigma_log_delta:   ", round(median(samples$sigma_log_delta), 3), "\n")
  cat("sigma_theta1_stand:", round(median(samples$sigma_theta1_stand), 3), "\n")
  cat("sigma_theta2_stand:", round(median(samples$sigma_theta2_stand), 3), "\n")
  
  cat("\n=== Overdispersion (phi) ===\n")
  cat("phi low state:  ", round(median(exp(samples$log_phi_state[,1])), 3), "\n")
  cat("phi high state: ", round(median(exp(samples$log_phi_state[,2])), 3), "\n")
  
  results[[sp]] <- data.frame(
    species            = sp,
    low_state_mean     = median(low_mean),
    high_state_fold_change = median(fold_change),
    high_state_mean    = median(high_mean),
    p_stay_low         = median(plogis(samples$grand_logit_theta1)),
    p_stay_high        = median(plogis(samples$grand_logit_theta2)),
    sigma_low_stand    = median(samples$sigma_low_stand),
    sigma_log_delta    = median(samples$sigma_log_delta),
    sigma_theta1_stand = median(samples$sigma_theta1_stand),
    sigma_theta2_stand = median(samples$sigma_theta2_stand),
    phi_low            = median(exp(samples$log_phi_state[,1])),
    phi_high           = median(exp(samples$log_phi_state[,2]))
  )
}

summary_table <- do.call(rbind, results)
print(summary_table, digits = 3)


# Diagnostic plots --------------------------------------------------------
# Divergence --------------------------------------------------------------


#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_CANO, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_TSME, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")



# General diagnostics -----------------------------------------------------


#ABAM
diagnostics <- util$extract_hmc_diagnostics(fit_ABAM)
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit_ABAM)
base_samples <- util$filter_expectands(samples,c(paste0('alpha_theta1_stand[',1:N_stands_abam, ']'),paste0('alpha_theta2_stand[',   1:N_stands_abam, ']'),'sigma_theta1_stand','sigma_theta2_stand','grand_logit_theta1','grand_logit_theta2','grand_mean_low','log_delta_high_grand_mean'))
util$check_all_expectand_diagnostics(base_samples)
base_samples2 <- util$filter_expectands(samples,c('log_delta_high_grand_mean','sigma_log_delta_high_stand',paste0('log_phi_state[1]'),paste0('log_phi_state[2]'),'grand_mean_low','sigma_low_stand'))
util$check_all_expectand_diagnostics(base_samples2)

#ABLA
diagnostics <- util$extract_hmc_diagnostics(fit_ABLA)
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit_ABLA)
base_samples <- util$filter_expectands(samples,c(paste0('alpha_theta1_stand[',1:N_stands_abla, ']'),paste0('alpha_theta2_stand[',   1:N_stands_abla, ']'),'sigma_theta1_stand','sigma_theta2_stand','grand_logit_theta1','grand_logit_theta2','grand_mean_low','log_delta_high_grand_mean'))
util$check_all_expectand_diagnostics(base_samples)
base_samples2 <- util$filter_expectands(samples,c('log_delta_high_grand_mean','sigma_log_delta_high_stand',paste0('log_phi_state[1]'),paste0('log_phi_state[2]'),'grand_mean_low','sigma_low_stand'))
util$check_all_expectand_diagnostics(base_samples2)


#Checking the geometry
print(fit_ABAM, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))


print(fit_ABAM, pars = c("grand_logit_theta1",
                          "grand_logit_theta2","mu_log_low",
                          "mu_log_delta",
                          "sigma_log_delta"))
print(fit_ABLA, pars = c("grand_logit_theta1",
                         "grand_logit_theta2","mu_log_low",
                         "mu_log_delta",
                         "sigma_log_delta"))
print(fit_CANO, pars = c("grand_logit_theta1",
                         "grand_logit_theta2","mu_log_low",
                         "mu_log_delta",
                         "sigma_log_delta"))
print(fit_PSME, pars = c("grand_logit_theta1",
                         "grand_logit_theta2","mu_log_low",
                         "mu_log_delta",
                         "sigma_log_delta"))
print(fit_TSHE, pars = c("grand_logit_theta1",
                         "grand_logit_theta2","mu_log_low",
                         "mu_log_delta",
                         "sigma_log_delta"))
print(fit_TSME, pars = c("grand_logit_theta1",
                         "grand_logit_theta2","mu_log_low",
                         "mu_log_delta",
                         "sigma_log_delta"))
print(fit_THPL, pars = c("grand_logit_theta1",
                         "grand_logit_theta2","mu_log_low",
                         "mu_log_delta",
                         "sigma_log_delta"))



# Plots for git issue #67 -------------------------------------------------

ls(util)
util$plot_line_hist
util$plot_line_hists
util$plot_disc_pushforward_quantiles



plot_predicted_vs_observed <- function(spp, bin_max_zoom = 200, bin_delta = 2, y_max = NULL) {
  
  fit      <- fits[[spp]]
  stan_dat <- stan_data_list[[spp]]
  
  state_draws <- as.matrix(fit, pars = "state")
  y_rep_draws <- as.matrix(fit, pars = "y_rep")
  
  y_rep_flat <- as.vector(y_rep_draws)
  state_flat <- as.vector(state_draws)
  
  y_rep_low  <- y_rep_flat[state_flat == 1]
  y_rep_high <- y_rep_flat[state_flat == 2]
  
  util$plot_line_hists(
    values1   = y_rep_low,
    values2   = y_rep_high,
    bin_min   = 0,
    bin_max   = bin_max_zoom,
    bin_delta = bin_delta,
    prob      = TRUE,
    xlab      = "#seeds",
    main      = paste0(spp, " — predicted vs observed"),
    col1      = "lightblue",
    col2      = "darkred"
  )
  
  legend("topright",
         legend = c("Predicted low state", "Predicted high state"),
         col    = c("lightblue", "darkred"),
         lwd    = 2, bty = "n")
}

par(mfrow = c(1, 1))
# Call with fixed y axis
plot_predicted_vs_observed("ABAM", bin_max_zoom = 200, bin_delta = 20)
plot_predicted_vs_observed("ABLA", bin_max_zoom = 200, bin_delta = 20)
plot_predicted_vs_observed("CANO", bin_max_zoom = 200, bin_delta = 20)
plot_predicted_vs_observed("PSME", bin_max_zoom = 200, bin_delta = 20)
plot_predicted_vs_observed("TSHE", bin_max_zoom = 600, bin_delta = 20)
plot_predicted_vs_observed("TSME", bin_max_zoom = 200, bin_delta = 20)
plot_predicted_vs_observed("THPL", bin_max_zoom = 200, bin_delta = 20)



#Seed production over the years
plot_time_predictive <- function(spp) {
  
  fit      <- fits[[spp]]
  stan_dat <- stan_data_list[[spp]]
  meta     <- metas[[spp]]
  
  samples  <- util$extract_expectand_vals(fit)
  n_obs    <- stan_dat$N
  y_names  <- paste0("y_rep[", 1:n_obs, "]")
  
  par(mfrow = c(3, 2), mar = c(5, 4, 4, 1))
  
  for (f in 1:stan_dat$F) {
    idx      <- stan_dat$start_idxs[f]:stan_dat$end_idxs[f]
    y_obs_f  <- stan_dat$y[idx]
    
    util$plot_disc_pushforward_quantiles(
      samples         = samples,
      names           = y_names[idx],
      baseline_values = y_obs_f,
      baseline_col    = "black",
      xlab            = "Year index",
      ylab            = "#seeds",
      main            = paste0(spp, " — ", meta$stand[f])
    )
  }
}

plot_time_predictive("ABAM")
plot_time_predictive("ABLA")
plot_time_predictive("CANO")
plot_time_predictive("PSME")
plot_time_predictive("TSHE")
plot_time_predictive("TSME")
plot_time_predictive("THPL")



#Lines 
plot_time_predictive_lines <- function(spp, n_lines = 50) {
  fit      <- fits[[spp]]
  stan_dat <- stan_data_list[[spp]]
  meta     <- metas[[spp]]
  y_rep_draws <- as.matrix(fit, pars = "y_rep")  # [n_draws x N]
  n_draws     <- nrow(y_rep_draws)
  # Randomly sample a subset of draws to plot
  draw_idx <- sample(1:n_draws, n_lines)
  par(mfrow = c(3, 2), mar = c(5, 4, 4, 1))
  for (f in 1:stan_dat$F) {
    idx     <- stan_dat$start_idxs[f]:stan_dat$end_idxs[f]
    y_obs_f <- stan_dat$y[idx]
    x       <- seq_along(idx)
    # Set up empty plot
    plot(NULL,
         xlim = c(1, length(idx)),
         ylim = c(0, max(c(y_obs_f, y_rep_draws[draw_idx, idx])) * 1.1),
         xlab = "Year index",
         ylab = "#seeds",
         main = paste0(spp, " — ", meta$stand[f]))
    # Draw posterior predictive lines
    for (d in draw_idx) {
      lines(x, y_rep_draws[d, idx], col = adjustcolor("darkred", alpha = 0.1), lwd = 0.8)
    }
    # Observed data on top
    #lines(x, y_obs_f, col = "black", lwd = 2)
    points(x, y_obs_f, pch = 19, col = "black", cex = 0.8)
  }
}
plot_time_predictive_lines("ABAM")


plot_time_predictive_lines <- function(spp) {
  
  fit      <- fits[[spp]]
  stan_dat <- stan_data_list[[spp]]
  meta     <- metas[[spp]]
  
  samples <- util$extract_expectand_vals(fit)
  y_names <- paste0("y_rep[", 1:stan_dat$N, "]")
  
  par(mfrow = c(3, 2), mar = c(5, 4, 4, 1))
  
  for (f in 1:stan_dat$F) {
    idx     <- stan_dat$start_idxs[f]:stan_dat$end_idxs[f]
    y_obs_f <- stan_dat$y[idx]
    plot_xs <- seq_along(idx)  # x axis = year index within stand
    
    util$plot_realizations(
      samples         = samples,
      names           = y_names[idx],
      plot_xs         = plot_xs,
      baseline_values = y_obs_f,
      baseline_col    = "black",
      xlab            = "Year index",
      ylab            = "#seeds",
      main            = paste0(spp, " — ", meta$stand[f])
    )
  }
}

plot_time_predictive_lines("ABAM")


# PPC -------------------------------------------------------------

#ABAM
draws       <- rstan::extract(fit_ABAM)
state_draws <- draws$state
samples <- util$extract_expectand_vals(fit_ABAM)
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_abam$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_abam$y, title("ABAM"))
#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_abam) {
  idx    <- start_idxs_abam[f]:end_idxs_abam[f]
  years_f <- stand_year_abam$year[idx]  
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,          
    xlab = paste("Series", f),
    ylab = "y ABAM"
  )
  points(years_f, stan_data_abam$y[idx], pch = 16, cex = 0.8)
}


#ABLA
samples <- util$extract_expectand_vals(fit_ABLA)
names_yrep <- paste0("y_rep[", 1:stan_data_abla$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_abla$y, title("ABLA"))
#Per series 
par(mfrow = c(2, 3))
for (f in 1:G_abla) {
  idx    <- start_idxs_abla[f]:end_idxs_abla[f]
  years_f <- stand_year_abla$year[idx]  
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           
    xlab = paste("Series", f),
    ylab = "y", 
    title("ABLA")
  )
  points(years_f, stan_data_abla$y[idx], pch = 16, cex = 0.8)
}


#CANO
samples <- util$extract_expectand_vals(fit_CANO)
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_cano$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_cano$y, title ("CANO"))
#Per series 
par(mfrow = c(3, 4))
for (f in 1:G_cano) {
  idx    <- start_idxs_cano[f]:end_idxs_cano[f]
  years_f <- stand_year_cano$year[idx]  
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           
    xlab = paste("Series", f),
    ylab = "y", title ("CANO")
  )
  points(years_f, stan_data_cano$y[idx], pch = 16, cex = 0.8)
}


#PSME
samples <- util$extract_expectand_vals(fit_PSME)
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_psme$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_psme$y, title ("PSME"))
#Per series 
par(mfrow = c(3, 4))
for (f in 1:G_psme) {
  idx    <- start_idxs_psme[f]:end_idxs_psme[f]
  years_f <- stand_year_psme$year[idx]  
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           
    xlab = paste("Series", f),
    ylab = "y", title ("PSME")
  )
  points(years_f, stan_data_psme$y[idx], pch = 16, cex = 0.8)
}


#TSHE
samples <- util$extract_expectand_vals(fit_TSHE)
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_tshe$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 2000, bin_delta = 20,
                         baseline_values = stan_data_tshe$y, title ("TSHE"))
#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_tshe) {
  idx    <- start_idxs_tshe[f]:end_idxs_tshe[f]
  years_f <- stand_year_tshe$year[idx]  
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           
    xlab = paste("Series", f),
    ylab = "y", title ("TSHE")
  )
  points(years_f, stan_data_tshe$y[idx], pch = 16, cex = 0.8)
}

#TSME
samples <- util$extract_expectand_vals(fit_TSME)
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_tsme$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 4000, bin_delta = 20,
                         baseline_values = stan_data_tsme$y, title ("TSME"))
#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_tsme) {
  idx    <- start_idxs_tsme[f]:end_idxs_tsme[f]
  years_f <- stand_year_tsme$year[idx]  
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           
    xlab = paste("Series", f),
    ylab = "y", title ("TSME")
  )
  points(years_f, stan_data_tsme$y[idx], pch = 16, cex = 0.8)
}

#THPL
samples <- util$extract_expectand_vals(fit_THPL)
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_thpl$N, "]")
# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep", title ("THPL"),
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_thpl$y)
#Per series 
par(mfrow = c(3, 4))
for (f in 1:G_thpl) {
  idx    <- start_idxs_thpl[f]:end_idxs_thpl[f]
  years_f <- stand_year_thpl$year[idx]  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           
    xlab = paste("Series", f),
    ylab = "y", title ("THPL")
  )
  points(years_f, stan_data_thpl$y[idx], pch = 16, cex = 0.8)
}


# State ID ----------------------------------------------------------------

#ABAM
samples <- util$extract_expectand_vals(fit_ABAM)
par(mfrow = c(2,1))

#checking for a stand year combination
f <- 1
start_id <- start_idxs_abam [f];
end_id <- end_idxs_abam [f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_abam$y[start_id:end_id], pch = 20)

par(mfrow = c(4, 4))
for (f in 1:length(start_idxs_abam)) {
  
  start_id <- start_idxs_abam[f]
  end_id   <- end_idxs_abam[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_abam$y[start_id:end_id], pch = 20, title("ABAM"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_abam), " -- Press [Enter] for next..."))
}


#ABLA
samples <- util$extract_expectand_vals(fit_ABLA)
par(mfrow = c(3, 4))
for (f in 1:length(start_idxs_abla)) {
  
  start_id <- start_idxs_abla[f]
  end_id   <- end_idxs_abla[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_abla$y[start_id:end_id], pch = 20, title("ABLA"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_abla), " -- Press [Enter] for next..."))
}

#CANO
samples <- util$extract_expectand_vals(fit_CANO)
par(mfrow = c(4,4))
for (f in 1:length(start_idxs_cano)) {
  
  start_id <- start_idxs_cano[f]
  end_id   <- end_idxs_cano[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_cano$y[start_id:end_id], pch = 20, title ("CANO"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_cano), " -- Press [Enter] for next..."))
}


#PSME
samples <- util$extract_expectand_vals(fit_PSME)
par(mfrow = c(4, 4))
for (f in 1:length(start_idxs_psme)) {
  
  start_id <- start_idxs_psme[f]
  end_id   <- end_idxs_psme[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_psme$y[start_id:end_id], pch = 20, title ("PSME"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_psme), " -- Press [Enter] for next..."))
}


#TSHE
samples <- util$extract_expectand_vals(fit_TSHE)
par(mfrow = c(4,4))
for (f in 1:length(start_idxs_tshe)) {
  
  start_id <- start_idxs_tshe[f]
  end_id   <- end_idxs_tshe[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_tshe$y[start_id:end_id], pch = 20, title ("TSHE"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_tshe), " -- Press [Enter] for next..."))
}



#TSME
samples <- util$extract_expectand_vals(fit_TSME)
par(mfrow = c(4, 4))
for (f in 1:length(start_idxs_tsme)) {
  
  start_id <- start_idxs_tsme[f]
  end_id   <- end_idxs_tsme[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_tsme$y[start_id:end_id], pch = 20, title ("TSME"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_tsme), " -- Press [Enter] for next..."))
}


#THPL
samples <- util$extract_expectand_vals(fit_THPL)
par(mfrow = c(4, 4))
for (f in 1:length(start_idxs_thpl)) {
  
  start_id <- start_idxs_thpl[f]
  end_id   <- end_idxs_thpl[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_thpl$y[start_id:end_id], pch = 20, title ("THPL"))
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_thpl), " -- Press [Enter] for next..."))
}





# Prior vs posterior ------------------------------------------------------


fitSP<- fit_ABAM
fitSP<- fit_ABLA
fitSP<- fit_CANO
fitSP<- fit_PSME
fitSP<- fit_TSHE
fitSP<- fit_TSME
fitSP<- fit_THPL

samples <- rstan::extract(fitSP)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2, 1.5),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 0.5, 1),
  grand_logit_theta2   = rnorm(n, -1, 0.5),
  sigma_low_stand      = abs(rnorm(n, 1, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 1.5, 0.7)),
  sigma_theta1_stand   = abs(rnorm(n, 1, 0.3)),
  sigma_theta2_stand   = abs(rnorm(n, 1, 0.5)),
  phi_low              = exp(rnorm(n, log(4), 0.6)),
  phi_high             = exp(rnorm(n, log(4), 0.6))
)


# --- Extract posteriors ---
posteriors <- data.frame(
  mu_log_low           = samples$mu_log_low,
  mu_log_delta         = samples$mu_log_delta,
  grand_logit_theta1   = samples$grand_logit_theta1,
  grand_logit_theta2   = samples$grand_logit_theta2,
  sigma_low_stand      = samples$sigma_low_stand,
  sigma_log_delta      = samples$sigma_log_delta,
  sigma_theta1_stand   = samples$sigma_theta1_stand,
  sigma_theta2_stand   = samples$sigma_theta2_stand,
  phi_low              = exp(samples$log_phi_state[, 1]),
  phi_high             = exp(samples$log_phi_state[, 2])
)
# --- Labels ---
param_labels <- c(
  mu_log_low          = "Grand mean Low state (log)",
  mu_log_delta        = "Diff. between Low/High (log)",
  grand_logit_theta1  = "P stay low (logit)",
  grand_logit_theta2  = "P stay high (logit)",
  sigma_low_stand     = "Low state stand variation",
  sigma_log_delta     = "Stand Fold change variation",
  sigma_theta1_stand  = "Low state trans- Stand var",
  sigma_theta2_stand  = "High state trans- Stand var",
  phi_low             = "Overdisperion (low state)",
  phi_high            = "Overdisperion (high state)"
)

label_order <- c(
  "Grand mean Low state (log)",
  "Diff. between Low/High (log)",
  "Low state stand variation",
  "Stand Fold change variation",
  "P stay low (logit)",
  "P stay high (logit)",
  "Low state trans- Stand var",
  "High state trans- Stand var",
  "Overdisperion (low state)",
  "Overdisperion (high state)"
)

# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
  mutate(
    label = factor(
      param_labels[param],
      levels = label_order),
    type  = factor(type, levels = c("Prior", "Posterior"))
  )
# --- Plot ---
ggplot(df, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_fill_manual(values  = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
  scale_color_manual(values = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
  coord_cartesian(ylim = c(0,2.5))+
  labs(
    title = "Prior vs posterior — single species HMM THPL",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )




# Victor's Graphs ---------------------------------------------------------



#Low seed production state

#Common things in the plot
elevation_df <- tibble(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13","AO03","AG05",
            "AV06","AX15","AB08","PP17","AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600, 668, 700, 850, 1450, 850, 900, 950,
                1060, 1090, 1100, 1150, 1150, 1200, 1450, 1600, 1700, 1800)
)
stand_levels <- elevation_df %>%
  arrange(elevation) %>%
  pull(stand)
ylim_common <- c(-4, 8)



# Selecting the draws from a specific specie
draws <- as.matrix(fit_ABAM, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta_abam <- years_per_series_abam 
abam_summary <- tibble(
  stand = meta_abam$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9))
abam_summary <- abam_summary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
abam_summary <- tibble(stand = stand_levels) %>%
  left_join(abam_summary, by = "stand") %>%
  mutate(x = 1:length(stand_levels))
par(mfrow=c(1,2))
#par(mfrow=c(3,3))
plot(NULL,
     xlim = c(0.5, nrow(abam_summary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "ABAM low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = stand_levels,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(abam_summary$x,
         abam_summary$q10,
         abam_summary$x,
         abam_summary$q90,
         lwd = 1,
         col = "steelblue")
points(abam_summary$x,
       abam_summary$median,
       pch = 19,
       col = "navy")



# Extract draws
draws <- as.matrix(fit_ABLA, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta_abla <- years_per_series_abla 
abla_summary <- tibble(
  stand = meta_abla$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9))
abla_summary <- abla_summary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
abla_summary <- tibble(stand = stand_levels) %>%
  left_join(abla_summary, by = "stand") %>%
  mutate(x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(abla_summary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "ABLa low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = stand_levels,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(abla_summary$x,
         abla_summary$q10,
         abla_summary$x,
         abla_summary$q90,
         lwd = 1,
         col = "steelblue")
points(abla_summary$x,
       abla_summary$median,
       pch = 19,
       col = "navy")


# Extract draws
draws <- as.matrix(fit_CANO, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta <- years_per_series_cano
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9))
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble(stand = stand_levels) %>%
  left_join(Isummary, by = "stand") %>%
  mutate(x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "CANO low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")


# Extract draws
draws <- as.matrix(fit_PSME, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta <- years_per_series_psme
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary<-tibble (stand = stand_levels) %>%
  left_join(Isummary, by = "stand") %>%
  mutate (x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "PSME low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")





# Extract draws
draws <- as.matrix(fit_TSHE, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta <- years_per_series_tshe
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble (stand = stand_levels) %>%
  left_join(Isummary, by= "stand") %>%
  mutate ( x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "TSHE low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")



# Extract draws
draws <- as.matrix(fit_TSME, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta <- years_per_series_tsme
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble (stand = stand_levels)%>%
  left_join (Isummary, by= "stand") %>%
  mutate (x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "TSME low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")



# Extract draws
draws <- as.matrix(fit_THPL, pars = "log_alpha_low")
# Metadata used to fit ABAM
meta <- years_per_series_thpl
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble (stand = stand_levels)%>%
  left_join (Isummary, by= "stand")%>%
  mutate (x= 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_low",
     main = "THPL low-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")



#High seed production state
#Common things in the plot
stand_levels <- elevation_df %>%
  arrange(elevation) %>%
  pull(stand)
ylim_common <- c(-4, 8)


# Selecting the draws from a specific specie
draws <- as.matrix(fit_ABAM, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta_abam <- years_per_series_abam 
abam_summary <- tibble(
  stand = meta_abam$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9))
abam_summary <- abam_summary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
abam_summary <- tibble(stand = stand_levels) %>%
  left_join(abam_summary, by = "stand") %>%
  mutate(x = 1:length(stand_levels))
#par(mfrow=c(3,3))
plot(NULL,
     xlim = c(0.5, nrow(abam_summary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "ABAM high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = stand_levels,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(abam_summary$x,
         abam_summary$q10,
         abam_summary$x,
         abam_summary$q90,
         lwd = 1,
         col = "steelblue")
points(abam_summary$x,
       abam_summary$median,
       pch = 19,
       col = "navy")



# Extract draws
draws <- as.matrix(fit_ABLA, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta_abla <- years_per_series_abla 
abla_summary <- tibble(
  stand = meta_abla$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9))
abla_summary <- abla_summary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
abla_summary <- tibble(stand = stand_levels) %>%
  left_join(abla_summary, by = "stand") %>%
  mutate(x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(abla_summary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "ABLa high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = stand_levels,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(abla_summary$x,
         abla_summary$q10,
         abla_summary$x,
         abla_summary$q90,
         lwd = 1,
         col = "steelblue")
points(abla_summary$x,
       abla_summary$median,
       pch = 19,
       col = "navy")


# Extract draws
draws <- as.matrix(fit_CANO, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta <- years_per_series_cano
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9))
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble(stand = stand_levels) %>%
  left_join(Isummary, by = "stand") %>%
  mutate(x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "CANO high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")


# Extract draws
draws <- as.matrix(fit_PSME, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta <- years_per_series_psme
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary<-tibble (stand = stand_levels) %>%
  left_join(Isummary, by = "stand") %>%
  mutate (x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "PSME high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")





# Extract draws
draws <- as.matrix(fit_TSHE, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta <- years_per_series_tshe
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble (stand = stand_levels) %>%
  left_join(Isummary, by= "stand") %>%
  mutate ( x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "TSHE high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")



# Extract draws
draws <- as.matrix(fit_TSME, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta <- years_per_series_tsme
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble (stand = stand_levels)%>%
  left_join (Isummary, by= "stand") %>%
  mutate (x = 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "TSME high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")



# Extract draws
draws <- as.matrix(fit_THPL, pars = "log_alpha_high")
# Metadata used to fit ABAM
meta <- years_per_series_thpl
Isummary <- tibble(
  stand = meta$stand,
  median = apply(draws, 2, median),
  q10 = apply(draws, 2, quantile, probs = 0.1),
  q90 = apply(draws, 2, quantile, probs = 0.9)
)
Isummary <- Isummary %>%
  left_join(elevation_df, by = "stand") %>%
  arrange(elevation) %>%
  mutate(x = row_number())
Isummary <- tibble (stand = stand_levels)%>%
  left_join (Isummary, by= "stand")%>%
  mutate (x= 1:length(stand_levels))
#par(mfrow=c(1,1))
plot(NULL,
     xlim = c(0.5, nrow(Isummary) + 0.5),
     ylim = ylim_common,
     xaxt = "n",
     xlab = "Stand (low → high elevation)",
     ylab = "log_alpha_high",
     main = "THPL high-state seed production")
axis(1,
     at = 1:length(stand_levels),
     labels = Isummary$stand,
     las = 2)
abline(h = 0, lty = 2, col = "grey70")
segments(Isummary$x,
         Isummary$q10,
         Isummary$x,
         Isummary$q90,
         lwd = 1,
         col = "steelblue")
points(Isummary$x,
       Isummary$median,
       pch = 19,
       col = "navy")


library(tidyverse)

# ── Setup ────────────────────────────────────────────────────────────────────
species_list <- c("ABAM", "ABLA", "CANO", "PSME", "TSHE", "TSME", "THPL")

fits <- list(
  ABAM = fit_ABAM, ABLA = fit_ABLA, CANO = fit_CANO,
  PSME = fit_PSME, TSHE = fit_TSHE, TSME = fit_TSME, THPL = fit_THPL
)

metas <- list(
  ABAM = years_per_series_abam, ABLA = years_per_series_abla, CANO = years_per_series_cano,
  PSME = years_per_series_psme, TSHE = years_per_series_tshe, TSME = years_per_series_tsme,
  THPL = years_per_series_thpl
)

stand_levels <- elevation_df %>% arrange(elevation) %>% pull(stand)
colors       <- setNames(RColorBrewer::brewer.pal(7, "Dark2"), species_list)

# ── Helper: extract summary for one species + one parameter ──────────────────
extract_summary <- function(spp, par) {
  draws <- as.matrix(fits[[spp]], pars = par)
  meta  <- metas[[spp]]
  
  tibble(
    species = spp,
    stand   = meta$stand,
    median  = apply(draws, 2, median),
    q10     = apply(draws, 2, quantile, probs = 0.1),
    q90     = apply(draws, 2, quantile, probs = 0.9)
  ) %>%
    left_join(elevation_df, by = "stand")
}

# ── Build combined data for both states ──────────────────────────────────────
all_low <- map_dfr(species_list, extract_summary, par = "log_alpha_low") %>%
  mutate(stand = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand))

all_high <- map_dfr(species_list, extract_summary, par = "log_alpha_high") %>%
  mutate(stand = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand))

# ── Jitter x positions so species don't overlap ──────────────────────────────
n_spp    <- length(species_list)
offsets  <- setNames(seq(-0.3, 0.3, length.out = n_spp), species_list)

all_low  <- all_low  %>% mutate(x = x_base + offsets[species])
all_high <- all_high %>% mutate(x = x_base + offsets[species])

# ── Plot function ─────────────────────────────────────────────────────────────
plot_state <- function(dat, title, ylim = c(-4, 8)) {
  n_stands <- length(stand_levels)
  
  plot(NULL,
       xlim = c(0.5, n_stands + 0.5),
       ylim = ylim,
       xaxt = "n",
       xlab = "Stand (low → high elevation)",
       ylab = "log alpha",
       main = title)
  
  axis(1, at = 1:n_stands, labels = stand_levels, las = 2, cex.axis = 0.75)
  abline(h = 0, lty = 2, col = "grey70")
  abline(v = 1:n_stands, col = "grey92", lwd = 0.5)  # subtle grid
  
  for (spp in species_list) {
    d <- filter(dat, species == spp)
    segments(d$x, d$q10, d$x, d$q90, lwd = 1.2, col = colors[spp])
    points(d$x, d$median, pch = 19, cex = 0.7, col = colors[spp])
  }
  
  legend("topright",
         legend = species_list,
         col    = colors,
         pch    = 19,
         lwd    = 1.2,
         bty    = "n",
         cex    = 0.8, 
         horiz= TRUE)
}

# ── Render ────────────────────────────────────────────────────────────────────
par(mfrow = c(1, 1), mar = c(6, 4, 3, 1))

plot_state(all_low,  "Low-state seed production — all species")

plot_state(all_high, "High-state seed production — all species")




library(tidyverse)
library(RColorBrewer)

# Plots for Victor
species_list <- c("ABAM", "ABLA", "CANO", "PSME", "TSHE", "TSME", "THPL")

fits <- list(
  ABAM = fit_ABAM, ABLA = fit_ABLA, CANO = fit_CANO,
  PSME = fit_PSME, TSHE = fit_TSHE, TSME = fit_TSME, THPL = fit_THPL
)

metas <- list(
  ABAM = years_per_series_abam, ABLA = years_per_series_abla, CANO = years_per_series_cano,
  PSME = years_per_series_psme, TSHE = years_per_series_tshe, TSME = years_per_series_tsme,
  THPL = years_per_series_thpl
)

stand_levels <- elevation_df %>% arrange(elevation) %>% pull(stand)
colors       <- setNames(RColorBrewer::brewer.pal(7, "Dark2"), species_list)

n_spp   <- length(species_list)
offsets <- setNames(seq(-0.3, 0.3, length.out = n_spp), species_list)

#Absolute log alpha (low or high state) 
extract_summary <- function(spp, par) {
  draws <- as.matrix(fits[[spp]], pars = par)
  meta  <- metas[[spp]]
  tibble(
    species = spp,
    stand   = meta$stand,
    median  = apply(draws, 2, median),
    q10     = apply(draws, 2, quantile, probs = 0.1),
    q90     = apply(draws, 2, quantile, probs = 0.9)
  ) %>%
    left_join(elevation_df, by = "stand")
}

#Extracting what I want
extract_deviation <- function(spp) {
  draws_alpha <- as.matrix(fits[[spp]], pars = "log_alpha_low")
  draws_mu    <- as.matrix(fits[[spp]], pars = "mu_log_low")
  # stand deviation in the baseline for the low
  draws_dev <- sweep(draws_alpha, 1, draws_mu[, 1], FUN = "-")
  meta <- metas[[spp]]
  tibble(
    species = spp,
    stand   = meta$stand,
    median  = apply(draws_dev, 2, median),
    q10     = apply(draws_dev, 2, quantile, probs = 0.1),
    q90     = apply(draws_dev, 2, quantile, probs = 0.9)
  ) %>%
    left_join(elevation_df, by = "stand")
}
extract_deviation2 <- function(spp) {
  draws_alpha <- as.matrix(fits[[spp]], pars = "log_alpha_high")
  draws_mu_delta    <- as.matrix(fits[[spp]], pars = "mu_log_delta")
  draws_mu_low <- as.matrix(fits[[spp]], pars= "mu_log_low")
  #  stand deviation in the baseline for the high
  draws_dev <- sweep(draws_alpha, 1, draws_mu_low [,1]+draws_mu_delta[, 1], FUN = "-")
  meta <- metas[[spp]]
  tibble(
    species = spp,
    stand   = meta$stand,
    median  = apply(draws_dev, 2, median),
    q10     = apply(draws_dev, 2, quantile, probs = 0.1),
    q90     = apply(draws_dev, 2, quantile, probs = 0.9)
  ) %>%
    left_join(elevation_df, by = "stand")
}
extract_delta_deviation <-function(spp){
  draws_sigma <-as.matrix(fits[[spp]], pars= "sigma_log_delta")
  draws_tilde <-as.matrix(fits[[spp]],pars= "log_delta_tilde")
  #stand deviation in the gap
  draws_dev <-draws_sigma[,1]* draws_tilde
  meta<-metas[[spp]]
  tibble(
    species = spp,
    stand   = meta$stand,
    median  = apply(draws_dev, 2, median),
    q10     = apply(draws_dev, 2, quantile, probs = 0.1),
    q90     = apply(draws_dev, 2, quantile, probs = 0.9)
  ) %>%
    left_join(elevation_df, by = "stand")
  
}

#builsing the datasets
all_low <- map_dfr(species_list, extract_summary, par = "log_alpha_low") %>%
  mutate(stand  = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand),
         x      = x_base + offsets[species])

all_high <- map_dfr(species_list, extract_summary, par = "log_alpha_high") %>%
  mutate(stand  = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand),
         x      = x_base + offsets[species])

all_dev <- map_dfr(species_list, extract_deviation) %>%
  mutate(stand  = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand),
         x      = x_base + offsets[species])

all_dev2 <- map_dfr(species_list, extract_deviation2) %>%
  mutate(stand  = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand),
         x      = x_base + offsets[species])

all_delta_dev <- map_dfr(species_list, extract_delta_deviation) %>%
  mutate(stand  = factor(stand, levels = stand_levels),
         x_base = as.numeric(stand),
         x      = x_base + offsets[species])


#Plotting a certain way
plot_state <- function(dat, title, ylab = "log alpha", ylim = c(-4, 8)) {
  n_stands <- length(stand_levels)
  
  plot(NULL,
       xlim = c(0.5, n_stands + 0.5),
       ylim = ylim,
       xaxt = "n",
       xlab = "Stand (low → high elevation)",
       ylab = ylab,
       main = title)
  
  axis(1, at = 1:n_stands, labels = stand_levels, las = 2, cex.axis = 0.75)
  abline(h = 0, lty = 2, col = "grey70")
  abline(v = 1:n_stands, col = "grey92", lwd = 0.5)
  
  for (spp in species_list) {
    d <- filter(dat, species == spp)
    segments(d$x, d$q10, d$x, d$q90, lwd = 1.2, col = colors[spp])
    points(d$x, d$median, pch = 19, cex = 0.7, col = colors[spp])
  }
  
  legend("top",
         legend = species_list,
         col    = colors,
         pch    = 19,
         lwd    = 1.2,
         bty    = "n",
         cex    = 0.8,
         horiz  = TRUE)
}

#Results
par(mfrow = c(1, 1), mar = c(6, 4, 4, 1))

plot_state(all_low,
           title = "Low-state seed production — all species",
           ylab  = "log alpha (low)",
           ylim  = c(-4, 8))

plot_state(all_high,
           title = "High-state seed production — all species",
           ylab  = "log alpha (high)",
           ylim  = c(-4, 8))

plot_state(all_dev,
           title = "Stand deviations from grand mean (low state) — all species",
           ylab  = "log alpha (low) - mu_log_low",
           ylim  = c(-8, 4))

plot_state(all_dev2,
           title = "Stand deviations from grand mean (high state) — all species",
           ylab  = "log alpha (high) - mu_log_low - mu_log_delta",
           ylim  = c(-6, 6))


plot_state(all_delta_dev,
           title = "Stand deviations in gap between states",
           ylab  = "sigma_log_delta * log_delta_tilde",
           ylim  = c(-4, 4))



#New plot for putting all togther
plot_dev_both <- function(dat_low, dat_high, title, ylim = c(-8, 4)) {
  n_stands <- length(stand_levels)
  
  plot(NULL,
       xlim = c(0.5, n_stands + 0.5),
       ylim = ylim,
       xaxt = "n",
       xlab = "Stand (low → high elevation)",
       ylab = "Stand deviation from grand mean",
       main = title)
  
  axis(1, at = 1:n_stands, labels = stand_levels, las = 2, cex.axis = 0.75)
  abline(h = 0, lty = 2, col = "grey70")
  abline(v = 1:n_stands, col = "grey92", lwd = 0.5)
  
  # small extra offset to separate low vs high within each species
  state_offset <- 0.08
  
  for (spp in species_list) {
    d_low  <- filter(dat_low,  species == spp)
    d_high <- filter(dat_high, species == spp)
    
    # low state: circles, solid
    segments(d_low$x - state_offset, d_low$q10,
             d_low$x - state_offset, d_low$q90,
             lwd = 1.2, col = colors[spp])
    points(d_low$x - state_offset, d_low$median,
           pch = 19, cex = 0.7, col = colors[spp])
    
    # high state: triangles, slightly transparent
    segments(d_high$x + state_offset, d_high$q10,
             d_high$x + state_offset, d_high$q90,
             lwd = 1.2, col = adjustcolor(colors[spp], alpha.f = 0.5))
    points(d_high$x + state_offset, d_high$median,
           pch = 17, cex = 0.7, col = adjustcolor(colors[spp], alpha.f = 0.5))
  }
  
  # species legend
  legend("top",
         legend = species_list,
         col    = colors,
         pch    = 19,
         lwd    = 1.2,
         bty    = "n",
         cex    = 0.8,
         horiz  = TRUE)
  
  # state legend
  legend("bottomleft",
         legend = c("Low state", "High state"),
         pch    = c(19, 17),
         col    = c("grey40", adjustcolor("grey40", alpha.f = 0.5)),
         lwd    = 1.2,
         bty    = "n",
         cex    = 0.8)
}

#Results
par(mfrow = c(1, 1), mar = c(6, 4, 4, 1))

plot_dev_both(
  all_dev,
  all_dev2,
  title = "Stand deviations from grand mean — low (circle) vs high (triangle)",
  ylim  = c(-8, 4)
)

