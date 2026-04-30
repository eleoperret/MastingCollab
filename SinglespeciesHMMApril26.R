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

stand_year_all <- stand_year_all %>%
  filter(!(stand == "SUNR" & year == 2024))


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
  area       = stand_year_abam$area
)

print(stan_data_abam)

#to verify
ABAM_data<-stand_year_all%>%
  filter(spp=="ABAM")%>%
  filter(year != 2009) %>% 
  group_by(stand, year) %>%
  summarise(
    total_viable_sds    = sum(y, na.rm = TRUE),
    area = sum(area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stand, year)
  
unique(ABAM_data$stand);length(unique(ABAM_data$stand))
print(ABAM_data$total_viable_sds); length(ABAM_data$total_viable_sds)
#works

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
  area       = stand_year_abla$area
)

print(stan_data_abla)

#to verify
ABLA_data<-seed_filtered%>%
  filter(spp=="ABLA")%>%
  filter(year != 2009) %>% 
  group_by(stand, year) %>%
  summarise(
    total_viable_sds    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stand, year)

unique(ABLA_data$stand);length(unique(ABLA_data$stand))
print(ABLA_data$total_viable_sds); length(ABLA_data$total_viable_sds)
#works
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
  area       = stand_year_cano$area
)

print(stan_data_cano)

#to verify
CANO_data<-seed_filtered%>%
  filter(spp=="CANO")%>%
  filter(year != 2009) %>% 
  group_by(stand, year) %>%
  summarise(
    total_viable_sds    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stand, year)

unique(CANO_data$stand);length(unique(CANO_data$stand))
print(CANO_data$total_viable_sds); length(CANO_data$total_viable_sds)
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
  area       = stand_year_psme$area
)

print(stan_data_psme)

#to verify
PSME_data<-seed_filtered%>%
  filter(spp=="PSME")%>%
  filter(year != 2009) %>% 
  group_by(stand, year) %>%
  summarise(
    total_viable_sds    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stand, year)

unique(PSME_data$stand);length(unique(PSME_data$stand))
print(PSME_data$total_viable_sds); length(PSME_data$total_viable_sds)
#works

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
  area       = stand_year_tshe$area
)

print(stan_data_tshe)

#to verify
TSHE_data<-seed_filtered%>%
  filter(spp=="TSHE")%>%
  filter(year != 2009) %>% 
  group_by(stand, year) %>%
  summarise(
    total_viable_sds    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stand, year)

unique(TSHE_data$stand);length(unique(TSHE_data$stand))
print(TSHE_data$total_viable_sds); length(TSHE_data$total_viable_sds)
#works

# THPL  --------------------------------------------------------------------
#Per species
stand_year_thpl <- stand_year_all %>%
  filter(spp== "THPL")

# Creating my stan data list
years_per_series_thpl <- stand_year_thpl %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)))  # NEW

N_stands_thpl <- length(unique(years_per_series_tshe$stand_id))  # NEW

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
  area       = stand_year_thpl$area
)

print(stan_data_thpl)

#to verify
THPL_data<-seed_filtered%>%
  filter(spp=="THPL")%>%
  filter(year != 2009) %>% 
  group_by(stand, year) %>%
  summarise(
    total_viable_sds    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(stand, year)

unique(THPL_data$stand);length(unique(THPL_data$stand))
print(THPL_data$total_viable_sds); length(THPL_data$total_viable_sds)
#works

# Fitting Model -----------------------------------------------------------

fit_ABAM <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbis.stan",
  data    = stan_data_abam,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_ABLA <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbis.stan",
  data    = stan_data_abla,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_CANO <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbis.stan",
  data    = stan_data_cano,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_PSME <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbis.stan",
  data    = stan_data_psme,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_TSHE <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbis.stan",
  data    = stan_data_tshe,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_THPL <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbis.stan",
  data    = stan_data_thpl,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# Diagnostic plots per species --------------------------------------------


species_list <- c("ABAM", "ABLA", "CANO", "PSME", "TSHE", "THPL")

fits_list <- list(
  ABAM = fit_ABAM,
  ABLA = fit_ABLA,
  CANO = fit_CANO,
  PSME = fit_PSME,
  TSHE = fit_TSHE,
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


# Diagnostic plots ABAM --------------------------------------------------------


#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_ABAM, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_ABAM, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")

#Model : Multispecies Paritial Pooling Delta Phi 6

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_ABAM)
util$check_all_hmc_diagnostics(diagnostics)
#there is 1% of divergence

samples <- util$extract_expectand_vals(fit_ABAM)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_stand[',   1:N_stands_abam, ']'),
                                         paste0('alpha_theta2_stand[',   1:N_stands_abam, ']'),
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_stand',
                                          # Phi 
                                          paste0('log_phi_state[1]'),
                                          paste0('log_phi_state[2]'),
                                          
                                          # Low state means
                                          'grand_mean_low','sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)




print(fit_ABAM, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))



#Extracting some infos
draws       <- rstan::extract(fit_ABAM)
state_draws <- draws$state



#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_abam$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_abam$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_abam$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_abam) {
  idx    <- start_idxs_abam[f]:end_idxs_abam[f]
  years_f <- stand_year_abam$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_abam$y[idx], pch = 16, cex = 0.8)
}






par(mfrow = c(4, 4))

for (f in 1:G_abam) {
  idx     <- start_idxs_abam[f]:end_idxs_abam[f]
  T_f     <- length(idx)
  years_f <- stand_year_abam$year[idx]
  
  util$plot_disc_pushforward_quantiles(
    samples,
    paste0("y_rep[", idx, "]")
  )
  points(x = 1:T_f, y = stan_data_abam$y[idx], pch = 16, cex = 0.8)
  
  # Add year labels on x axis
  axis(1, at = 1:T_f, labels = years_f, las = 2, cex.axis = 0.6)
  
  title(main = paste0(years_per_series_abam$stand[f],
                      " (", min(years_f), "-", max(years_f), ")"),
        cex.main = 0.8)
}

par(mfrow = c(1, 1))










# State identification 

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




for (f in 1:length(start_idxs_abam)) {
  
  start_id <- start_idxs_abam[f]
  end_id   <- end_idxs_abam[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_abam$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_abam), " -- Press [Enter] for next..."))
}



# Prior vs posterior 


#partial pooling effect: 
samples <- rstan::extract(fit_ABAM)

# Observed mean per stand
obs_mean <- sapply(1:length(start_idxs_abam), function(f) {
  idx <- start_idxs_abam[f]:end_idxs_abam[f]
  mean(stan_data_abam$y[idx] / stan_data_abam$area[idx])
})

# Posterior median of low and high state means per stand
post_low  <- apply(exp(samples$log_alpha_low),  2, median)
post_high <- apply(exp(samples$log_alpha_high), 2, median)

# Plot
par(mfrow = c(1, 2))

plot(obs_mean, post_low,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "Low state", pch = 19, col = "steelblue")
text(obs_mean, post_low, labels = unique(ABAM_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

plot(obs_mean, post_high,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "High state", pch = 19, col = "darkorange")
text(obs_mean, post_high, labels = unique(ABAM_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)






# Diagnostic plots ABLA --------------------------------------------------------

#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_ABLA, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_ABLA, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_ABLA)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_ABLA)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_stand[',   1:N_stands_abla, ']'),
                                         paste0('alpha_theta2_stand[',   1:N_stands_abla, ']'),
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_stand',
                                          # Phi 
                                          paste0('log_phi_state[1]'),
                                          paste0('log_phi_state[2]'),
                                          
                                          # Low state means
                                          'grand_mean_low','sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)

#Extracting some infos
draws       <- rstan::extract(fit_ABLA)
state_draws <- draws$state



#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_abla$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_abla$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_abla$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_abla) {
  idx    <- start_idxs_abla[f]:end_idxs_abla[f]
  years_f <- stand_year_abla$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_abla$y[idx], pch = 16, cex = 0.8)
}



# print(fit_ABLA, pars = c("grand_mean_low",
#                          "log_delta_high_grand_mean", 
#                          "sigma_theta2_stand",
#                          "sigma_low_stand",
#                          "sigma_log_delta_high_stand",
#                          "log_phi_state"))

print(fit_ABLA, pars = c("mu_log_low",
                         "mu_log_delta",
                         "sigma_theta1_stand",
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta",
                         "log_phi_state"))

# State identification 

samples <- util$extract_expectand_vals(fit_ABLA)
par(mfrow = c(2,1))

#checking for a stand year combination
f <- 1
start_id <- start_idxs_abla [f];
end_id <- end_idxs_abla [f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_abla$y[start_id:end_id], pch = 20)




for (f in 1:length(start_idxs_abla)) {
  
  start_id <- start_idxs_abla[f]
  end_id   <- end_idxs_abla[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_abla$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_abla), " -- Press [Enter] for next..."))
}


par(mfrow = c(4, 4))
for (f in 1:G_abla) {
  idx     <- start_idxs_abla[f]:end_idxs_abla[f]
  T_f     <- length(idx)
  years_f <- stand_year_abla$year[idx]
  
  util$plot_disc_pushforward_quantiles(
    samples,
    paste0("y_rep[", idx, "]")
  )
  points(x = 1:T_f, y = stan_data_abla$y[idx], pch = 16, cex = 0.8)
  
  # Add year labels on x axis
  axis(1, at = 1:T_f, labels = years_f, las = 2, cex.axis = 0.6)
  
  title(main = paste0(years_per_series_abla$stand[f],
                      " (", min(years_f), "-", max(years_f), ")"),
        cex.main = 0.8)
}

par(mfrow = c(1, 1))



# Prior vs posterior 


#partial pooling effect: 
samples <- rstan::extract(fit_ABLA)

# Observed mean per stand
obs_mean <- sapply(1:length(start_idxs_abla), function(f) {
  idx <- start_idxs_abla[f]:end_idxs_abla[f]
  mean(stan_data_abla$y[idx] / stan_data_abla$area[idx])
})

# Posterior median of low and high state means per stand
post_low  <- apply(exp(samples$log_alpha_low),  2, median)
post_high <- apply(exp(samples$log_alpha_high), 2, median)

# Plot
par(mfrow = c(1, 2))

plot(obs_mean, post_low,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "Low state", pch = 19, col = "steelblue")
text(obs_mean, post_low, labels = unique(ABLA_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

plot(obs_mean, post_high,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "High state", pch = 19, col = "darkorange")
text(obs_mean, post_high, labels = unique(ABLA_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)





# Diagnostic plots CANO --------------------------------------------------------

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
pairs(fit_CANO, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")

#checking where the divergence are still from (after changing the prior in : SinglespeciesPartialPooling_deltaPhiNC2)
pairs(fit_CANO, pars = c("grand_mean_low", 
                         "log_delta_high_grand_mean",
                         "sigma_log_delta_high_stand",
                         "sigma_low_stand"))

print(fit_CANO, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))
##

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_CANO)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_CANO)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_stand[',   1:N_stands_cano, ']'),
                                         paste0('alpha_theta2_stand[',   1:N_stands_cano, ']'),
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_stand',
                                          # Phi 
                                          paste0('log_phi_state[1]'),
                                          paste0('log_phi_state[2]'),
                                          
                                          # Low state means
                                          'grand_mean_low','sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)


print(fit_CANO, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))

#Extracting some infos
draws       <- rstan::extract(fit_CANO)
state_draws <- draws$state



#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_cano$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_cano$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_cano$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_cano) {
  idx    <- start_idxs_cano[f]:end_idxs_cano[f]
  years_f <- stand_year_cano$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_cano$y[idx], pch = 16, cex = 0.8)
}



# State identification 

samples <- util$extract_expectand_vals(fit_CANO)
par(mfrow = c(2,1))

#checking for a stand year combination
f <- 1
start_id <- start_idxs_cano [f];
end_id <- end_idxs_cano [f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_cano$y[start_id:end_id], pch = 20)




for (f in 1:length(start_idxs_cano)) {
  
  start_id <- start_idxs_cano[f]
  end_id   <- end_idxs_cano[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_cano$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_cano), " -- Press [Enter] for next..."))
}


par(mfrow = c(4, 4))
for (f in 1:G_cano) {
  idx     <- start_idxs_cano[f]:end_idxs_cano[f]
  T_f     <- length(idx)
  years_f <- stand_year_cano$year[idx]
  
  util$plot_disc_pushforward_quantiles(
    samples,
    paste0("y_rep[", idx, "]")
  )
  points(x = 1:T_f, y = stan_data_cano$y[idx], pch = 16, cex = 0.8)
  
  # Add year labels on x axis
  axis(1, at = 1:T_f, labels = years_f, las = 2, cex.axis = 0.6)
  
  title(main = paste0(years_per_series_cano$stand[f],
                      " (", min(years_f), "-", max(years_f), ")"),
        cex.main = 0.8)
}

par(mfrow = c(1, 1))

# Prior vs posterior 



#partial pooling effect: 
samples <- rstan::extract(fit_CANO)

# Observed mean per stand
obs_mean <- sapply(1:length(start_idxs_cano), function(f) {
  idx <- start_idxs_cano[f]:end_idxs_cano[f]
  mean(stan_data_cano$y[idx] / stan_data_cano$area[idx])
})

# Posterior median of low and high state means per stand
post_low  <- apply(exp(samples$log_alpha_low),  2, median)
post_high <- apply(exp(samples$log_alpha_high), 2, median)

# Plot
par(mfrow = c(1, 2))

plot(obs_mean, post_low,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "Low state", pch = 19, col = "steelblue")
text(obs_mean, post_low, labels = unique(CANO_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

plot(obs_mean, post_high,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "High state", pch = 19, col = "darkorange")
text(obs_mean, post_high, labels = unique(CANO_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)



# Diagnostic plots PSME --------------------------------------------------------

#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_PSME, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_PSME, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_PSME)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_PSME)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_stand[',   1:N_stands_psme, ']'),
                                         paste0('alpha_theta2_stand[',   1:N_stands_psme, ']'),
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_stand',
                                          # Phi 
                                          paste0('log_phi_state[1]'),
                                          paste0('log_phi_state[2]'),
                                          
                                          # Low state means
                                          'grand_mean_low','sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)


print(fit_PSME, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))


#Extracting some infos
draws       <- rstan::extract(fit_PSME)
state_draws <- draws$state



#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_psme$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_psme$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_psme$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_psme) {
  idx    <- start_idxs_psme[f]:end_idxs_psme[f]
  years_f <- stand_year_psme$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_psme$y[idx], pch = 16, cex = 0.8)
}



# State identification 

samples <- util$extract_expectand_vals(fit_PSME)
par(mfrow = c(2,1))

#checking for a stand year combination
f <- 1
start_id <- start_idxs_psme [f];
end_id <- end_idxs_psme [f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_psme$y[start_id:end_id], pch = 20)



tcheck<-PSME_data %>% mutate(f = row_number()) %>% filter(stand %in% c("AM16", "AR07"))



for (f in 1:length(start_idxs_psme)) {
  
  start_id <- start_idxs_psme[f]
  end_id   <- end_idxs_psme[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_psme$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_psme), " -- Press [Enter] for next..."))
}

#all in one
par(mfrow = c(4, 4))
for (f in 1:G_psme) {
  idx     <- start_idxs_psme[f]:end_idxs_psme[f]
  T_f     <- length(idx)
  years_f <- stand_year_psme$year[idx]
  
  util$plot_disc_pushforward_quantiles(
    samples,
    paste0("y_rep[", idx, "]")
  )
  points(x = 1:T_f, y = stan_data_psme$y[idx], pch = 16, cex = 0.8)
  
  # Add year labels on x axis
  axis(1, at = 1:T_f, labels = years_f, las = 2, cex.axis = 0.6)
  
  title(main = paste0(years_per_series_psme$stand[f],
                      " (", min(years_f), "-", max(years_f), ")"),
        cex.main = 0.8)
}

par(mfrow = c(1, 1))

# Prior vs posterior 



#partial pooling effect: 
samples <- rstan::extract(fit_PSME)

# Observed mean per stand
obs_mean <- sapply(1:length(start_idxs_psme), function(f) {
  idx <- start_idxs_psme[f]:end_idxs_psme[f]
  mean(stan_data_psme$y[idx] / stan_data_psme$area[idx])
})

# Posterior median of low and high state means per stand
post_low  <- apply(exp(samples$log_alpha_low),  2, median)
post_high <- apply(exp(samples$log_alpha_high), 2, median)

# Plot
par(mfrow = c(1, 2))

plot(obs_mean, post_low,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "Low state", pch = 19, col = "steelblue")
text(obs_mean, post_low, labels = unique(PSME_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

plot(obs_mean, post_high,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "High state", pch = 19, col = "darkorange")
text(obs_mean, post_high, labels = unique(PSME_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)


# Diagnostic plots TSHE --------------------------------------------------------

#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_TSHE, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_TSHE, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_TSHE)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_TSHE)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_stand[',   1:N_stands_tshe, ']'),
                                         paste0('alpha_theta2_stand[',   1:N_stands_tshe, ']'),
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_stand',
                                          # Phi 
                                          paste0('log_phi_state[1]'),
                                          paste0('log_phi_state[2]'),
                                          
                                          # Low state means
                                          'grand_mean_low','sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)



print(fit_TSHE, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))

#Extracting some infos
draws       <- rstan::extract(fit_TSHE)
state_draws <- draws$state



#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_tshe$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_tshe$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_tshe$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_tshe) {
  idx    <- start_idxs_tshe[f]:end_idxs_tshe[f]
  years_f <- stand_year_tshe$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_tshe$y[idx], pch = 16, cex = 0.8)
}



# State identification 

samples <- util$extract_expectand_vals(fit_TSHE)
par(mfrow = c(2,1))

#checking for a stand year combination
f <- 1
start_id <- start_idxs_tshe [f];
end_id <- end_idxs_tshe [f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_tshe$y[start_id:end_id], pch = 20)




for (f in 1:length(start_idxs_tshe)) {
  
  start_id <- start_idxs_tshe[f]
  end_id   <- end_idxs_tshe[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_tshe$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_tshe), " -- Press [Enter] for next..."))
}

#all in one
par(mfrow = c(4, 4))
for (f in 1:G_tshe) {
  idx     <- start_idxs_tshe[f]:end_idxs_tshe[f]
  T_f     <- length(idx)
  years_f <- stand_year_tshe$year[idx]
  
  util$plot_disc_pushforward_quantiles(
    samples,
    paste0("y_rep[", idx, "]")
  )
  points(x = 1:T_f, y = stan_data_tshe$y[idx], pch = 16, cex = 0.8)
  
  # Add year labels on x axis
  axis(1, at = 1:T_f, labels = years_f, las = 2, cex.axis = 0.6)
  
  title(main = paste0(years_per_series_tshe$stand[f],
                      " (", min(years_f), "-", max(years_f), ")"),
        cex.main = 0.8)
}

par(mfrow = c(1, 1))

# Prior vs posterior 




#partial pooling effect: 
samples <- rstan::extract(fit_TSHE)

# Observed mean per stand
obs_mean <- sapply(1:length(start_idxs_tshe), function(f) {
  idx <- start_idxs_tshe[f]:end_idxs_tshe[f]
  mean(stan_data_tshe$y[idx] / stan_data_tshe$area[idx])
})

# Posterior median of low and high state means per stand
post_low  <- apply(exp(samples$log_alpha_low),  2, median)
post_high <- apply(exp(samples$log_alpha_high), 2, median)

# Plot
par(mfrow = c(1, 2))

plot(obs_mean, post_low,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "Low state", pch = 19, col = "steelblue")
text(obs_mean, post_low, labels = unique(TSHE_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

plot(obs_mean, post_high,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "High state", pch = 19, col = "darkorange")
text(obs_mean, post_high, labels = unique(TSHE_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

# Diagnostic plots THPL --------------------------------------------------------

#If there is a divergence : 

#where is this divergence?
divergent <- get_sampler_params(fit_THPL, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)

cat("Divergence in chain:", divergent$chain, "\n")
cat("At draw:", divergent$draw, "\n")

#Checking the pair plots : 
pairs(fit_THPL, 
      pars = c("sigma_theta1_stand", "sigma_theta2_stand", "sigma_low_stand",  "sigma_log_delta"),
      condition = "accept_stat__")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_THPL)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_THPL)

base_samples <- util$filter_expectands(samples,
                                       c(paste0('alpha_theta1_stand[',   1:N_stands_thpl, ']'),
                                         paste0('alpha_theta2_stand[',   1:N_stands_thpl, ']'),
                                         'sigma_theta1_stand',   'sigma_theta2_stand',
                                         'grand_logit_theta1',   'grand_logit_theta2',
                                         'grand_mean_low',       'log_delta_high_grand_mean'
                                       ))
util$check_all_expectand_diagnostics(base_samples)

base_samples2 <- util$filter_expectands(samples,
                                        c(# Delta parameters
                                          'log_delta_high_grand_mean',
                                          'sigma_log_delta_high_stand',
                                          # Phi 
                                          paste0('log_phi_state[1]'),
                                          paste0('log_phi_state[2]'),
                                          
                                          # Low state means
                                          'grand_mean_low','sigma_low_stand'
                                        ))

util$check_all_expectand_diagnostics(base_samples2)


print(fit_THPL, pars = c("grand_mean_low",
                         "log_delta_high_grand_mean", 
                         "sigma_theta2_stand",
                         "sigma_low_stand",
                         "sigma_log_delta_high_stand",
                         "log_phi_state"))

#Extracting some infos
draws       <- rstan::extract(fit_THPL)
state_draws <- draws$state



#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_thpl$N, "]")

# Plot PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = stan_data_thpl$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = stan_data_thpl$y)


#Per series 
par(mfrow = c(4, 4))
for (f in 1:G_thpl) {
  idx    <- start_idxs_thpl[f]:end_idxs_thpl[f]
  years_f <- stand_year_thpl$year[idx]  # actual calendar years
  
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  
  util$plot_conn_pushforward_quantiles(
    samples,
    pred_names_f,
    years_f,           # use actual years on x-axis
    xlab = paste("Series", f),
    ylab = "y"
  )
  points(years_f, stan_data_thpl$y[idx], pch = 16, cex = 0.8)
}



# State identification 

samples <- util$extract_expectand_vals(fit_THPL)
par(mfrow = c(2,1))

#checking for a stand year combination
f <- 1
start_id <- start_idxs_thpl [f];
end_id <- end_idxs_thpl [f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_thpl$y[start_id:end_id], pch = 20)




for (f in 1:length(start_idxs_thpl)) {
  
  start_id <- start_idxs_thpl[f]
  end_id   <- end_idxs_thpl[f]
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_thpl$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(start_idxs_thpl), " -- Press [Enter] for next..."))
}


#all in one
par(mfrow = c(4, 4))
for (f in 1:G_thpl) {
  idx     <- start_idxs_thpl[f]:end_idxs_thpl[f]
  T_f     <- length(idx)
  years_f <- stand_year_thpl$year[idx]
  
  util$plot_disc_pushforward_quantiles(
    samples,
    paste0("y_rep[", idx, "]")
  )
  points(x = 1:T_f, y = stan_data_thpl$y[idx], pch = 16, cex = 0.8)
  
  # Add year labels on x axis
  axis(1, at = 1:T_f, labels = years_f, las = 2, cex.axis = 0.6)
  
  title(main = paste0(years_per_series_thpl$stand[f],
                      " (", min(years_f), "-", max(years_f), ")"),
        cex.main = 0.8)
}

par(mfrow = c(1, 1))


# Prior vs posterior



#partial pooling effect: 
samples <- rstan::extract(fit_THPL)

# Observed mean per stand
obs_mean <- sapply(1:length(start_idxs_thpl), function(f) {
  idx <- start_idxs_thpl[f]:end_idxs_thpl[f]
  mean(stan_data_thpl$y[idx] / stan_data_thpl$area[idx])
})

# Posterior median of low and high state means per stand
post_low  <- apply(exp(samples$log_alpha_low),  2, median)
post_high <- apply(exp(samples$log_alpha_high), 2, median)

# Plot
par(mfrow = c(1, 2))

plot(obs_mean, post_low,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "Low state", pch = 19, col = "steelblue")
text(obs_mean, post_low, labels = unique(THPL_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)

plot(obs_mean, post_high,
     xlab = "Observed mean density", ylab = "Posterior median",
     main = "High state", pch = 19, col = "darkorange")
text(obs_mean, post_high, labels = unique(THPL_data$stand), pos = 3, cex = 0.7)
abline(0, 1, lty = 2)
