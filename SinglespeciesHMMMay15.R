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
  area       = stand_year_tsme$area
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
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_abam,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_ABLA <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_abla,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_CANO <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_cano,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_PSME <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_psme,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_TSHE <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_tshe,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_TSME <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_tsme,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


fit_THPL <- stan(
  file    = "Stan_code/Species_Stan_Model/SinglespeciesPartialPooling2_delta_NCbisFULL.stan",
  data    = stan_data_thpl,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# Diagnostic plots per species --------------------------------------------


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
divergent <- get_sampler_params(fit_ABAM2, inc_warmup = FALSE) |>
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



# General diagnostics -----------------------------------------------------
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

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












#ABAM
samples <- rstan::extract(fit_ABAM)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM ABAM",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )



#ABLA
samples <- rstan::extract(fit_ABLA)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM ABLA",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )



#CANO
samples <- rstan::extract(fit_CANO)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM CANO",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )



#PSME
samples <- rstan::extract(fit_PSME)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM PSME",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )



#TSHE
samples <- rstan::extract(fit_TSHE)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM TSHE",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )



#TSME
samples <- rstan::extract(fit_TSME)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM TSME",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )


#THPL
samples <- rstan::extract(fit_THPL)  
n <- length(samples$mu_log_low)
# --- Simulate priors ---
priors <- data.frame(
  mu_log_low           = rnorm(n, 2.6, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.0, 0.7),
  grand_logit_theta2   = rnorm(n, 0.0, 0.7),
  sigma_low_stand      = abs(rnorm(n, 0, 0.5)),
  sigma_log_delta      = abs(rnorm(n, 0, 1.0)),
  sigma_theta1_stand   = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand   = abs(rnorm(n, 0, 0.7)),
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
  mu_log_low          = "µ log low (grand mean, low state)",
  mu_log_delta        = "µ log delta (fold-change)",
  grand_logit_theta1  = "logit θ₁ (P stay low)",
  grand_logit_theta2  = "logit θ₂ (P stay high)",
  sigma_low_stand     = "σ low stand",
  sigma_log_delta     = "σ log delta stand",
  sigma_theta1_stand  = "σ θ₁ stand",
  sigma_theta2_stand  = "σ θ₂ stand",
  phi_low             = "φ low state",
  phi_high            = "φ high state"
)
# --- Reshape ---
df <- bind_rows(
  priors     %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Prior"),
  posteriors %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>% mutate(type = "Posterior")
) %>%
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
    title = "Prior vs posterior — single species HMM THPL",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
