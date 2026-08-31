##Code for multispecies HMM model 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

##Next things to do : 

#1 : Prepare plots for manuscript
  #Figure of the time series of seed production overall
  #Figure of the trace plots of chain convergence
  #Fifure of state probability for every stands and species
  #Figure with posteriror distribution of parameters (1:1?)
#2: Think about how to calculate the interspecies synchrony
#3: 


#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)
library(tidyverse)
library(posterior)


library(purrr)


options(mc.cores = parallel::detectCores())




#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

list.files()
list.files("Data/")
list.files("Rscript_AnalysisEleonore/")

#Adding Mike's code
setwd("Rscript_AnalysisEleonore/")
util <- new.env()
source("mcmc_visualization_tools.R", local=util)
source("mcmc_analysis_tools_rstan.R", local=util)
ls(util)

setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

# Data import + modification -------------------------------------------------------------
seed_data<-read.csv("Data/SeedData_all.csv")

# 1) Keeping species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","ABLA","CANO","PSME","TSHE","TSME","THPL"))
unique(seed_data$stand)


# 2) I want only possible combination of stands (ie. for CANO I don't want it at low elevation stands)
# a) Extra stand ( to add manually as there is no map/records of the species present at SPRY, PARA and CANO as it was not installed by the PSP crew but Janneke)
# Add SPRY, SUNR, PARA for CANO, ABAM, ABLA
extra_stands <- seed_data %>%
  filter(spp %in% c("CANO", "ABAM", "ABLA"),
         stand %in% c("SPRY", "SUNR", "PARA")) %>%
  filter(year != 2009) %>%
  group_by(spp, stand, year) %>%
  summarise(y    = sum(total_viable_sds, na.rm = TRUE),area = sum(size, na.rm = TRUE),.groups = "drop") %>%
  arrange(spp, stand, year)
# Checking
extra_stands %>%
  group_by(spp, stand) %>%
  summarise(n_years = n(), mean_density = mean(y/area), .groups = "drop")

# b) Getting the information from the PSP data
mapRainier<-read.csv("data/rawdata/Cleaned_mapping_2017 Rainier.csv", sep=";")

#Species per stand if present 
stands_per_species <- mapRainier %>%
  group_by(species) %>%  
  summarise(stands = list(unique(stand_id)),.groups = "drop") 

# Saving the object (I prefer this way but can be done also directly)
saveRDS(stands_per_species, file = "data/stands_per_species.rds")

#c) Adding only the possible combination of species x stands
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R (orignal script but put there here for better readability)
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)
#Seeds_filtered contains now all the species x stand combination 

#d)Stand elevation
stand_elev <- data.frame(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13",
            "AO03","AG05","AV06","AX15","AB08","PP17",
            "AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600,668,700,850,1450,850,
                900,950,1060,1090,1100,1150,
                1150,1200,1450,1600,1700,1800)
)
saveRDS(stand_elev, file= "data/stand_elevation_table.rds")


# Data prep for Stan ---------------------------------------------------------------

#Removing year 2009, and merging all seed counts
stand_year_all <- seed_filtered %>%
  filter(year != 2009) %>% #because like this they start all the same year
  group_by(spp, stand, year) %>%
  summarise(y= sum(total_viable_sds, na.rm = TRUE),area = sum(size, na.rm = TRUE),.groups = "drop") %>%
  arrange(spp, stand, year)

#Adding the data from SUNR, SPRY and PARA
stand_year_all <- bind_rows(stand_year_all, extra_stands) %>%
  arrange(spp, stand, year)

#As I have a year missing for SUNR. 
missing_rows <- tibble(
  spp   = c("ABAM", "ABLA", "CANO"),
  stand = "SUNR",
  year  = 2023,
  y     = NA_integer_,
  area  = NA_real_
)

#Adding the missing year to the dataset
stand_year_all <- bind_rows(stand_year_all, missing_rows) %>%
  group_by(spp, stand) %>%
  mutate(area = ifelse(is.na(area), mean(area, na.rm = TRUE), area)) %>%
  ungroup() %>%
  arrange(spp, stand, year)

#Checking where I have missing "y" (sum of seed count)
obs_missing <- as.integer(is.na(stand_year_all$y))
#Replacing those NA values by 0 as I cannot have NA values in my dataset
stand_year_all <- stand_year_all %>%
  mutate(y = ifelse(is.na(y), 0L, as.integer(y)))

#Defining my stan list
# Species index
stand_year_all$species_id <- as.numeric(as.factor(stand_year_all$spp))
S <- length(unique(stand_year_all$spp)) #which is 7 at the moment

# Checking my species order before running
print(levels(as.factor(stand_year_all$spp)))

# Creating my stan data list
years_per_series <- stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(T_i = n(), .groups = "drop")%>%
  mutate(stand_id = as.numeric(as.factor(stand)),
         species_id = as.numeric(as.factor(spp)) #NEW
         ) 

N_stands <- length(unique(years_per_series$stand_id))  

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
  N_stands   = N_stands, 
  stand_id   = years_per_series$stand_id,  
  species_id = years_per_series$species_id, 
  start_idxs = start_idxs,
  end_idxs   = end_idxs,
  y          = stand_year_all$y,
  area       = stand_year_all$area,
  obs_missing = obs_missing
)

print (stan_data_all)

# Fitting Model -----------------------------------------------------------
#Current model used after discussion in git
fit_75 <- stan(
  file    = "Stan_code_Eleonore/Species_Stan_Model/MultispeciesStanCodeLast.stan",
  data    = stan_data_all,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# Model diagnostic --------------------------------------------------------

summary(fit_75)
print(fit_75)

diagnostics <- util$extract_hmc_diagnostics(fit_75)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_75)

names(samples)
base_samples1 <- util$filter_expectands(samples,c("rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2"))
base_samples2 <- util$filter_expectands(samples,c(paste0('alpha_theta1_species_nc[',1:S, ']'),paste0('alpha_theta2_species_nc[',   1:S, ']'), "sigma_theta1_species","sigma_theta2_species"))
base_samples3 <- util$filter_expectands(samples,c("mu_log_low",paste0('alpha_low_species[',1:S, ']'),paste0('alpha_stand[',   1:S*N_stands, ']'), "sigma_low_species","sigma_stand"))
base_samples4<-util$filter_expectands(samples,c("mu_log_delta",paste0('log_delta_species_nc[',1:S, ']')),"sigma_log_delta_species","sigma_delta_stand")
base_samples5 <- util$filter_expectands(samples,c(paste0('log_phi_species[',1:S, ']'),"log_phi_state[1]","log_phi_state[2]"))
util$check_all_expectand_diagnostics(base_samples1);util$check_all_expectand_diagnostics(base_samples2)

base_samples <- util$filter_expectands(c(paste0('alpha_theta1_species_nc[',1:S, ']'),paste0('alpha_theta2_species_nc[',   1:S, ']'),"rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2",'sigma_theta1_species','sigma_theta2_species','grand_logit_theta1','grand_logit_theta2','mu_log_low'))

base_samples <- util$filter_expectands(c(paste0('alpha_theta1_species_nc[',1:S, ']'),paste0('alpha_theta2_species_nc[',   1:S, ']'),"rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2",'sigma_theta1_species','sigma_theta2_species','grand_logit_theta1','grand_logit_theta2','mu_log_low'))

# PPC ---------------------------------------------------------------------

#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")

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

#Single observation example
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


# Prior vs posterior ------------------------------------------------

samples <- rstan::extract(fit_75)
n <- nrow(samples$mu_log_low)

priors <- data.frame(
  mu_log_low                = rnorm(n, 2, 1.5),
  mu_log_delta              = rnorm(n, 1.5, 1.5),
  grand_logit_theta1        = rnorm(n, 0.5, 1),
  grand_logit_theta2        = rnorm(n, -1, 0.5),
  sigma_low_species         = abs(rnorm(n, 0.5, 1)),
  sigma_stand               = abs(rnorm(n,1,0.5)),
  sigma_log_delta_species   = abs(rnorm(n, 1.5, 0.7)),
  sigma_delta_stand         = abs(rnorm(n,1,0.5)),
  sigma_theta1_species      = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_species      = abs(rnorm(n, 0, 0.7)),
  phi_low                   = exp(rnorm(n,0,0.3)),
  phi_high                  = exp(rnorm(n,0,0.3))
)

posteriors <- data.frame(
  mu_log_low                = samples$mu_log_low,
  mu_log_delta              = samples$mu_log_delta,
  grand_logit_theta1        = samples$grand_logit_theta1,
  grand_logit_theta2        = samples$grand_logit_theta2,
  sigma_low_species         = samples$sigma_low_species,
  sigma_stand               = samples$sigma_stand,
  sigma_log_delta_species   = samples$sigma_log_delta_species,
  sigma_delta_stand         = samples$sigma_delta_stand,
  sigma_theta1_species      = samples$sigma_theta1_species,
  sigma_theta2_species      = samples$sigma_theta2_species,
  phi_low                   = exp(samples$log_phi_state[,1]),
  phi_high                  = exp(samples$log_phi_state[,2])
)

#labels 
param_labels <- c(
  mu_log_low              = "low state mean",
  mu_log_delta            = "fold-change",
  grand_logit_theta1      = "P stay low",
  grand_logit_theta2      = "P stay high",
  sigma_low_species       = "Species variation (low)",
  sigma_stand             = "Stands variation (low)",
  sigma_log_delta_species = "Species variation (delta)",
  sigma_delta_stand       = "Stands variation (delta)",
  sigma_theta1_species    = "Species transition variation (k=1)",
  sigma_theta2_species    = "Species transition variation (k=2)",
  phi_species             = "φ species",   
  phi_low                 = "φ low state ",
  phi_high                = "φ high state "
)

#Reshape to long 
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

#Plot
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



# Model results -----------------------------------------------------------

print(levels(as.factor(stand_year_all$spp)))

species<- c("ABAM" ,"ABLA" ,"CANO" ,"PSME", "THPL" ,"TSHE" ,"TSME")

results <-list ()

for (s in 1:S){
  #seed count
  low_mean <- exp(samples$mu_log_low + samples[[paste0("alpha_low_species[",s,"]")]])
  delta <- samples$mu_log_delta + samples[[paste0("log_delta_species[",s,"]")]]
  fold_change <- 1 + exp(delta)
  mean_high <- low_mean * fold_change
  
  #transition probability
  theta_1<- plogis(samples$grand_logit_theta1 + samples[[paste0("alpha_theta1_species[",s,"]")]])
  theta_2<- plogis(samples$grand_logit_theta2 + samples[[paste0("alpha_theta2_species[",s,"]")]])
  
  #overdispersion
  phi_species <- exp(samples[[paste0("log_phi_species[",s,"]")]])


results[[s]] <- data.frame(
  species            = species[s],
  low_state_mean     = median(low_mean),
  high_state_fold_change = median(fold_change),
  high_state_mean    = median(mean_high),
  p_stay_low         = median(theta_1),
  p_stay_high        = median(theta_2),
  sigma_low_specie   = median(samples$sigma_low_species),
  sigma_log_delta    = median(samples$sigma_log_delta_species),
  sigma_theta1_stand = median(samples$sigma_theta1_species),
  sigma_theta2_stand = median(samples$sigma_theta2_species),
  phi_low            = median(exp(samples[[paste0("log_phi_state[",1,"]")]])),
  phi_high           = median(exp(samples[[paste0("log_phi_state[",2,"]")]])),
  phi_species        = median (phi_species)
)
}

results

summary<-do.call(rbind, results)
summary


# Intraspecific synchrony -------------------------------------------------


#pulling ths state for each draws
state_draws <- rstan::extract(fit_75, pars = "state")$state

meta_all <- stand_year_all %>%
  mutate(row_id = row_number()) %>%
  rename(species = spp)

get_unit_year_state <- function(draw_idx) {
  tibble(
    species = meta_all$species,
    stand   = meta_all$stand,
    year    = meta_all$year,
    state   = state_draws[draw_idx, ]
  )
}

# Victor's synchrony. Within one year, for one draw: the fraction of stands in whichever state is more common that year 
majority_state <- function(df) {
  state_counts <- table(df$state)
  tibble(majority = max(state_counts) / sum(state_counts))
}

summarise_draws <- function(df, value_col, group_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      median = median(.data[[value_col]], na.rm = TRUE),
      lo90   = quantile(.data[[value_col]], 0.05, na.rm = TRUE),
      hi90   = quantile(.data[[value_col]], 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

#ABAM
target_species <- "ABAM"
abam_meta <- meta_all %>%
  filter(species == target_species)

abam_state_draws <- state_draws[, abam_meta$row_id]

n_draws_total <-nrow(state_draws) #or abam state draws but I keep it like this so I can reuse it for the next species
draw_idx_use<- seq(1,n_draws_total, by=4)

abam_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

abam_summary <- summarise_draws(abam_draws, "majority", "year")

ggplot(abam_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = NULL, y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "ABAM"
  ) +
  theme_minimal(base_size = 13)


#ABLA
target_species <- "ABLA"
target_meta <- meta_all %>%
  filter(species == target_species)

target_state_draws <- state_draws[, target_meta$row_id]

n_draws_total <-nrow(target_state_draws) 

abla_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

target_summary <- summarise_draws(abla_draws, "majority", "year")

ggplot(target_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  labs(
    x = "year", y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "ABLA"
  ) +
  theme_minimal(base_size = 13)


#CANO
target_species <- "CANO"
target_meta <- meta_all %>%
  filter(species == target_species)

target_state_draws <- state_draws[, target_meta$row_id]

n_draws_total <-nrow(target_state_draws) 

cano_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

target_summary <- summarise_draws(cano_draws, "majority", "year")

ggplot(target_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  labs(
    x = "year", y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "CANO"
  ) +
  theme_minimal(base_size = 13)


#PSME
target_species <- "PSME"
target_meta <- meta_all %>%
  filter(species == target_species)

target_state_draws <- state_draws[, target_meta$row_id]

n_draws_total <-nrow(target_state_draws) 

psme_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

target_summary <- summarise_draws(psme_draws, "majority", "year")

ggplot(target_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = "year", y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "PSME"
  ) +
  theme_minimal(base_size = 13)


#TSHE
target_species <- "TSHE"
target_meta <- meta_all %>%
  filter(species == target_species)

target_state_draws <- state_draws[, target_meta$row_id]

n_draws_total <-nrow(target_state_draws) 

tshe_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

target_summary <- summarise_draws(tshe_draws, "majority", "year")

ggplot(target_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  labs(
    x = "year", y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "TSHE"
  ) +
  theme_minimal(base_size = 13)


#TSME
target_species <- "TSME"
target_meta <- meta_all %>%
  filter(species == target_species)

target_state_draws <- state_draws[, target_meta$row_id]

n_draws_total <-nrow(target_state_draws) 

tsme_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

target_summary <- summarise_draws(tsme_draws, "majority", "year")

ggplot(target_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = "year", y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "TSME"
  ) +
  theme_minimal(base_size = 13)


#THPL
target_species <- "THPL"
target_meta <- meta_all %>%
  filter(species == target_species)

target_state_draws <- state_draws[, target_meta$row_id]

n_draws_total <-nrow(target_state_draws) 

thpl_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})

target_summary <- summarise_draws(thpl_draws, "majority", "year")

ggplot(target_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = "year", y = "Posterior probability of being in the majority state (%)",
    title = "Percentage of stands within the majority state per year",
    subtitle = "THPL"
  ) +
  theme_minimal(base_size = 13)
