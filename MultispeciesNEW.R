##Code for multispecies HMM model 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

##Next things to do : 

#1 : Add the 2009 year
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

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)


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

obs_missing <- as.integer(is.na(stand_year_all$y))

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


#Final Stan data list #added a stand parameter 
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

fit_all <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesNEWMikeFull_10062026.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_multi <- stan(
  file    = "Stan_code/Species_Stan_Model/Multispecies_05062026.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# Diagnostic plots Centered --------------------------------------------------------


#If there is a divergence : 
#where is this divergence?
divergent <- get_sampler_params(fit_all2, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)
##Add the divergence color 

pairs(
  fit_all2,
  pars = c(
    "log_phi_state[1]",
    "log_phi_state[2]"
  ),
  condition = "divergent__"
)

pairs(
  fit_all2,
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
  fit_all2,
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
  fit_all2,
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
#If there is a divergence : 
#where is this divergence?

divergent <- get_sampler_params(fit_all, inc_warmup = FALSE)
lapply(divergent, function(x) sum(x[,"divergent__"]))


# Extract chain 4 samples with divergence info
sampler_params_ch4 <- get_sampler_params(fit_all2, inc_warmup = FALSE)[[4]]

# Find which iteration diverged
div_iter <- which(sampler_params_ch4[, "divergent__"] == 1)
cat("Divergence at iteration:", div_iter, "\n")

# Check the energy and stepsize at that iteration
sampler_params_ch4[div_iter, ]

divergent <- get_sampler_params(fit_all2, inc_warmup = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "chain") |>
  mutate(draw = row_number()) |>
  filter(divergent__ == 1)
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
    "log_phi_species[3]",
    "log_phi_species[2]"
  ),
  condition = "divergent__"
)


pairs(
  fit_all2,
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
  fit_all2,
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


print(fit_all2, pars = c("grand_logit_theta1", "grand_logit_theta2",
                        "mu_log_low", "mu_log_delta",
                        "sigma_theta1_species", "sigma_theta2_species",
                        "sigma_low_species", "sigma_log_delta_species",
                        "sigma_theta1_stand", "sigma_theta2_stand",
                        "sigma_low_stand", "sigma_log_delta_stand",
                        "log_phi_state"))


pairs(fit_all, pars=c("grand_logit_theta1",
                  "grand_logit_theta2",
                  "mu_log_delta"))

pairs(fit_all,
      pars=c("mu_log_delta",
             "log_phi[1,1]",
             "log_phi[1,2]"))
mcmc_trace(
  as.array(fit_all),
  pars = c(
    "mu_log_low",
    "mu_log_delta",
    "grand_logit_theta1",
    "grand_logit_theta2",
    "sigma_log_delta_species",
    "sigma_log_delta_stand"
  )
)

draws <- as_draws_df(fit_all)

draws |>
  dplyr::group_by(.chain) |>
  dplyr::summarise(
    mu_delta = mean(mu_log_delta),
    theta1 = mean(grand_logit_theta1),
    theta2 = mean(grand_logit_theta2)
  )

mcmc_rank_overlay(
  as.array(fit_all),
  pars = c(
    "mu_log_delta",
    "grand_logit_theta1",
    "sigma_log_delta_species"
  )
)

summ <- summarise_draws(as_draws(fit_all))

summ |>
  dplyr::arrange(desc(rhat)) |>
  dplyr::select(variable, rhat, ess_bulk, ess_tail) |>
  head(30)


pairs(
  fit_all,
  pars = c(
    "grand_logit_theta1",
    "grand_logit_theta2",
    "mu_log_delta",
    "log_phi[1,1]",
    "log_phi[1,2]"
  ),
  condition = "chain"
)

mcmc_trace(
  as.array(fit_all),
  pars = c(
    "log_phi[1,1]",
    "log_phi[1,2]"
  )
)

mcmc_scatter(
  as.array(fit_all),
  pars = c("grand_logit_theta1", "mu_log_delta"),
  np = nuts_params(fit_all)
)


mcmc_scatter(
  as.array(fit_all),
  pars = c("mu_log_delta", "log_phi[1,2]"),
  np = nuts_params(fit_all)
)

draws <- as_draws_df(fit_all)

draws %>%
  mutate(
    prop_state2 = rowMeans(
      select(., starts_with("state[")) == 2
    )
  ) %>%
  group_by(.chain) %>%
  summarise(mean_occ = mean(prop_state2))

# Results -----------------------------------------------------------------

summary(fit_all)

s <- summary(fit_all)$summary

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

##1) checking the fit of my model
diagnostics <- util$extract_hmc_diagnostics(fit_all)
util$check_all_hmc_diagnostics(diagnostics)
#there is 1% of divergence

samples <- util$extract_expectand_vals(fit_all)

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
                                          'log_phi'
                                        ),check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples2)


util$plot_pairs_by_chain(samples[["log_alpha_low[20]"]], "log_alpha_low[20]",
                         samples[["log_alpha_high[20]"]], "log_alpha_high[20]")

# PPC ---------------------------------------------------------------------

#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")

samples <- util$extract_expectand_vals(fit_all)

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
samples <- util$extract_expectand_vals(fit_all)
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



samples <- rstan::extract(fit_all)

n <- length(samples$mu_log_low)

priors <- data.frame(
  mu_log_low           = rnorm(n, 2.8, 1.0),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 1.2, 0.7),
  grand_logit_theta2   = rnorm(n, 0.5, 0.7),
  sigma_low_species    = abs(rnorm(n, 0.5, 1)),
  sigma_low_stand      = abs(rnorm(n, 0.5, 1)),
  sigma_log_delta_sp   = abs(rnorm(n, 0, 1)),
  sigma_log_delta_stand= abs(rnorm(n, 0, 1)),
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




# Posteriror vs prior multioverdispersion ---------------------------------

samples <- rstan::extract(fit_all)
n <- length(samples$mu_log_low)
species_names <- c("ABAM", "ABLA", "CANO", "PSME", "TSHE", "TSME", "THPL")

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

# Priors are the same for all species
priors <- data.frame(
  mu_log_low            = rnorm(n, 2.8, 1.0),
  mu_log_delta          = rnorm(n, 1.5, 1.5),
  grand_logit_theta1    = rnorm(n, 1.2, 0.7),
  grand_logit_theta2    = rnorm(n, 0.5, 0.7),
  sigma_low_species     = abs(rnorm(n, 0.5, 1)),
  sigma_low_stand       = abs(rnorm(n, 0.5, 1)),
  sigma_log_delta_sp    = abs(rnorm(n, 0, 1)),
  sigma_log_delta_stand = abs(rnorm(n, 0, 1)),
  sigma_theta1_species  = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_species  = abs(rnorm(n, 0, 0.7)),
  sigma_theta1_stand    = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_stand    = abs(rnorm(n, 0, 0.7)),
  phi_low               = exp(rnorm(n, log(4), 0.6)),
  phi_high              = exp(rnorm(n, log(4), 0.6))
)

prior_long <- priors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(type = "Prior")

# Loop over species and save one plot per species
for (sp_idx in 1:length(species_names)) {
  
  posteriors <- data.frame(
    mu_log_low            = samples$mu_log_low,
    mu_log_delta          = samples$mu_log_delta,
    grand_logit_theta1    = samples$grand_logit_theta1,
    grand_logit_theta2    = samples$grand_logit_theta2,
    sigma_low_species     = samples$sigma_low_species,
    sigma_low_stand       = samples$sigma_low_stand,
    sigma_log_delta_sp    = samples$sigma_log_delta_species,
    sigma_log_delta_stand = samples$sigma_log_delta_stand,
    sigma_theta1_species  = samples$sigma_theta1_species,
    sigma_theta2_species  = samples$sigma_theta2_species,
    sigma_theta1_stand    = samples$sigma_theta1_stand,
    sigma_theta2_stand    = samples$sigma_theta2_stand,
    phi_low               = exp(samples$log_phi[, sp_idx, 1]),  # species-specific
    phi_high              = exp(samples$log_phi[, sp_idx, 2])   # species-specific
  )
  
  post_long <- posteriors %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    mutate(type = "Posterior")
  
  df <- bind_rows(prior_long, post_long) %>%
    mutate(
      label = param_labels[param],
      type  = factor(type, levels = c("Prior", "Posterior"))
    )
  
  p <- ggplot(df, aes(x = value, fill = type, color = type)) +
    geom_density(alpha = 0.4, linewidth = 0.7) +
    facet_wrap(~ label, scales = "free", ncol = 3) +
    scale_fill_manual(values  = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
    scale_color_manual(values = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
    labs(
      title = paste("Prior vs posterior —", species_names[sp_idx]),
      x = NULL, y = "Density", fill = NULL, color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "top",
      strip.text       = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}


# Chain subset retrodictive check -----------------------------------------

#subset
subset_chains <- function(expectand_vals, chains) {
  expectand_vals[chains, , drop = FALSE]
}

plot_retrodictive_by_chain_group <- function(
    f,                                     
    samples,                                
    stan_data_all,                              
    chain_groups = list("Chains 1&3" = c(1, 3),
                        "Chains 2&4" = c(2, 4))) {
  
  t_idx   <- stan_data_all$start_idxs[f]:stan_data_all$end_idxs[f]
  sp_f    <- stan_data_all$sp[stan_data_all$start_idxs[f]]
  stand_f <- stan_data_all$stand_id[f]
  y_obs   <- stan_data_all$y[t_idx]
  T_len   <- length(t_idx)
  years   <- seq_along(t_idx)            
  
  par(mfrow = c(length(chain_groups), 1),
      mar   = c(4, 4, 3, 1))
  
  for (grp_name in names(chain_groups)) {
    chains <- chain_groups[[grp_name]]
    
    # y rep quantiles
    yrep_lo  <- numeric(T_len)
    yrep_med <- numeric(T_len)
    yrep_hi  <- numeric(T_len)
    
    for (i in seq_along(t_idx)) {
      t      <- t_idx[i]
      col_nm <- paste0("y_rep[", t, "]")
      vals   <- subset_chains(samples[[col_nm]], chains)   
      flat   <- as.vector(vals)
      yrep_lo[i]  <- quantile(flat, 0.05)
      yrep_med[i] <- quantile(flat, 0.50)
      yrep_hi[i]  <- quantile(flat, 0.95)
    }
    
    # Plot
    ylim_max <- max(yrep_hi, y_obs) * 1.05
    
    plot(years, yrep_med,
         type = "n",
         ylim = c(0, ylim_max),
         xlab = "Time step",
         ylab = "Count",
         main = paste0(grp_name,
                       "  |  f=", f,
                       "  species=", sp_f,
                       "  stand=",   stand_f))
    
    # 90% posterior predictive 
    polygon(c(years, rev(years)),
            c(yrep_lo, rev(yrep_hi)),
            col = util$c_light, border = NA)
    
    # Median prediction
    lines(years, yrep_med, col = util$c_dark, lwd = 2)
    
    # Observed data
    points(years, y_obs, pch = 16, cex = 0.8, col = "black")
  }
}


# samples comes from Mike's extract_expectand_vals()
samples   <- util$extract_expectand_vals(fit_all)

# Plot series f = 20, chains 1&3 vs 2&4
plot_retrodictive_by_chain_group(
  f         = 20,
  samples   = samples,
  stan_data_all = stan_data_all
)

# Loop over all series
for (f in seq_len(stan_data_all$F)) {
  plot_retrodictive_by_chain_group(f, samples, stan_data_all)
}



# Figuring out where the issue is  ----------------------------------------

# Extract the relevant parameters
par(mfrow = c(2, 2))

# Grand mean delta
util$plot_expectand_pushforward(
  samples[["mu_log_delta"]], 30,
  display_name = "mu_log_delta",
  main = "grand mean delta"
)

# Species 1 effect
util$plot_expectand_pushforward(
  samples[["log_delta_species_nc[1]"]], 30,
  display_name = "log_delta_species_nc[1]",
  main = "species 1 delta (non-centered)"
)

# Stand 6 effect
util$plot_expectand_pushforward(
  samples[["log_delta_stand_nc[6]"]], 30,
  display_name = "log_delta_stand_nc[6]",
  main = "stand 6 delta (non-centered)"
)

# sigma for stand delta
util$plot_expectand_pushforward(
  samples[["sigma_log_delta_stand"]], 30,
  display_name = "sigma_log_delta_stand",
  main = "sigma log delta stand"
)



util$plot_pairs_by_chain(
  samples[["log_alpha_low[6]"]],  "log_alpha_low[6]",
  samples[["log_alpha_high[6]"]], "log_alpha_high[6]"
)



# log_alpha_high - log_alpha_low = the implied log gap
# On the count scale that's exp(log_alpha_high) / exp(log_alpha_low)

delta_vals <- samples[["log_alpha_high[6]"]] - samples[["log_alpha_low[6]"]]

util$plot_expectand_pushforward(
  delta_vals, 30,
  display_name = "log_alpha_high[6] - log_alpha_low[6]",
  main = "implied log gap for f=6"
)


util$plot_expectand_pushforward(
  samples[["alpha_low_stand[6]"]], 30,
  display_name = "alpha_low_stand[6]",
  main = "stand 6 baseline effect"
)



# Compute the ratio of upper ribbon between chain groups
# as a systematic way to find the most bimodal series
bimodality_score <- numeric(stan_data_all$F)

for (f in seq_len(stan_data_all$F)) {
  t_idx <- stan_data_all$start_idxs[f]:stan_data_all$end_idxs[f]
  
  hi_chains13 <- mean(sapply(t_idx, function(t)
    quantile(as.vector(subset_chains(samples[[paste0("y_rep[",t,"]")]], c(1,3))), 0.95)))
  
  hi_chains24 <- mean(sapply(t_idx, function(t)
    quantile(as.vector(subset_chains(samples[[paste0("y_rep[",t,"]")]], c(2,4))), 0.95)))
  
  bimodality_score[f] <- abs(log(hi_chains24 / hi_chains13))
}

# Most problematic series
order(bimodality_score, decreasing = TRUE)[1:10]
# [1]  6 71 18 72 22 21 30 31 13 27

# What stand does f=71 belong to?
stan_data_all$stand_id[6] #6
stan_data_all$stand_id[71]#12
stan_data_all$stand_id[18]#2
stan_data_all$stand_id[72]#17
stan_data_all$stand_id[22]#14
stan_data_all$stand_id[21]#13
stan_data_all$stand_id[30]#11
stan_data_all$stand_id[31]#12
stan_data_all$stand_id[13]#13
stan_data_all$stand_id[27]#8

par(mfrow = c(3, 4))

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[6]"]], 30,
  display_name = "alpha_low_stand[6]",
  main = "Stand 6"
)


util$plot_expectand_pushforward(
  samples[["alpha_low_stand[17]"]], 30,
  display_name = "alpha_low_stand[72]",
  main = "Stand 17"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[2]"]], 30,
  display_name = "alpha_low_stand[18]",
  main = "Stand 2"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[14]"]], 30,
  display_name = "alpha_low_stand[22]",
  main = "Stand 14"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[13]"]], 30,
  display_name = "alpha_low_stand[21]",
  main = "Stand 13"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[11]"]], 30,
  display_name = "alpha_low_stand[30]",
  main = "Stand 11"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[12]"]], 30,
  display_name = "alpha_low_stand[31]",
  main = "Stand 12"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[8]"]], 30,
  display_name = "alpha_low_stand[27]",
  main = "Stand 8"
)

util$plot_expectand_pushforward(
  samples[["alpha_low_stand[13]"]], 30,
  display_name = "alpha_low_stand[13]",
  main = "Stand 13"
)


draws <- as_draws_df(fit_all)

chain_13 <- subset(draws, chain %in% c(1,3))
chain_24 <- subset(draws, chain %in% c(2,4))

summary_by_group <- function(param) {
  c(
    c13 = mean(chain_13[[param]]),
    c24 = mean(chain_24[[param]])
  )
}

summary_by_group("log_phi[2,1]")
summary_by_group("log_phi[2,2]")
summary_by_group("log_delta_species[2]")



# Checking chains ---------------------------------------------------------
# ── helpers ──────────────────────────────────────────────────────────────────
get_chain_draws <- function(fit_all, param, chain_ids) {
  all_chains <- as.array(fit_all)[, , param]  # iter x chain
  all_chains[, chain_ids, drop = FALSE]
}

# ── identify series (f) for species 1 ────────────────────────────────────────
# assumes you have your data list in memory as stan_data
sp_vec      <- stan_data_all$sp
start_idxs  <- stan_data_all$start_idxs
end_idxs    <- stan_data_all$end_idxs
stand_id    <- stan_data_all$stand_id

f_sp2 <- which(sp_vec[start_idxs] == 2)   # series indices for species 1

# ── extract y_rep draws (iterations x chains x N) ───────────────────────────
y_rep_array <- as.array(fit_all, pars = "y_rep")  # iter x chain x N
# dim: [iter, chain, N]

# ── build plot data ───────────────────────────────────────────────────────────
chain_pairs <- list(c(1, 3), c(2, 4))
pair_labels <- c("Chains 1&3", "Chains 2&4")

plot_list <- list()

for (f in f_sp2) {
  s_id   <- stand_id[f]
  t_idx  <- start_idxs[f]:end_idxs[f]
  t_local <- seq_along(t_idx)
  y_obs  <- stan_data_all$y[t_idx]
  
  for (p in seq_along(chain_pairs)) {
    chains <- chain_pairs[[p]]
    # subset to these chains and these time points
    # y_rep_array dims: [iter, chain, N] -> flatten iter x chain for these t
    draws <- y_rep_array[, chains, t_idx, drop = FALSE]  
    # reshape to (iter*chain) x length(t_idx)
    draws_mat <- matrix(draws, 
                        nrow = dim(draws)[1] * length(chains), 
                        ncol = length(t_idx))
    
    ribbon_df <- data.frame(
      t      = t_local,
      lo     = apply(draws_mat, 2, quantile, 0.025),
      hi     = apply(draws_mat, 2, quantile, 0.975),
      med    = apply(draws_mat, 2, median),
      y_obs  = y_obs,
      pair   = pair_labels[p],
      stand  = s_id,
      series = f
    )
    plot_list[[length(plot_list) + 1]] <- ribbon_df
  }
}

plot_df <- bind_rows(plot_list)

# ── plot: one page per stand, two panels (chain pairs) ───────────────────────
par(mfrow= c(4,4))
stands_sp2 <- stand_id[f_sp2]

for (st in unique(stands_sp2)) {
  df_st <- filter(plot_df, stand == st)
  
  p <- ggplot(df_st, aes(x = t)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.25) +
    geom_line(aes(y = med), colour = "firebrick", linewidth = 0.8) +
    geom_point(aes(y = y_obs), size = 1.5) +
    facet_wrap(~pair, ncol = 1, scales = "free_y") +
    labs(
      title = paste0("Species 2 | Stand ", st),
      x = "Time step", y = "Count"
    ) +
    theme_bw(base_size = 12)
  
  print(p)
  # or save: ggsave(paste0("stand_", st, "_sp1.png"), p, width=7, height=6)
}

par(mfrow= c(4,4))
ggplot(filter(plot_df, grepl("1&3|2&4", pair)), aes(x = t)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.25) +
  geom_line(aes(y = med), colour = "firebrick", linewidth = 0.8) +
  geom_point(aes(y = y_obs), size = 1.5) +
  facet_wrap(stand ~ pair, scales = ("free_y"), ncol = 4) +
  coord_cartesian(ylim= c(0,600))+
  labs(x = "Time step", y = "Count", title = "Species 2 — all stands") +
  theme_bw(base_size = 8)

# Checking bimodality  ----------------------------------------------------
plot_bimodal_check <- function(fit_all, par, warmup_frac = 0.5) {
  arr      <- as.array(fit_all)   # iter x chain x param
  n_iter   <- dim(arr)[1]
  n_chains <- dim(arr)[2]
  warmup   <- floor(n_iter * warmup_frac)
  
  df <- do.call(rbind, lapply(seq_len(n_chains), function(ch) {
    data.frame(
      value  = arr[, ch, par],
      iter   = seq_len(n_iter),
      phase  = ifelse(seq_len(n_iter) <= warmup, "warmup", "sampling"),
      chain  = paste0("Chain ", ch)
    )
  }))
  
  # colour: grey = warmup, dark red = sampling
  chain_colours <- c("Chain 1" = "#E41A1C", "Chain 2" = "#377EB8",
                     "Chain 3" = "#4DAF4A", "Chain 4" = "#984EA3")
  
  # top: density per chain (sampling only)
  p_density <- df %>%
    filter(phase == "sampling") %>%
    ggplot(aes(x = value, colour = chain, fill = chain)) +
    geom_density(alpha = 0.15, linewidth = 0.8) +
    scale_colour_manual(values = chain_colours) +
    scale_fill_manual(values = chain_colours) +
    labs(title = par, x = par, y = "Density") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  # bottom: trace coloured by chain
  p_trace <- ggplot(df, aes(x = iter, y = value, colour = chain)) +
    geom_line(alpha = 0.5, linewidth = 0.3) +
    geom_vline(xintercept = warmup, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = chain_colours) +
    labs(x = "Iteration", y = par, colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  
  p_density / p_trace   # patchwork stacking
}

library(patchwork)

# one parameter at a time
plot_bimodal_check(fit_all, "mu_log_low")
plot_bimodal_check(fit_all, "mu_log_delta")
plot_bimodal_check(fit_all, "grand_logit_theta1")
plot_bimodal_check(fit_all, "grand_logit_theta2")

# or loop over a set and save
key_params <- c(
  "mu_log_low", "mu_log_delta",
  "grand_logit_theta1", "grand_logit_theta2",
  "sigma_low_species", "sigma_log_delta_species"
)

for (par in key_params) {
  p <- plot_bimodal_check(fit_all, par)
  print(p)
  # ggsave(paste0("bimodal_", par, ".png"), p, width = 7, height = 6)
}















plot_bimodal_scatter <- function(fit_all, par, warmup_frac = 0.5) {
  arr      <- as.array(fit_all)   # iter x chain x param
  n_iter   <- dim(arr)[1]
  n_chains <- dim(arr)[2]
  warmup   <- floor(n_iter * warmup_frac)
  
  df <- do.call(rbind, lapply(seq_len(n_chains), function(ch) {
    data.frame(
      value = arr[, ch, par],
      iter  = seq_len(n_iter),
      chain = paste0("Chain ", ch)
    )
  }))
  
  ggplot(df, aes(x = iter, y = value, colour = iter)) +
    geom_point(size = 0.6, alpha = 0.5, shape = 18) +   # diamond shape like yours
    geom_vline(xintercept = warmup, linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    scale_colour_gradient(low = "grey80", high = "darkred") +
    facet_wrap(~chain, ncol = 2) +
    labs(
      title  = par,
      x      = "Iteration",
      y      = par,
      colour = "Iteration"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold")
    )
}


library(patchwork)

key_params <- c(
  "mu_log_low", "mu_log_delta",
  "grand_logit_theta1", "grand_logit_theta2",
  "sigma_low_species", "sigma_log_delta_species"
)

for (par in key_params) {
  print(plot_bimodal_scatter(fit_all, par))
}












# Rhat is Na --------------------------------------------------------------

#This means that the sampler is stuck with similar values for a or more parameters and that the variance cannot be defined. 

# Extract the summary
fit_summary <- summary(fit_all)$summary

# Find parameters with NA Rhat
na_rhat_params <- rownames(fit_summary)[is.na(fit_summary[, "Rhat"])]
na_rhat_params

stan_data_all$y[c(183, 289, 432)] 

# Also look at zero-variance parameters (n_eff will be NA or 0 too)
fit_summary[is.na(fit_summary[, "Rhat"]), ]


# Extract posterior draws for a suspicious parameter
draws <- rstan::extract(fit_all)

# Check variance of specific parameters - if 0, it's stuck
sapply(draws, function(x) var(as.vector(x)))


# Check if the parameter hits 0, 1, or extreme boundary
for (p in na_rhat_params) {
  cat(p, ":", range(as.vector(rstan::extract(fit_all, p)[[1]])), "\n")
}

# Check your data for NAs, zeros where log() is applied, etc.
lapply(stan_data_all, function(x) sum(is.na(x)))
lapply(stan_data_all, summary)


# Are some species/groups completely missing data?
# (common in multispecies models)
table(stan_data_all$S)  # or however species is indexed

# Check n_eff alongside Rhat - both NA = never sampled properly
fit_summary[fit_summary[,"n_eff"] < 10, ]


# Find continuous params with low ESS (the ones that actually matter)
fit_summary_clean <- fit_summary[!is.na(fit_summary[,"Rhat"]), ]

# Parameters with poor mixing
fit_summary_clean[fit_summary_clean[,"Rhat"] > 1.05, 
                  c("mean","sd","n_eff","Rhat")]

# Lowest ESS continuous params
head(fit_summary_clean[order(fit_summary_clean[,"n_eff"]), 
                       c("mean","sd","n_eff","Rhat")], 20)


# Which species have fewest detections?
table(stan_data_all$sp)
# You already know: species 2 (n=65) and 7 (n=67) are smallest

# Check detections (not just presences) per species
tapply(stan_data_all$y, stan_data_all$sp, sum)


# Extract state draws
state_draws <- rstan::extract(fit_all, "state")[[1]]  # shape: [samples x N]

# Find sites that DO vary (i.e., the sampler sees transitions)
state_sd <- apply(state_draws, 2, sd)
cat("Sites with no state variation:", sum(state_sd == 0), "\n")
cat("Sites with some variation:", sum(state_sd > 0), "\n")

# For the varying sites, what's the mean state?
par(mfrow= c(1,1))
hist(apply(state_draws, 2, mean), 
     main = "Mean state across sites", xlab = "Mean state (1=absent, 2=present)")

# Cross-reference stuck sites with species and year
stuck_indices <- which(state_sd == 0)

# What species do the stuck sites belong to?
table(stan_data_all$sp[stuck_indices])

# What years?
table(stan_data_all$start_idxs[stuck_indices])  # adjust field name

# Any detections at stuck sites?
tapply(stan_data_all$y[stuck_indices], 
       stan_data_all$sp[stuck_indices], sum)

dim(Gamma_draws)  # what are the actual dimensions?
Gamma_draws <- rstan::extract(fit_all, "Gamma")[[1]]
# shape likely: [samples, n_sites, 2, 2]

# Persistence probabilities (staying in same state)
# Gamma[s, 1, 1] = P(stay absent | absent)
# Gamma[s, 2, 2] = P(stay present | present)

# Are persistence probs near 1.0 for all species?
for (f in 1:stan_data_all$F) {
  s <- stan_data_all$sp[stan_data_all$start_idxs[f]]
  cat("Series", f, "- Species", s,
      "- P(stay absent):",  round(mean(Gamma_draws[, f, 1, 1]), 3),
      "- P(stay present):", round(mean(Gamma_draws[, f, 2, 2]), 3), "\n")
}

# How long is the state vector vs. what you'd expect?
length(stan_data_all$y)              # total observations
stan_data_all$stand_id * stan_data_all$F  # expected if state[site,year]
# (adjust field names to yours)


# Correct check: how many unique series do you have?
cat("Total observations:", stan_data_all$N, "\n")               # 950
cat("Number of stands:", stan_data_all$N_stands, "\n")
cat("Number of species:", stan_data_all$S, "\n")
cat("Number of years (F):", stan_data_all$F, "\n")

# How long is each time series?
series_lengths <- stan_data_all$end_idxs - stan_data_all$start_idxs + 1
summary(series_lengths)
table(series_lengths)  # are all series 15 years? or variable length?

# How many series total?
length(stan_data_all$start_idxs)


# Check missing data pattern for stuck series
stuck_start_positions <- c(87, 399, 449, 490, 516, 635, 661, 
                           702, 713, 728, 799, 821, 847, 884, 914)

for (s in stuck_start_positions) {
  end_s <- stan_data_all$end_idxs[which(stan_data_all$start_idxs == s)]
  missing_rate <- mean(stan_data_all$obs_missing[s:end_s])
  y_sum <- sum(stan_data_all$y[s:end_s])
  cat("start:", s, "- missing rate:", round(missing_rate, 2), 
      "- total detections:", y_sum, "\n")
}


#what series are causing the issues
# Which series contains observations 785-839?
for (f in 1:stan_data_all$F) {
  s <- stan_data_all$start_idxs[f]
  e <- stan_data_all$end_idxs[f]
  if (any(785:839 >= s & 785:839 <= e)) {
    cat("Series", f, "- species", stan_data_all$sp[s], 
        "- stand", stan_data_all$stand_id[f],
        "- obs", s, "to", e, "\n")
  }
}


gap <- colMeans(log_alpha_high_draws) - colMeans(log_alpha_low_draws)
summary(exp(gap))  # compare to before: was min=8x, median=19x


mu_delta_draws <- rstan::extract(fit_all, "mu_log_delta")[[1]]
cat("mu_log_delta posterior mean:", mean(mu_delta_draws), "\n")
cat("Implied grand mean count ratio:", mean(exp(log1p(exp(mu_delta_draws)))), "\n")

# And the full gap distribution now
hist(exp(gap), main = "Count ratio high/low after fix",
     xlab = "Multiplicative factor", breaks = 30)

# Prior checking ----------------------------------------------------------


# What does sigma_theta_stand = 2.5 mean for transition probabilities?
sigma_prior_draws <- abs(rnorm(10000, 2.5, 0.7))
stand_effects <- rnorm(10000, 0, 1) * sigma_prior_draws  # non-centered
grand_mean <- rnorm(10000, 1.2, 0.7)  # theta1 grand mean
theta_prior <- plogis(grand_mean + stand_effects)
hist(theta_prior, main = "Prior on theta1 across stands",
     xlab = "Transition probability")

# What does sigma ~ normal(0.5, 0.3) imply instead?
sigma_prior_draws <- abs(rnorm(10000, 0, 0.5))  # half-normal centered at 0
stand_effects <- rnorm(10000, 0, 1) * sigma_prior_draws
grand_mean <- rnorm(10000, 1.2, 0.7)
theta_prior <- plogis(grand_mean + stand_effects)
hist(theta_prior, main = "Revised prior on theta1 across stands",
     xlab = "Transition probability")

#New proposed prior on sigma_theta 1 and 2 of normal (0,0.5)


# Check what different priors imply on the count ratio scale
par(mfrow = c(1,3))

# Current prior
mu <- rnorm(10000, 1.5, 1.5)
hist(exp(log1p(exp(mu))), xlim = c(0,30), breaks = 50,
     main = "Current: normal(1.5, 1.5)", xlab = "Count ratio high/low")

# More constrained
mu <- rnorm(10000, 0.5, 0.7)
hist(exp(log1p(exp(mu))), xlim = c(0,30), breaks = 50,
     main = "Proposed: normal(0.5, 0.7)", xlab = "Count ratio high/low")

# Even tighter
mu <- rnorm(10000, 0.0, 0.5)
hist(exp(log1p(exp(mu))), xlim = c(0,30), breaks = 50,
     main = "Tight: normal(0.0, 0.5)", xlab = "Count ratio high/low")

# What does the grand mean gap imply alone?
mu_draws <- rnorm(10000, 1.5, 1.5)  # your current prior
hist(exp(log1p(exp(mu_draws))), 
     main = "Grand mean count ratio (high/low)",
     xlab = "Multiplicative factor")

# Species/stand effects now freely shrink or stretch this

# Divergence --------------------------------------------------------------

pairs(fit_all, pars = c("sigma_theta1_stand", "sigma_theta2_stand",
                        "grand_logit_theta1", "grand_logit_theta2",
                        "sigma_theta1_species", "sigma_theta2_species"))


# What does your current prior imply vs a more moderate one?
# Current: sigma_theta1_stand ~ normal(2.5, 0.7)
hist(plogis(rnorm(10000, 1.2, 0.7) + rnorm(10000, 0, abs(rnorm(10000, 2.5, 0.7)))),
     main = "Current prior", xlab = "theta1")

# Alternative: sigma_theta1_stand ~ normal(0, 1) — half-normal
hist(plogis(rnorm(10000, 1.2, 0.7) + rnorm(10000, 0, abs(rnorm(10000, 0, 1)))),
     main = "Alternative prior", xlab = "theta1")


fit_summary <- summary(fit_all)$summary
# Check which params still have NA Rhat
na_params <- rownames(fit_summary)[is.na(fit_summary[,"Rhat"])]
head(na_params, 10)

# Check if they are all state/log_omega (benign) or something new
table(gsub("\\[.*\\]", "", na_params))  # strips indices to show param names

# Only look at non-missing y_rep
fit_summary_yrep <- fit_summary[grep("y_rep", rownames(fit_summary)), ]
fit_summary_yrep_obs <- fit_summary_yrep[fit_summary_yrep[,"mean"] != 0, ]
summary(fit_summary_yrep_obs[, c("mean", "sd", "n_eff", "Rhat")])


# Are ALL log_omega stuck at zero (missing) or at real values?
fit_summary[grep("log_omega", rownames(fit_summary))[1:10], 
            c("mean", "sd", "n_eff", "Rhat")]

# How many log_omega are zero vs non-zero?
lo_means <- fit_summary[grep("log_omega", rownames(fit_summary)), "mean"]
table(lo_means == 0)

# Same for y_rep
yr_means <- fit_summary[grep("y_rep", rownames(fit_summary)), "mean"]
table(yr_means == 0)

table(stan_data_all$obs_missing)
# If most values are 1, the encoding is still wrong
mean(stan_data_all$obs_missing)  # should be a small fraction, not ~1.0


# Check current encoding
table(stan_data_all$obs_missing)

# If 1 = observed (wrong), flip it:
stan_data_all$obs_missing2 <- 1L - stan_data_all$obs_missing

# Verify: missing should be a small fraction
mean(stan_data_all$obs_missing2)  # expect something like 0.1–0.3, not 0.99


# What does obs_missing mean in your raw data?
# Check the 3 observations where obs_missing = 0
idx_observed <- which(stan_data_all$obs_missing == 0)
stan_data_all$y[idx_observed]        # what are their counts?
stan_data_all$sp[idx_observed]       # which species?

# And a sample of the 947 where obs_missing = 1
idx_missing <- which(stan_data_all$obs_missing == 1)
head(stan_data_all$y[idx_missing], 20)   # do these have real counts?


table(stan_data_all$obs_missing)


# State ID issue ----------------------------------------------------------
log_alpha_low_draws  <- rstan::extract(fit_all, "log_alpha_low")[[1]]
log_alpha_high_draws <- rstan::extract(fit_all, "log_alpha_high")[[1]]

# Mean gap between states per series
gap <- colMeans(log_alpha_high_draws) - colMeans(log_alpha_low_draws)
summary(gap)
par(mfrow= c(1,1))
hist(gap, main = "Log-scale gap between states", xlab = "log_alpha_high - log_alpha_low")

# What does this mean on the count scale?
summary(exp(gap))  # ratio of high to low state mean counts



draws <- as_draws_df(as.array(fit_all))
low_cols  <- grep("^log_alpha_low\\[", names(draws), value = TRUE)
high_cols <- grep("^log_alpha_high\\[", names(draws), value = TRUE)
ratios <- exp(draws[, high_cols] - draws[, low_cols])

apply(ratios, 2, median)

par(mfrow=c(1,1))
boxplot(ratios,
        las = 2,
        ylab = "High/Low mean ratio",
        main = "State separation")


draws <- as.data.frame(fit_all)

state_cols <- grep("^state\\[", names(draws), value = TRUE)

states <- draws[, state_cols]
state1_prop <- apply(states, 1, function(x) mean(x == 1))
state2_prop <- apply(states, 1, function(x) mean(x == 2))
summary(state1_prop)
summary(state2_prop)



theta1_cols <- grep("^theta1\\[", names(draws), value = TRUE)
theta2_cols <- grep("^theta2\\[", names(draws), value = TRUE)

theta1 <- draws[, theta1_cols]
theta2 <- draws[, theta2_cols]


diff_state <- diff(as.numeric(states[1, ]))
mean(diff_state != 0)


acf(as.numeric(states[1, ]))
rle(as.numeric(states[1, ]))$lengths
mean(rle(as.numeric(states[1, ]))$lengths)



state_matrix <- states  # draws x time OR summarized posterior

cor(state_series_i, state_series_j)
mean_state_t <- apply(states, 2, mean)
plot(mean_state_t, type="l")

mast_t <- colMeans(state_matrix == 2)
plot(mast_t, type="l")

cor(state_i == 2, state_j == 2)

ccf(state_i, state_j)



state_by_species <- t(sapply(1:S, function(s) {
  cols <- which(stan_data_all$sp == s)
  rowMeans(state_matrix[cols, , drop = FALSE] == 2)
}))

species_sync <- apply(state_by_species, 2, var)
species_sync2 <- apply(state_by_species, 2, function(x) {
  mean(abs(x - mean(x)))
})



# Victor's plot -----------------------------------------------------------

plot(stand_year_all$stand,stand_year_all$y, NA= TRUE)
# Order the stand factor by elevation
stand_year_all$stand <- factor(stand_year_all$stand, levels = stand_levels)

# Plot
plot(stand_year_all$stand, stand_year_all$y,
     xlab = "Stand (low → high elevation)",
     ylab = "Seed count",
     main = "Raw seed production per stand",
     las = 2)


# 1. Extract draws as matrix [n_draws x N_stands]
draws <- as.matrix(fit_all, pars = "alpha_low_stand")

# 2. Stand labels ordered by stand_id
stand_labels <- years_per_series %>%
  distinct(stand, stand_id) %>%
  arrange(stand_id) %>%
  pull(stand)

# 3. Compute quantiles [9 x N_stands]
quants <- apply(draws, 2, quantile, probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

# 4. Plot
N <- ncol(draws)
plot(NULL,
     xlim = c(0.5, N + 0.5), ylim = range(quants),
     xaxt = "n", xlab = "Stand", ylab = "alpha_low_stand (log scale)",
     main = "Shared stand effect on low-state seed production")
axis(1, at = 1:N, labels = stand_labels, las = 2)
abline(h = 0, lty = 2, col = "grey60")

for (k in 1:N) {
  lines(c(k, k), c(quants[1, k], quants[9, k]), lwd = 1, col = "steelblue")  # 80% CI
  lines(c(k, k), c(quants[3, k], quants[7, k]), lwd = 3, col = "steelblue")  # 60% CI
  points(k, quants[5, k], pch = 19, col = "navy")                             # median
}


# Series belonging to species 1, ordered by stand
meta_sp1 <- years_per_series %>%
  mutate(f = row_number()) %>%
  filter(spp == species_list[1]) %>%   # change index or use e.g. "ABAM"
  arrange(stand_id)

# Extract draws for those series only
draws_sp1 <- as.matrix(fit_all, pars = "log_alpha_low")[, meta_sp1$f]

# Quantiles
quants_sp1 <- apply(draws_sp1, 2, quantile,
                    probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

N <- ncol(draws_sp1)
plot(NULL,
     xlim = c(0.5, N + 0.5), ylim = range(quants_sp1),
     xaxt = "n", xlab = "Stand", ylab = "log_alpha_low (log scale)",
     main = paste("Implied stand effect —", species_list[1]))
axis(1, at = 1:N, labels = meta_sp1$stand, las = 2)
abline(h = 0, lty = 2, col = "grey60")

for (k in 1:N) {
  lines(c(k, k), c(quants_sp1[1, k], quants_sp1[9, k]), lwd = 1, col = "steelblue")
  lines(c(k, k), c(quants_sp1[3, k], quants_sp1[7, k]), lwd = 3, col = "steelblue")
  points(k, quants_sp1[5, k], pch = 19, col = "navy")
}


# Species list — in the same order as Stan's factor levels
species_list <- levels(as.factor(stand_year_all$spp))
# Double-check: should match what you printed during data prep
print(species_list)

# All series with species and stand info
meta_all <- years_per_series %>%
  mutate(f = row_number()) %>%
  arrange(stand_id)

# Extract all log_alpha_low draws [n_draws x F]
draws_all <- as.matrix(fit_all, pars = "log_alpha_low")

# Compute median per series
meta_all$median <- apply(draws_all[, meta_all$f], 2, median)
meta_all$q10    <- apply(draws_all[, meta_all$f], 2, quantile, probs = 0.1)
meta_all$q90    <- apply(draws_all[, meta_all$f], 2, quantile, probs = 0.9)

# Stand x-axis positions (consistent across species)
stand_levels <- years_per_series %>%
  distinct(stand, stand_id) %>%
  arrange(stand_id) %>%
  pull(stand)
N_stands <- length(stand_levels)

meta_all <- meta_all %>%
  mutate(x = match(stand, stand_levels))

# Colors per species
cols <- setNames(
  RColorBrewer::brewer.pal(length(species_list), "Dark2"),
  species_list
)

# Plot
par(mar = c(6, 4, 3, 12), xpd = TRUE)
plot(NULL,
     xlim = c(0.5, N_stands + 0.5),
     ylim = range(c(meta_all$q10, meta_all$q90)),
     xaxt = "n", xlab = "Stand", ylab = "log_alpha_low (log scale)",
     main = "Implied stand effects by species")
axis(1, at = 1:N_stands, labels = stand_levels, las = 2)
abline(h = 0, lty = 2, col = "grey70")
abline(v = 1:N_stands, col = "grey92", lwd = 1)

for (sp in species_list) {
  d <- meta_all %>% filter(spp == sp) %>% arrange(x)
  col <- cols[sp]
  
  # 80% CI as thin segment
  segments(d$x, d$q10, d$x, d$q90, col = adjustcolor(col, 0.4), lwd = 1)
  # Line connecting medians
  #lines(d$x, d$median, col = col, lwd = 1.5, lty = 2)
  # Median points
  points(d$x, d$median, col = col, pch = 19, cex = 1.1)
}

# Legend outside right margin
par(xpd = TRUE)
legend(x = N_stands + 1.5, y = mean(range(c(meta_all$q10, meta_all$q90))),
       legend = species_list, col = cols[species_list],
       lwd = 2, pch = 19, bty = "n", cex = 0.85,
       yjust = 0.5, title = "Species")



elevation_df <- tibble(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13","AO03","AG05",
            "AV06","AX15","AB08","PP17","AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600, 668, 700, 850, 1450, 850, 900, 950,
                1060, 1090, 1100, 1150, 1150, 1200, 1450, 1600, 1700, 1800)
)

stand_levels <- elevation_df %>%
  arrange(elevation) %>%
  pull(stand)

# Stand metadata with stand_id (used to match Stan's column order)
stand_meta <- years_per_series %>%
  distinct(stand, stand_id) %>%
  arrange(stand_id)   # this matches the column order in the draws matrix

# Reorder columns of draws to match elevation order
draws <- as.matrix(fit_all, pars = "alpha_low_stand")  # columns = stand_id order
col_order <- match(stand_levels, stand_meta$stand)     # elevation order → stand_id position
draws_ordered <- draws[, col_order]

# Quantiles on reordered draws
quants <- apply(draws_ordered, 2, quantile,
                probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

N <- ncol(draws_ordered)
plot(NULL,
     xlim = c(0.5, N + 0.5), ylim = range(quants),
     xaxt = "n", xlab = "Stand (low → high elevation)",
     ylab = "alpha_low_stand (log scale)",
     main = "Shared stand effect on low-state seed production")
axis(1, at = 1:N, labels = stand_levels, las = 2)
abline(h = 0, lty = 2, col = "grey60")

for (k in 1:N) {
  lines(c(k, k), c(quants[1, k], quants[9, k]), lwd = 1, col = "steelblue")
  lines(c(k, k), c(quants[3, k], quants[7, k]), lwd = 3, col = "steelblue")
  points(k, quants[5, k], pch = 19, col = "navy")
}


# Species list — in the same order as Stan's factor levels
species_list <- levels(as.factor(stand_year_all$spp))
# Double-check: should match what you printed during data prep
print(species_list)

# All series with species and stand info
meta_all <- years_per_series %>%
  mutate(f = row_number()) %>%
  arrange(stand_id)

# Extract all log_alpha_low draws [n_draws x F]
draws_all <- as.matrix(fit_all, pars = "log_alpha_low")

# Compute median per series
meta_all$median <- apply(draws_all[, meta_all$f], 2, median)
meta_all$q10    <- apply(draws_all[, meta_all$f], 2, quantile, probs = 0.1)
meta_all$q90    <- apply(draws_all[, meta_all$f], 2, quantile, probs = 0.9)

# Stand x-axis positions (consistent across species)
stand_levels <- elevation_df %>%
  arrange(elevation) %>%
  pull(stand)
N_stands <- length(stand_levels)

meta_all <- meta_all %>%
  mutate(x = match(stand, stand_levels))

# Colors per species
cols <- setNames(
  RColorBrewer::brewer.pal(length(species_list), "Dark2"),
  species_list
)

# Plot
par(mar = c(6, 4, 3, 12), xpd = TRUE)
plot(NULL,
     xlim = c(0.5, N_stands + 0.5),
     ylim = range(c(meta_all$q10, meta_all$q90)),
     xaxt = "n", xlab = "Stand", ylab = "log_alpha_low (log scale)",
     main = "Implied stand effects by species")
axis(1, at = 1:N_stands, labels = stand_levels, las = 2)
abline(h = 0, lty = 2, col = "grey70")
abline(v = 1:N_stands, col = "grey92", lwd = 1)

for (sp in species_list) {
  d <- meta_all %>% filter(spp == sp) %>% arrange(x)
  col <- cols[sp]
  
  # 80% CI as thin segment
  segments(d$x, d$q10, d$x, d$q90, col = adjustcolor(col, 0.4), lwd = 1)
  # Line connecting medians
  #lines(d$x, d$median, col = col, lwd = 1.5, lty = 2)
  # Median points
  points(d$x, d$median, col = col, pch = 19, cex = 1.1)
}

# Legend outside right margin
par(xpd = TRUE)
legend(x = N_stands + 1.5, y = mean(range(c(meta_all$q10, meta_all$q90))),
       legend = species_list, col = cols[species_list],
       lwd = 2, pch = 19, bty = "n", cex = 0.85,
       yjust = 0.5, title = "Species")


# Posterior median of y_rep per observation
draws_yrep <- as.matrix(fit_all, pars = "y_rep")  # [n_draws x N]
y_rep_median <- apply(draws_yrep, 2, median)

# Attach to data
stand_year_all$y_rep_median <- y_rep_median
stand_year_all$residual     <- stand_year_all$y - y_rep_median

stand_year_all$log_y      <- log(stand_year_all$y + 1)
stand_year_all$log_y_rep  <- log(y_rep_median + 1)
stand_year_all$resid_log  <- stand_year_all$log_y - stand_year_all$log_y_rep

resid_summary_log <- stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(mean_resid = mean(resid_log, na.rm = TRUE), .groups = "drop") %>%
  mutate(x = match(stand, stand_levels))

# Summarise residuals by species x stand
resid_summary <- stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(
    mean_resid = mean(resid_log, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(elevation_df, by = "stand") %>%
  mutate(stand = factor(stand, levels = stand_levels))  # elevation order

# Plot
cols <- setNames(
  RColorBrewer::brewer.pal(length(species_list), "Dark2"),
  species_list
)

par(mar = c(6, 4, 3, 12), xpd = FALSE)
plot(NULL,
     xlim = c(0.5, N_stands + 0.5),
     ylim = range(resid_summary$mean_resid, na.rm = TRUE),
     xaxt = "n", xlab = "Stand (low → high elevation)",
     ylab = "Mean residual (observed - predicted)",
     main = "Residuals by species and stand")
axis(1, at = 1:N_stands, labels = stand_levels, las = 2)
abline(h = 0, lty = 2, col = "grey60")
abline(v = 1:N_stands, col = "grey92", lwd = 1)

for (sp in species_list) {
  d <- resid_summary %>%
    filter(spp == sp) %>%
    mutate(x = match(stand, stand_levels))
  col <- cols[sp]
  points(d$x, d$mean_resid, col = col, pch = 19, cex = 1.1)
}

par(xpd = TRUE)
legend(x = N_stands + 1.5,
       y = mean(range(resid_summary$mean_resid, na.rm = TRUE)),
       legend = species_list, col = cols[species_list],
       pch = 19, bty = "n", cex = 0.85, yjust = 0.5, title = "Species")





# Get most probable state per observation
draws_state <- as.matrix(fit_all, pars = "state")
state_median <- apply(draws_state, 2, median)
stand_year_all$state <- round(state_median)

# Plot only low-state years
stand_year_low <- stand_year_all %>% filter(state == 1)
stand_year_low$stand <- factor(stand_year_low$stand, levels = stand_levels)

plot(stand_year_low$stand, stand_year_low$y,
     xlab = "Stand (low → high elevation)",
     ylab = "Seed count (low state years only)",
     main = "Raw seed production — low state years only",
     las = 2)

# diagnosis 04.06.2026 ----------------------------------------------------

# 1. What persistence probabilities did the model learn?
print(fit_all, pars = c("grand_logit_theta1", "grand_logit_theta2",
                        "sigma_theta1_species", "sigma_theta2_species"))

# 2. What is the state sequence for a problematic species x stand?
# e.g. TSHE at SUNR
meta_problem <- years_per_series %>%
  mutate(f = row_number()) %>%
  filter(spp == "TSHE", stand == "TO04")

states_problem <- stand_year_all %>%
  filter(spp == "TSHE", stand == "TO04") %>%
  mutate(state = state_median[row_number()])  

print(states_problem %>% select(year, y, state))

# 3. Check the delta (gap between states) per species
print(fit_all, pars = c("mu_log_delta", 
                        "sigma_log_delta_species",
                        "log_delta_species_nc"))


tshecheck<- stand_year_all %>%
  filter(spp == "TSHE", stand == "TO04") %>%
  arrange(y) %>%
  print(n = Inf)

plot(tshecheck$year,tshecheck$y)


# Extract prior predictive states
state_prior <- as.matrix(fit_prior, pars = "state")
state_prior_median <- apply(state_prior, 2, median)
stand_year_all$state_prior <- round(state_prior_median)

# Check switching for TSHE at TO04
stand_year_all %>%
  filter(spp == "TSHE", stand == "TO04") %>%
  select(year, y, state_prior)

# What fraction of years are assigned state 2 (mast) overall?
mean(stand_year_all$state_prior == 2)

# Check transitions: does the prior produce switching?
table(head(stand_year_all$state_prior, -1), 
      tail(stand_year_all$state_prior, -1))


print(fit_all, pars = c("alpha_low_species", "mu_log_low", "sigma_low_species"))

print(fit_TSHE, pars = c("grand_logit_theta1", "grand_logit_theta2"))
print(fit_TSME, pars = c("grand_logit_theta1", "grand_logit_theta2"))
print(fit_THPL, pars = c("grand_logit_theta1", "grand_logit_theta2"))
print(fit_ABAM, pars = c("grand_logit_theta1", "grand_logit_theta2"))
print(fit_ABLA, pars = c("grand_logit_theta1", "grand_logit_theta2"))
print(fit_CANO, pars = c("grand_logit_theta1", "grand_logit_theta2"))
print(fit_PSME, pars = c("grand_logit_theta1", "grand_logit_theta2"))

print(fit_all, pars = c("grand_logit_theta1", "grand_logit_theta2", "theta1[1]","theta1[2]","theta1[3]","theta1[4]","theta1[5]","theta1[6]","theta1[7]","theta2[1]","theta2[2]","theta2[3]","theta2[4]","theta2[5]","theta2[6]","theta2[7]"))


# Recompute state median from new fit
state_new <- as.matrix(fit_all, pars = "state")
state_new_median <- apply(state_new, 2, median)
stand_year_all$state_new <- round(state_new_median)

# Check TSHE at TO04
stand_year_all %>%
  filter(spp == "TSHE", stand == "TO04") %>%
  select(year, y, state_new)

# Get log_alpha_low and log_alpha_high for TSHE series
meta_tshe <- years_per_series %>%
  mutate(f = row_number()) %>%
  filter(spp == "TSHE")

draws_low  <- as.matrix(fit_all, pars = "log_alpha_low")
draws_high <- as.matrix(fit_all, pars = "log_alpha_high")

# Posterior median of implied means for TSHE
tshe_low  <- apply(draws_low[,  meta_tshe$f], 2, median)
tshe_high <- apply(draws_high[, meta_tshe$f], 2, median)

data.frame(
  stand     = meta_tshe$stand,
  low_mean  = round(exp(tshe_low)),
  high_mean = round(exp(tshe_high)),
  ratio     = round(exp(tshe_high - tshe_low), 1)
)

# What is mu_log_low + alpha_low_species[6] for TSHE?
mu_draws    <- as.matrix(fit_all, pars = "mu_log_low")[,1]
alpha_draws <- as.matrix(fit_all, pars = "alpha_low_species")[, 6]  # TSHE = species 6

tshe_baseline <- median(exp(mu_draws + alpha_draws))
cat("TSHE implied baseline (no stand effect):", round(tshe_baseline), "seeds\n")

# Compare to actual TSHE low years
stand_year_all %>%
  filter(spp == "TSHE") %>%
  summarise(
    min_y    = min(y),
    q10_y    = quantile(y, 0.1),
    q25_y    = quantile(y, 0.25),
    median_y = median(y)
  )


# Check state assignment for ALL TSHE stands
stand_year_all %>%
  filter(spp == "TSHE") %>%
  group_by(stand) %>%
  summarise(
    n_low  = sum(state_new == 1),
    n_high = sum(state_new == 2),
    min_y  = min(y),
    max_y  = max(y),
    .groups = "drop"
  ) %>%
  arrange(desc(n_low))


#works better when removing the partial pooling of species on the low 

# Which species are working correctly now (without alpha_low_species)?
stand_year_all %>%
  group_by(spp) %>%
  summarise(
    n_low  = sum(state_new == 1),
    n_high = sum(state_new == 2),
    pct_high = round(mean(state_new == 2) * 100),
    .groups = "drop"
  )

stand_year_all %>%
  filter(spp == "THPL") %>%
  group_by(stand) %>%
  summarise(
    n_low  = sum(state_new == 1),
    n_high = sum(state_new == 2),
    min_y  = min(y),
    max_y  = max(y),
    .groups = "drop"
  )


# Fit a version WITH alpha_low_species but print the implied
# low state means per species to see what's happening
print(fit_all, pars = c("mu_log_low",
                                     "alpha_low_species",
                                     "sigma_low_species",
                                     "mu_log_delta",
                                     "sigma_log_delta_species"))

print(fit_multi, pars = c("mu_log_low",
                        "alpha_low_species",
                        "sigma_low_species",
                        "mu_log_delta",
                        "sigma_log_delta_species"))

print(
  fit_multi,
  pars = c(
    "log_phi_state",
    "sigma_theta1_species",
    "sigma_theta2_species",
    "sigma_theta1_stand",
    "sigma_theta2_stand"
  )
)


print(fit_multi, pars = c("mu_log_delta", 
                        "sigma_log_delta_species",
                        "grand_logit_theta1",
                        "grand_logit_theta2"))

# And state assignments per species
stand_year_all %>%
  group_by(spp) %>%
  summarise(
    n_low    = sum(state_new == 1),
    n_high   = sum(state_new == 2),
    pct_high = round(mean(state_new == 2) * 100),
    .groups  = "drop"
  )


print(fit_multi, pars = c("logit_theta1_species",
                          "logit_theta2_species",
                          "mu_log_delta",
                          "sigma_log_delta_species"))

stand_year_all %>%
  group_by(spp) %>%
  summarise(
    n_low    = sum(state_new == 1),
    n_high   = sum(state_new == 2),
    pct_high = round(mean(state_new == 2) * 100),
    .groups  = "drop"
  )
