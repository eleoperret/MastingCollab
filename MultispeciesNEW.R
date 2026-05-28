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
  file    = "Stan_code/Species_Stan_Model/MultispeciesNEWMikeNewPriorsPhiSpecies.stan",
  data    = stan_data_all,
  iter    = 4000, #change based on how much iterations you need
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

divergent <- get_sampler_params(fit_all2, inc_warmup = FALSE)
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
