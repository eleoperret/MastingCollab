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


# Data import -------------------------------------------------------------


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

fit_all3 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesGitIssue74.stan",
  data    = stan_data_all,
  iter    = 1000, #change based on how much iterations you need
  warmup  = 500, #make the warmup longer 
  chains  = 1,
  seed    = 123,
)

fit_all2 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesGitIssue75.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

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

#Problematic things example

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


samples <- rstan::extract(fit_all2)

n <- length(samples$mu_log_low)

priors <- data.frame(
  mu_log_low           = rnorm(n, 2, 1.5),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 0.5, 1),
  grand_logit_theta2   = rnorm(n, -1, 0.5),
  sigma_low_species    = abs(rnorm(n, 0.5, 1)),
  sigma_stand          = abs(rnorm(n,1,0.5)),
  #sigma_low_stand      = abs(rnorm(n, 0.5, 1)),
  sigma_log_delta_sp   = abs(rnorm(n, 1.5, 0.7)),
  #sigma_log_delta_stand= abs(rnorm(n, 0, 1)),
  sigma_theta1_species = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_species = abs(rnorm(n, 0, 0.7)),
  sigma_theta1_stand   = abs(rnorm(n, 1, 0.3)),
  sigma_theta2_stand   = abs(rnorm(n, 1, 0.3)),
  phi_species          = exp(rnorm(n, log(4), 0.6)),
  phi_low              = exp(rnorm(n,0,0.3)),
  phi_high             = exp(rnorm(n,0,0.3))
)

posteriors <- data.frame(
  mu_log_low           = samples$mu_log_low,
  mu_log_delta         = samples$mu_log_delta,
  grand_logit_theta1   = samples$grand_logit_theta1,
  grand_logit_theta2   = samples$grand_logit_theta2,
  sigma_low_species    = samples$sigma_low_species,
  sigma_stand          = samples$sigma_stand,
  sigma_log_delta_sp   = samples$sigma_log_delta_species,
  sigma_theta1_species = samples$sigma_theta1_species,
  sigma_theta2_species = samples$sigma_theta2_species,
  sigma_theta1_stand   = samples$sigma_theta1_stand,
  sigma_theta2_stand   = samples$sigma_theta2_stand,
  phi_species          = exp(rowMeans(samples$log_phi_species)),
  phi_low              = exp(samples$log_phi_state[, 1]),
  phi_high             = exp(samples$log_phi_state[, 2])
)

# --- Nice labels ---
param_labels <- c(
  mu_log_low            = "grand mean seeds, low state",
  mu_log_delta          = "grand mean fold-change",
  grand_logit_theta1    = "P stay low",
  grand_logit_theta2    = "P stay high",
  sigma_low_species     = "σ low species",
  sigma_stand       = "σ low stand",
  sigma_log_delta_sp    = "σ log delta species",
  sigma_theta1_species  = "σ θ₁ species",
  sigma_theta2_species  = "σ θ₂ species",
  sigma_theta1_stand    = "σ θ₁ stand",
  sigma_theta2_stand    = "σ θ₂ stand",
  phi_species           = "φ species",   
  phi_low               = "φ low state ",
  phi_high              = "φ high state "
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




# Posteriror vs prior per SPECIES ---------------------------------

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



# Diagnostic plots Centered --------------------------------------------------------


#If there is a divergence : 
#where is this divergence?
divergent <- get_sampler_params(fit_all, inc_warmup = FALSE) |>
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

pairs(fit_all, pars=c("alpha_low_species[7]","sigma_low_species","alpha_stand[7,1]","sigma_stand"))
pairs(fit_all,pars=c("alpha_stand"))

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


# Diagnostics code line -------------------------------------------------------------


# Checking model samples --------------------------------------------------
diagnostics <- util$extract_hmc_diagnostics(fit_all)
util$check_all_hmc_diagnostics(diagnostics)


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


# Checking chains ---------------------------------------------------------
samples   <- util$extract_expectand_vals(fit_all)

#Parameters by chain
util$plot_pairs_by_chain(samples[["log_alpha_low[20]"]], "log_alpha_low[20]",
                         samples[["log_alpha_high[20]"]], "log_alpha_high[20]")


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


# To find the most bimodal series
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

# Checking parameters  ----------------------------------------

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

#Parameter : Seed mean production (Difference between the two states) 
delta_vals <- samples[["log_alpha_high[6]"]] - samples[["log_alpha_low[6]"]]

util$plot_expectand_pushforward(
  delta_vals, 30,
  display_name = "log_alpha_high[6] - log_alpha_low[6]",
  main = "implied log gap for f=6"
)










# Rhat is Na --------------------------------------------------------------

#This means that the sampler is stuck with similar values for a or more parameters and that the variance cannot be defined. 

# Find parameters with NA Rhat
na_rhat_params <- rownames(fit_summary)[is.na(fit_summary[, "Rhat"])]
na_rhat_params

#Checking what they are
stan_data_all$y[c(183, 289, 432)] 
fit_summary[is.na(fit_summary[, "Rhat"]), ]

# Find parameters with low ESS 
fit_summary_clean <- fit_summary[!is.na(fit_summary[,"Rhat"]), ]
# Lowest ESS continuous params
head(fit_summary_clean[order(fit_summary_clean[,"n_eff"]), 
                       c("mean","sd","n_eff","Rhat")], 20)

# Parameters with poor mixing
fit_summary_clean[fit_summary_clean[,"Rhat"] > 1.05, 
                  c("mean","sd","n_eff","Rhat")]



# Extract posterior draws for a suspicious parameter
draws <- rstan::extract(fit_all)

# Check variance of specific parameters - if 0, it's stuck
sapply(draws, function(x) var(as.vector(x)))

# Check if the parameter hits 0, 1, or extreme boundary
for (p in na_rhat_params) {
  cat(p, ":", range(as.vector(rstan::extract(fit_all, p)[[1]])), "\n")
}




# Check n_eff alongside Rhat - both NA = never sampled properly
fit_summary[fit_summary[,"n_eff"] < 10, ]





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
# Extract state draws
state_draws <- rstan::extract(fit_all, "state")[[1]]  # shape: [samples x N]

# For the varying sites, what's the mean state?
par(mfrow= c(1,1))
hist(apply(state_draws, 2, mean), 
     main = "Mean state across sites", xlab = "Mean state (1=absent, 2=present)")


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

# And state assignments per species
stand_year_all %>%
  group_by(spp) %>%
  summarise(
    n_low    = sum(state_new == 1),
    n_high   = sum(state_new == 2),
    pct_high = round(mean(state_new == 2) * 100),
    .groups  = "drop",
    min_y    = min(y),
    q10_y    = quantile(y, 0.1),
    q25_y    = quantile(y, 0.25),
    median_y = median(y)
  )





