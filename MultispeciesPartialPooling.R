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


# 
# The theta2 result is actually scientifically interesting — it suggests that once a species enters the high masting state, the probability of staying there doesn't differ much between species. That's a real ecological finding, not just a computational nuisance.

##Next things to do : 

#1 : Add the 2009 year, handle the SUNR missing year
#2 : Trap-level
#3 : Extract values for the synchrony.
#4 : Clean my code


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


# Fitting Model -----------------------------------------------------------

fit_all <- stan(
    file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_delta.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)



# Diagnostic plots --------------------------------------------------------

check_hmc_diagnostics(fit_all)

#Checking the divergence issues
# Extract posterior draws (excluding warmup)
draws <- as_draws_df(fit_all)

# Extract sampler params from rstan
sampler_params <- get_sampler_params(fit_all, inc_warmup = FALSE)

# Combine across chains into one vector
divergent <- do.call(rbind, sampler_params)[, "divergent__"]

#Stands in HIGH
# Extract all alpha_high_stand columns
alpha_cols_high <- grep("^alpha_high_stand_nc\\[", names(draws), value = TRUE)

plot_df_high <- draws[, alpha_cols_high] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_high_stand = log(draws$sigma_high_stand),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_high),
    names_to  = "stand",
    values_to = "alpha_high_stand_nc"
  )
ggplot(plot_df_high, aes(x = alpha_high_stand_nc, y = log_sigma_high_stand)) +
  geom_point(data = subset(plot_df_high, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_high, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ stand) +
  labs(
    x     = "alpha_high_stand_nc",
    y     = "log(sigma_high_stand)",
    title = "Divergence diagnostic: high state stands"
  ) +
  theme_bw()

#Stands in LOW
alpha_cols_low <- grep("^alpha_low_stand_nc\\[", names(draws), value = TRUE)
plot_df_low <- draws[, alpha_cols_low] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_low_stand = log(draws$sigma_low_stand),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_low),
    names_to  = "stand",
    values_to = "alpha_low_stand_nc"
  )
ggplot(plot_df_low, aes(x = alpha_low_stand_nc, y = log_sigma_low_stand)) +
  geom_point(data = subset(plot_df_low, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_low, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ stand) +
  labs(
    x     = "alpha_low_stand_nc",
    y     = "log(sigma_low_stand)",
    title = "Divergence diagnostic: low state stands"
  ) +
  theme_bw()


#Species in LOW
# Extract all alpha_low_species columns
alpha_cols_low_sp <- grep("^alpha_low_species_nc\\[", names(draws), value = TRUE)
plot_df_low_sp <- draws[, alpha_cols_low_sp] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_low_species = log(draws$sigma_low_species),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_low_sp),
    names_to  = "species",
    values_to = "alpha_low_species_nc"
  )
ggplot(plot_df_low_sp, aes(x = alpha_low_species_nc, y = log_sigma_low_species)) +
  geom_point(data = subset(plot_df_low_sp, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_low_sp, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ species) +
  labs(
    x     = "alpha_low_species_nc",
    y     = "log(sigma_low_species)",
    title = "Divergence diagnostic: low state species"
  ) +
  theme_bw()

#For all the species in HIGH
# Extract all alpha_low_species columns
alpha_cols_high_sp <- grep("^alpha_high_species_nc\\[", names(draws), value = TRUE)
plot_df_high_sp <- draws[, alpha_cols_high_sp] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_high_species = log(draws$sigma_high_species),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_high_sp),
    names_to  = "species",
    values_to = "alpha_high_species_nc"
  )
ggplot(plot_df_high_sp, aes(x = alpha_high_species_nc, y = log_sigma_high_species)) +
  geom_point(data = subset(plot_df_high_sp, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_high_sp, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ species) +
  labs(
    x     = "alpha_high_species_nc",
    y     = "log(sigma_high_species)",
    title = "Divergence diagnostic: high state species"
  ) +
  theme_bw()


#For theta 1
alpha_cols_high_theta1_sp <- grep("^alpha_theta1_species_nc\\[", names(draws), value = TRUE)
plot_df_theta1_sp <- draws[, alpha_cols_high_theta1_sp] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_theta1_species = log(draws$sigma_theta1_species),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_high_theta1_sp),
    names_to  = "species",
    values_to = "alpha_theta1_species_nc"
  )
ggplot(plot_df_theta1_sp, aes(x = alpha_theta1_species_nc, y = log_sigma_theta1_species)) +
  geom_point(data = subset(plot_df_theta1_sp, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_theta1_sp, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ species) +
  labs(
    x     = "alpha_theta1_species_nc",
    y     = "log(sigma_theta1_species)",
    title = "Divergence diagnostic: theta1 species"
  ) +
  theme_bw()

#Theta 1 stand
alpha_cols_high_theta1_st <- grep("^alpha_theta1_stand_nc\\[", names(draws), value = TRUE)
plot_df_theta1_st <- draws[, alpha_cols_high_theta1_st] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_theta1_stand = log(draws$sigma_theta1_stand),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_high_theta1_st),
    names_to  = "stand",
    values_to = "alpha_theta1_stand_nc"
  )
ggplot(plot_df_theta1_st, aes(x = alpha_theta1_stand_nc, y = log_sigma_theta1_stand)) +
  geom_point(data = subset(plot_df_theta1_st, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_theta1_st, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ stand) +
  labs(
    x     = "alpha_theta1_stand_nc",
    y     = "log(sigma_theta1_stand)",
    title = "Divergence diagnostic: theta1 stand"
  ) +
  theme_bw()

#For theta 2
alpha_cols_theta2_sp <- grep("^alpha_theta2_species_nc\\[", names(draws), value = TRUE)
plot_df_theta2_sp <- draws[, alpha_cols_theta2_sp] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_theta2_species = log(draws$sigma_theta2_species),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_theta2_sp),
    names_to  = "species",
    values_to = "alpha_theta2_species_nc"
  )
ggplot(plot_df_theta2_sp, aes(x = alpha_theta2_species_nc, y = log_sigma_theta2_species)) +
  geom_point(data = subset(plot_df_theta2_sp, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_theta2_sp, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ species) +
  labs(
    x     = "alpha_theta2_species_nc",
    y     = "log(sigma_theta2_species)",
    title = "Divergence diagnostic: theta2 species"
  ) +
  theme_bw()

#Theta 2 stand
alpha_cols_high_theta2_st <- grep("^alpha_theta2_stand_nc\\[", names(draws), value = TRUE)
plot_df_theta2_st <- draws[, alpha_cols_high_theta2_st] %>%
  as.data.frame() %>%
  mutate(
    log_sigma_theta2_stand = log(draws$sigma_theta2_stand),
    divergent = factor(divergent, levels = c(0, 1), labels = c("No", "Yes"))
  ) %>%
  pivot_longer(
    cols      = all_of(alpha_cols_high_theta2_st),
    names_to  = "stand",
    values_to = "alpha_theta2_stand_nc"
  )
ggplot(plot_df_theta2_st, aes(x = alpha_theta2_stand_nc, y = log_sigma_theta2_stand)) +
  geom_point(data = subset(plot_df_theta2_st, divergent == "No"),
             colour = "grey60", alpha = 0.3, size = 0.5) +
  geom_point(data = subset(plot_df_theta2_st, divergent == "Yes"),
             colour = "red", alpha = 0.8, size = 1.5) +
  facet_wrap(~ stand) +
  labs(
    x     = "alpha_theta2_stand_nc",
    y     = "log(sigma_theta2_stand)",
    title = "Divergence diagnostic: theta2 stand"
  ) +
  theme_bw()


#Checking the overdispersion geometry OLD
plot_df <- data.frame(
  Low  = draws$phi1,
  High = draws$phi2
) %>%
  pivot_longer(everything(), names_to = "state", values_to = "phi")

ggplot(plot_df, aes(x = phi, fill = state, colour = state)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values   = c("Low" = "steelblue", "High" = "orange")) +
  scale_colour_manual(values = c("Low" = "steelblue", "High" = "orange")) +
  labs(
    x      = "phi (dispersion)",
    y      = "Density",
    title  = "Posterior dispersion: low vs high state",
    fill   = "State",
    colour = "State"
  ) +
  theme_bw()

ggplot(as.data.frame(draws), aes(x = phi1, y = phi2)) +
  geom_point(alpha = 0.2, size = 0.5, colour = "grey40") +
  geom_abline(linetype = "dashed", colour = "red") +
  labs(
    x     = "phi1 (low state)",
    y     = "phi2 (high state)",
    title = "Joint posterior of phi1 and phi2"
  ) +
  theme_bw()

# Diagnostics -------------------------------------------------------------

util<- new.env()
#creating a source file to make all the plots. 
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit_all)
util$check_all_hmc_diagnostics(diagnostics)


# extraire les posterior values
samples <- util$extract_expectand_vals(fit_all)

# diagnostics parametre par parametre
base_samples <- util$filter_expectands(samples,
                                       c("rho", "grand_logit_theta1","grand_logit_theta2",
                                         "alpha_theta1_species_nc", "alpha_theta2_species_nc",
                                         "sigma_theta1_species", "sigma_theta2_species","alpha_theta1_stand_nc", "sigma_theta1_stand", "alpha_theta2_stand_nc","sigma_theta1_stand", "sigma_theta2_stand", "alpha_low_species","alpha_high_species","sigma_low_species", "sigma_high_species","alpha_low_stand","alpha_high_stand_nc", "sigma_low_stand","sigma_high_stand","phi1","phi2"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)


base_samples <- util$filter_expectands(samples,
                                       c("log_delta_high_stand", "sigma_log_delta_high_stand"), check_arrays = TRUE)


util$plot_div_pairs("log_delta_high_stand[1]", "sigma_log_delta_high_stand", samples,
                    diagnostics, transforms = list("sigma_log_delta_high_stand" =1))

util$plot_div_pairs("log_delta_high_species[1]", "sigma_log_delta_high_species", samples,
                    diagnostics)




util$plot_pairs_by_chain(samples[["sigma_high_stand"]], "sigma_high_stand",
                         samples[["sigma_high_species"]], "sigma_high_species")

util$plot_expectand_pushforward(samples[["sigma_high_stand"]], 30)

util$plot_pairs_by_chain(samples[["phi1"]], "phi1",
                         samples[["phi2"]], "phi2")

util$plot_pairs_by_chain(samples[["log_alpha_low[10]"]], "log_alpha_low[10]",
                         samples[["log_alpha_high[10]"]], "log_alpha_high[10]")


util$plot_pairs_by_chain(samples[["log_alpha_low[45]"]], "log_alpha_low[45]",
                         samples[["log_alpha_high[45]"]], "log_alpha_high[45]")


par(mfrow = c(2,1))
f <- 40
start_id <- start_idxs[f];
end_id <- end_idxs[f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"),
                                     display_ylim = c(0,400))
points(x = 1:length(start_id:end_id), y = stan_data_all$y[start_id:end_id], pch = 20)


# PPC plots ---------------------------------------------------------------
# Extract y_rep and observed y
y_obs  <- stan_data_all$y
y_rep  <- as.matrix(fit_all, pars = "y_rep")

# ── 1. Global distribution ────────────────────────────────────────────────────
ppc_dens_overlay(y_obs, y_rep[1:200, ]) +
  labs(title = "PPC: overall distribution") +
  xlim(0,500)+
  theme_bw()

# log scale version — often more readable for seed count data
ppc_dens_overlay(log1p(y_obs), log1p(y_rep[1:200, ])) +
  labs(title = "PPC: overall distribution (log1p scale)") +
  theme_bw()

# ── 2. Summary statistics ─────────────────────────────────────────────────────
# Mean
ppc_stat(y_obs, y_rep, stat = "mean") +
  labs(title = "PPC: mean") + theme_bw()

# Standard deviation
ppc_stat(y_obs, y_rep, stat = "sd") +
  labs(title = "PPC: sd") + theme_bw()

# Proportion of zeros
prop_zero <- function(x) mean(x == 0)
ppc_stat(y_obs, y_rep, stat = prop_zero) +
  labs(title = "PPC: proportion of zeros") + theme_bw()
# Bayesian p-value for proportion of zeros
# Should be between 0.05 and 0.95 to be acceptable
prop_zero <- function(x) mean(x == 0)
p_value <- mean(apply(y_rep, 1, prop_zero) <= prop_zero(y_obs))
cat("Bayesian p-value for prop zeros:", round(p_value, 3), "\n")

# Max — checks for extreme values
ppc_stat(y_obs, y_rep, stat = "max") +
  labs(title = "PPC: max") + theme_bw()

# ── 3. By species ─────────────────────────────────────────────────────────────
species_id <- stan_data_all$sp
species_names <- levels(as.factor(stand_year_all$spp))

ppc_dens_overlay_grouped(
  y        = log1p(y_obs),
  yrep     = log1p(y_rep[1:200, ]),
  group    = species_names[species_id]
) +
  labs(title = "PPC by species (log1p scale)") +
  theme_bw()

ppc_stat_grouped(
  y      = y_obs,
  yrep   = y_rep,
  group  = species_names[species_id],
  stat   = "mean"
) +
  labs(title = "PPC: mean by species") +
  theme_bw()

ppc_stat_grouped(
  y      = y_obs,
  yrep   = y_rep,
  group  = species_names[species_id],
  stat   = prop_zero
) +
  labs(title = "PPC: proportion of zeros by species") +
  theme_bw()

# ── 4. By stand ───────────────────────────────────────────────────────────────
# Build stand label per observation (length N)
stand_labels <- years_per_series$stand[
  findInterval(seq_along(y_obs), stan_data_all$start_idxs)
]

ppc_stat_grouped(
  y      = y_obs,
  yrep   = y_rep,
  group  = stand_labels,
  stat   = "mean"
) +
  labs(title = "PPC: mean by stand") +
  theme_bw()

ppc_stat_grouped(
  y      = y_obs,
  yrep   = y_rep,
  group  = stand_labels,
  stat   = prop_zero
) +
  labs(title = "PPC: proportion of zeros by stand") +
  theme_bw()

stand_id_obs <- rep(NA_integer_, stan_data_all$N)

for (f in seq_len(stan_data_all$F)) {
  idx <- stan_data_all$start_idxs[f]:stan_data_all$end_idxs[f]
  stand_id_obs[idx] <- stan_data_all$stand_id[f]
}
stand_labels <- paste0("Stand_", stand_id_obs)
y_obs <- stan_data_all$y  # or your y_obs object

ppc_dens_overlay_grouped(
  y     = log1p(y_obs),
  yrep  = log1p(y_rep[1:200, ]),
  group = stand_labels
) +
  labs(title = "PPC by stand (log1p scale)") +
  theme_bw()


# Prior vs posterior ------------------------------------------------------


n_prior <- 1e4

half_norm <- function(n, sd) abs(rnorm(n, 0, sd))


# 2. SAMPLE FROM PRIORS

prior_list <- list()

# Grand means (emissions)
prior_list$grand_mean_low  <- rnorm(n_prior, 2.6, 1.0)
prior_list$grand_mean_high <- rnorm(n_prior, 5.7, 1.0)

# Transition logits
prior_list$grand_logit_theta1 <- rnorm(n_prior, 1.4, 1)
prior_list$grand_logit_theta2 <- rnorm(n_prior, -1.4, 1)

# Variance parameters (half-normal)
prior_list$sigma_theta1_species <- half_norm(n_prior, 0.5)
prior_list$sigma_theta2_species <- half_norm(n_prior, 0.5)
prior_list$sigma_theta1_stand   <- half_norm(n_prior, 0.5)
prior_list$sigma_theta2_stand   <- half_norm(n_prior, 0.5)

prior_list$sigma_low_species  <- half_norm(n_prior, 0.5)
prior_list$sigma_high_species <- half_norm(n_prior, 0.5)
prior_list$sigma_low_stand    <- half_norm(n_prior, 0.5)
prior_list$sigma_high_stand   <- half_norm(n_prior, 0.5)

# Dispersion
prior_list$phi1 <- rgamma(n_prior, 4.0, 0.6)
prior_list$phi2 <- rgamma(n_prior, 4.0, 0.6)

# Convert to data frame
prior_df <- as.data.frame(prior_list)

# Transform logits → probabilities (interpretability)
prior_df$theta1 <- plogis(prior_df$grand_logit_theta1)
prior_df$theta2 <- plogis(prior_df$grand_logit_theta2)


# 3. EXTRACT POSTERIORS

post <- rstan::extract(fit_all)

post_list <- list()

# Grand means
post_list$grand_mean_low  <- post$grand_mean[, 1]
post_list$grand_mean_high <- post$grand_mean[, 2]

# Transition logits
post_list$grand_logit_theta1 <- post$grand_logit_theta1
post_list$grand_logit_theta2 <- post$grand_logit_theta2

# Variance parameters
post_list$sigma_theta1_species <- post$sigma_theta1_species
post_list$sigma_theta2_species <- post$sigma_theta2_species
post_list$sigma_theta1_stand   <- post$sigma_theta1_stand
post_list$sigma_theta2_stand   <- post$sigma_theta2_stand

post_list$sigma_low_species  <- post$sigma_low_species
post_list$sigma_high_species <- post$sigma_high_species
post_list$sigma_low_stand    <- post$sigma_low_stand
post_list$sigma_high_stand   <- post$sigma_high_stand

# Dispersion
post_list$phi1 <- post$phi1
post_list$phi2 <- post$phi2

posterior_df <- as.data.frame(post_list)

# Transform logits → probabilities
posterior_df$theta1 <- plogis(posterior_df$grand_logit_theta1)
posterior_df$theta2 <- plogis(posterior_df$grand_logit_theta2)


# 4. RESHAPE TO LONG FORMAT

to_long <- function(df, label) {
  stack_df <- utils::stack(df)
  names(stack_df) <- c("value", "parameter")
  stack_df$distribution <- label
  stack_df
}

prior_long     <- to_long(prior_df, "Prior")
posterior_long <- to_long(posterior_df, "Posterior")

combined <- rbind(prior_long, posterior_long)


ggplot(combined, aes(x = value, colour = distribution)) +
  geom_density(linewidth = 0.7) +
  facet_wrap(~parameter, scales = "free") +
  scale_colour_manual(values = c("Prior" = "steelblue", "Posterior" = "tomato")) +
  labs(
    title  = "Prior vs Posterior — Hierarchical HMM",
    x      = "Parameter value",
    y      = "Density",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text      = element_text(size = 8),
    legend.position = "top"
  )



# New ---------------------------------------------------------------------

# Extract grand_mean posteriors
grand_mean_draws <- draws %>%
  select(`grand_mean[1]`, `grand_mean[2]`) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value")

# Get summary stats to paste into the widget
grand_mean_draws %>%
  group_by(parameter) %>%
  summarise(
    post_mean = mean(value),
    post_sd   = sd(value),
    post_q025 = quantile(value, 0.025),
    post_q975 = quantile(value, 0.975)
  )


# 1. Posterior means of theta1 and theta2
mcmc_intervals(fit_all, pars = c("grand_logit_theta1", "grand_logit_theta2"))

# 2. State separation per series
mcmc_intervals(fit_all, regex_pars = "log_alpha_high|log_alpha_low")

# 3. Pairs plot for the ridge
mcmc_pairs(fit_all, pars = c("grand_logit_theta2", "log_delta_high_grand_mean"))
