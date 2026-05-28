##Code for multispecies HMM model 
# Stan code: 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

####Doesn't include PARA and SUNR as I do not have a map of the species. 
####Issues with ABLA so I will remove it. 


##Partially pooling also species
##theta partially pooled by species and stand
##make my priors wider. 
##check if the stands is all stands or just the selected stands. 
##high mu for the higher state
##add the new parameters as centered
##Check the centered (on the psme data). 


##Next things to do : 

#1 : Add the 2009 year, ABLA and handle the SUNR missing year
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
  filter(spp %in% c("ABAM","CANO","PSME","TSHE","THPL"))
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

# # Overlook at the data
# stand_year_all2<- seed_filtered %>%
#   group_by(spp, stand) %>%
#   summarise(
#     max_seed  = max(total_viable_sds),
#     mean_seed = mean(total_viable_sds),
#     zeros  = mean(total_viable_sds == 0),
#     n_obs  = n(),
#     .groups = "drop"
#   )
# write.csv(stand_year_all, file = "C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting/OverlookOfData.csv")

#Not sure if I should do that yet or not. 
# seed_filtered_clean <- seed_filtered %>%
#   filter(!(spp == "CANO" & stand %in% c("AB08", "AV06", "AV14", "AX15", "PP17")))
# seed_filtered_clean <- seed_filtered %>%
#   filter(!(spp == "ABAM" & stand %in% c("AB08", "TO04", "TB13", "AX15", "PP17", "TA01")))


# #Total seeds per year
# total_stand<-seed_filtered%>%
#   group_by(spp,year)%>%
#   summarise(total_seeds=sum(total_viable_sds))

#which stand share which species?
stand_year_all %>%
  group_by(stand) %>%
  summarise(
    n_species = n_distinct(spp),
    species   = paste(sort(unique(spp)), collapse = ", ")
  ) %>%
  arrange(desc(n_species))

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
  summarise(T_i = n(), .groups = "drop")

G          <- nrow(years_per_series) #56 rows
T_i        <- years_per_series$T_i #different for every stands and species (56 rows)
start_idxs <- cumsum(c(1, T_i[-G])) #cumulative sum: ragged vector 
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
#add a stand parameter 

# Fitting Model -----------------------------------------------------------

fit_all <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesC2.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)





# 4. Diagnostics --------------------------------------------------------

posterior_array <- as.array(fit_all)
divergences     <- nuts_params(fit_all)
np_style        <- pairs_style_np(div_color = "red", div_shape = 16, div_size = 1.5)

# Convert Stan fit to array
posterior_array <- as.array(fit_all)
# Get divergences
divergences <- nuts_params(fit_all)
pars_low <- c(
  "log_means[1,1]", "log_means[2,1]", "log_means[3,1]",      # species-level
  "log_alpha_low[1]", "log_alpha_low[2]", "log_alpha_low[3]", # stand-level
  "sigma_low"                                                 # SD controlling shrinkage
)
pars_low2 <- c(
  "log_means[4,1]", "log_means[5,1]",      # species-level
  "log_alpha_low[4]", "log_alpha_low[5]", # stand-level
  "sigma_low"                                                 # SD controlling shrinkage
)

np_style <- pairs_style_np(
  div_color = "red",  # divergent transitions highlighted in red
  div_shape = 16,     # point shape for divergences
  div_size = 1.5      # size of divergent points
)
mcmc_pairs(
  posterior_array,
  pars = pars_low2,
  np = divergences,
  off_diag_args = list(size = 1.5, alpha = 0.5),
  diag_fun = "dens",
  np_style = np_style
)


pars_high <- c(
  "log_means[1,2]", "log_means[2,2]", "log_means[3,2]",      # species-level
  "log_alpha_high[1]", "log_alpha_high[2]", "log_alpha_high[3]", # stand-level
  "sigma_high"                                                 # SD controlling shrinkage
)
pars_high2 <- c(
  "log_means[4,2]", "log_means[5,2]",      # species-level
  "log_alpha_high[4]", "log_alpha_high[5]", # stand-level
  "sigma_low"                                                 # SD controlling shrinkage
)

np_style <- pairs_style_np(
  div_color = "red",  # divergent transitions highlighted in red
  div_shape = 16,     # point shape for divergences
  div_size = 1.5      # size of divergent points
)
mcmc_pairs(
  posterior_array,
  pars = pars_high2,
  np = divergences,
  off_diag_args = list(size = 1.5, alpha = 0.5),
  diag_fun = "dens",
  np_style = np_style
)
mcmc_areas(
  posterior_array,
  pars = pars_high2,
  np = divergences,
  off_diag_args = list(size = 1.5, alpha = 0.5),
  diag_fun = "dens",
  np_style = np_style
)

pars_others <- c(
  "sigma_low", "sigma_high", "mu_log_phi1", "mu_log_phi2", 
  "sigma_log_phi2"                                              
)

mcmc_pairs(
  posterior_array,
  pars = pars_others,
  np = divergences,
  off_diag_args = list(size = 1.5, alpha = 0.5),
  diag_fun = "dens",
  np_style = np_style
)




pars_phi <- c(
  "mu_log_phi1", "log_phi1", "log_phi2","sigma_log_phi1"
)

np_style <- pairs_style_np(
  div_color = "red",  # divergent transitions highlighted in red
  div_shape = 16,     # point shape for divergences
  div_size = 1.5      # size of divergent points
)
mcmc_pairs(
  posterior_array,
  pars = pars_phi,
  np = divergences,
  off_diag_args = list(size = 1.5, alpha = 0.5),
  diag_fun = "dens",
  np_style = np_style
)


#Checking nESS ratio (efficiency of sampler -- values <1% bad as it means that draws were redundant or hiigly correlated)
cat("\nESS ratio < 0.1 (excluding state and log_omega):\n")
neff_vals <- neff_ratio(fit_all)
bad_neff  <- neff_vals[!is.na(neff_vals) & neff_vals < 0.1]
print(bad_neff[!grepl("state|log_omega|y_rep", names(bad_neff))])


#Old style Based on Mike's code 

util<- new.env()
#creating a source file to make all the plots. 
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)
samples <- util$extract_expectand_vals(fit_all)
par(mfrow = c(1, 1))
util$plot_hist_quantiles(
  samples, "y_rep",
  bin_min        = 0,
  bin_max        = 600,
  bin_delta      = 20,
  baseline_values = stan_data_all$y
)

#Probability of being in mast or not
par(mfrow = c(4, 3), mar = c(5, 5, 2, 1))

for (i in 1:stan_data_all$F) {
  idxs <- stan_data_all$start_idxs[i]:stan_data_all$end_idxs[i]  # i not S
  names <- paste0("y_rep[", idxs, "]")
  
  # Get species and stand name for the title
  stand_info <- years_per_series[i, ]
  title <- paste0(stand_info$spp, " — ", stand_info$stand)
  
  util$plot_disc_pushforward_quantiles(
    samples,
    names,
    baseline_values = stan_data_all$y[idxs],
    main = title                               # adds informative title per panel
  )
}

#Scaled
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data_all$F){
  s <- perm[i]
  idxs <- stan_data_all$start_idxs[s]:stan_data_all$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples,names,baseline_values = stan_data_all$y[idxs],display_ylim = c(0,400),main = paste0(stand_order[i], " (",       stand_elev$elevation[stand_elev$stand == stand_order[i]]," m)")
  )
}


# New plotting ------------------------------------------------------------

#Extracting everything
posterior_draws  <- rstan::extract(fit_all)
y_rep            <- posterior_draws$y_rep          # [draws × N]
state_draws      <- posterior_draws$state          # [draws × N]
log_means_draws  <- posterior_draws$log_means      # [draws × S × 2]
theta1_draws     <- posterior_draws$theta1         # [draws × S]
theta2_draws     <- posterior_draws$theta2         # [draws × S]
sigma_low_draws      <- posterior_draws$sigma_low          # [draws]
sigma__high_draws      <- posterior_draws$sigma_high         # [draws]
rho_draws        <- posterior_draws$rho            # [draws × 2]
log_phi1_draws   <- posterior_draws$log_phi1       # [draws × S]
log_phi2_draws   <- posterior_draws$log_phi2       # [draws × S]
log_alpha_draws  <- posterior_draws$log_alpha      # [draws × F]

species_names <- levels(as.factor(stand_year_all$spp))
# Should be: ABAM CANO PSME THPL TSHE

#PPC
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

#MAST VS NON_MAST 
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


#MASTING TIMELINE 
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

#TRANSITIONS
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


#Seed production
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






# For Lizzie --------------------------------------------------------------

# ── Prior vs Posterior: multispecies HMM ──────────────────────────────────────
species_names <- c("ABAM", "CANO", "PSME", "THPL", "TSHE")   # adjust if needed
S <- length(species_names)
n_prior <- 1e4

#Samples from priors
prior_list <- list()

# To match based on priors in .stan
log_means_priors <- list(
  ABAM = list(low = list(mu = 1.6,   sd = 0.8), high = list(mu = 4.1, sd = 0.8)),
  CANO = list(low = list(mu = 2.3, sd = 0.8),  high = list(mu = 5.7, sd = 0.8)),
  PSME = list(low = list(mu = 2.3,   sd = 1), high = list(mu = 5, sd = 1)),
  THPL = list(low = list(mu = 3,   sd = 0.8),  high = list(mu = 6.7, sd = 0.8)),
  TSHE = list(low = list(mu = 3.9,   sd = 1),  high = list(mu = 6.9, sd = 1))
)

for (s in seq_along(species_names)) {
  sp  <- species_names[s]
  pr  <- log_means_priors[[sp]]
  prior_list[[paste0("log_low_", sp)]]  <- rnorm(n_prior, pr$low$mu,  pr$low$sd)
  prior_list[[paste0("log_high_", sp)]] <- rnorm(n_prior, pr$high$mu, pr$high$sd)
}

# Thetha
for (s in seq_along(species_names)) {
  sp <- species_names[s]
  prior_list[[paste0("theta1_", sp)]] <- rbeta(n_prior, 4, 1)
  prior_list[[paste0("theta2_", sp)]] <- rbeta(n_prior, 1, 4)
}

# sigma
prior_list[["sigma_low"]] <- abs(rnorm(n_prior, 0, 0.5))   # half-normal (sigma > 0)
prior_list[["sigma_high"]] <- abs(rnorm(n_prior, 0, 0.5))   # half-normal (sigma > 0)

# # stand effects  (marginal: stand_effect_raw * sigma)
# prior_sigma <- abs(rnorm(n_prior, 0, 0.5))
# for (f in seq_len(stan_data_all$F)) {
#   prior_list[[paste0("stand_effect_", f)]] <- rnorm(n_prior, 0, 1) * prior_sigma
# }

# log_phi1, log_phi2 — stored on log scale in Stan, plot phi = exp(log_phi)
prior_list[["phi1"]] <- exp(rnorm(n_prior, log(6.5), 0.5))
prior_list[["phi2"]] <- exp(rnorm(n_prior, log(6.5), 0.5))

prior_df <- as.data.frame(prior_list)

#  2. Extract posteriors 
post_samples <- rstan::extract(fit_all)

post_list <- list()

# log_means[s, 1] = low state, log_means[s, 2] = high state
for (s in seq_along(species_names)) {
  sp <- species_names[s]
  post_list[[paste0("log_low_",  sp)]] <- post_samples$log_means[, s, 1]
  post_list[[paste0("log_high_", sp)]] <- post_samples$log_means[, s, 2]
}

# theta1[s], theta2[s]
for (s in seq_along(species_names)) {
  sp <- species_names[s]
  post_list[[paste0("theta1_", sp)]] <- post_samples$theta1[, s]
  post_list[[paste0("theta2_", sp)]] <- post_samples$theta2[, s]
}

post_list[["sigma_low"]] <- post_samples$sigma_low
post_list[["sigma_high"]] <- post_samples$sigma_high

# for (f in seq_len(stan_data_all$F)) {
#   post_list[[paste0("stand_effect_", f)]] <- post_samples$stand_effect_raw[, f] *
#     post_samples$sigma
# }

# Back-transform dispersion to natural scale for interpretability
post_list[["phi1"]] <- exp(post_samples$log_phi1[, 1])   # or use rowMeans if you want species mean
post_list[["phi2"]] <- exp(post_samples$log_phi2[, 1])

posterior_df <- as.data.frame(post_list)

# Reshape to long format
to_long <- function(df, label) {
  stack_df <- utils::stack(df)          # base R — no tidyr dependency
  names(stack_df) <- c("value", "parameter")
  stack_df$distribution <- label
  stack_df
}

prior_long     <- to_long(prior_df,     "Prior")
posterior_long <- to_long(posterior_df, "Posterior")
combined       <- rbind(prior_long, posterior_long)

# PLot
ggplot(combined, aes(x = value, colour = distribution)) +
  geom_density(linewidth = 0.7) +
  facet_wrap(~parameter, scales = "free") +
  scale_colour_manual(values = c("Prior" = "steelblue", "Posterior" = "tomato")) +
  labs(
    title  = "Prior vs Posterior — multispecies HMM",
    x      = "Parameter value",
    y      = "Density",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(size = 7),
    legend.position  = "top"
  )



# New diagnostics  --------------------------------------------------------


# See R-hat and ESS for every parameter
print(fit_all, pars = c("log_means", "theta1", "theta2", "sigma_low","sigma_high", "log_phi1", "log_phi2"))

# Which parameters have the worst R-hat?
summary_fit <- summary(fit_all)$summary
summary_fit[is.na(summary_fit[,"Rhat"]) | summary_fit[,"Rhat"] > 1.05, ]


# For presentation --------------------------------------------------------
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── 1. Extract draws ──────────────────────────────────────────────────────────

draws       <- as.data.frame(fit_all)  # fit = your stanfit object
sampler_par <- get_sampler_params(fit_all, inc_warmup = FALSE)

# Divergence flag — one row per draw (chains concatenated)
divergent <- do.call(rbind, sampler_par)[, "divergent__"] == 1

# Parameters of interest
mu_phi1         <- draws$mu_log_phi1
sigma_phi1      <- draws$sigma_log_phi1

# Species-level: log_phi1[1], log_phi1[2], ...  ── adjust S to your number of species
S <- 5
log_phi1_cols <- paste0("log_phi1[", 1:S, "]")
log_phi1_mat  <- draws[, log_phi1_cols]

# ── 2. Build long data frames ─────────────────────────────────────────────────

# Plot A: sigma_log_phi1 vs mu_log_phi1
df_hyper <- data.frame(
  mu_log_phi1    = mu_phi1,
  sigma_log_phi1 = sigma_phi1,
  divergent      = divergent
)

# Plot B: sigma_log_phi1 vs each log_phi1[s]
df_species <- log_phi1_mat |>
  mutate(sigma_log_phi1 = sigma_phi1,
         divergent      = divergent) |>
  pivot_longer(
    cols      = all_of(log_phi1_cols),
    names_to  = "species",
    values_to = "log_phi1"
  ) |>
  mutate(species = gsub("log_phi1\\[|\\]", "", species),
         species = paste0("Species ", species))

# ── 3. Shared theme ───────────────────────────────────────────────────────────

div_colors <- c("FALSE" = "steelblue", "TRUE" = "red")
div_labels <- c("FALSE" = "No divergence", "TRUE" = "Divergence")

base_theme <- theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

# ── 4. Plot A: funnel check — sigma vs mu ─────────────────────────────────────

p_hyper <- ggplot(df_hyper, aes(mu_log_phi1, sigma_log_phi1,
                                colour = factor(divergent))) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_colour_manual(values = div_colors, labels = div_labels,
                      name = NULL) +
  labs(
    title    = "Hyperparameter funnel: σ_log_phi1 vs μ_log_phi1",
    subtitle = "Divergences in orange — funnel neck at low σ is normal",
    x        = expression(mu[log~phi[1]]),
    y        = expression(sigma[log~phi[1]])
  ) +
  base_theme

# ── 5. Plot B: sigma vs each species log_phi1[s] ─────────────────────────────

p_species <- ggplot(df_species, aes(log_phi1, sigma_log_phi1,
                                    colour = factor(divergent))) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_colour_manual(values = div_colors, labels = div_labels,
                      name = NULL) +
  facet_wrap(~ species, nrow = 1) +
  labs(
    title    = "Species-level funnel: σ_log_phi1 vs log_phi1[s]",
    subtitle = "One panel per species — checking partial pooling geometry",
    x        = expression(log~phi[1][s]),
    y        = expression(sigma[log~phi[1]])
  ) +
  base_theme

# ── 6. Print ──────────────────────────────────────────────────────────────────

print(p_hyper)
print(p_species)


# New stuff ---------------------------------------------------------------

draws_df <- as.data.frame(fit_all)
div      <- do.call(rbind, get_sampler_params(fit_all, inc_warmup = FALSE))[, "divergent__"] == 1

# Check phi hierarchy
mcmc_scatter(as.matrix(fit_all),
             pars   = c("sigma_log_phi1", "log_phi1[1]"),
             np     = nuts_params(fit_all)) +
  ggtitle("phi1: centered funnel check")

# Check alpha hierarchy  
mcmc_scatter(as.matrix(fit_all),
             pars   = c("sigma_low", "log_alpha_low[1]"),
             np     = nuts_params(fit_all)) +
  ggtitle("alpha_low: centered funnel check")



# New 2 -------------------------------------------------------------------

draws <- as_draws_array(fit_all)


print(fit_all, pars = c("mu_log_phi1", "mu_log_phi2", "sigma_phi",
                    "logit_mu_theta1", "logit_mu_theta2", "sigma_theta",
                    "sigma_low", "sigma_high"))


print(fit_all)

check_hmc_diagnostics(fit_all)

traceplot(fit_all, pars = c("sigma_phi", "sigma_theta", "sigma_low", "sigma_high"))




pairs(fit_all, pars = c("sigma_phi", "mu_log_phi1[1]", "mu_log_phi2[1]", "phi1_raw[1]", "phi2_raw[1]"))
pairs(fit_all, pars = c("sigma_phi", "mu_log_phi1[2]", "mu_log_phi2[2]", "phi1_raw[2]", "phi2_raw[2]"))
pairs(fit_all, pars = c("sigma_theta", "logit_mu_theta1[1]","logit_mu_theta2[1]","theta1_raw[1]","theta2_raw[1]" ))
pairs(fit_all, pars = c("sigma_theta", "logit_mu_theta1[2]","logit_mu_theta2[2]","theta1_raw[2]","theta2_raw[2]" ))
pairs(fit_all, pars = c("sigma_low[1]","stand_effect_low_raw[1]","sigma_low[2]","stand_effect_low_raw[2]","sigma_low[3]","stand_effect_low_raw[3]"))
pairs(fit_all, pars = c("sigma_low[4]","stand_effect_low_raw[4]","sigma_low[5]","stand_effect_low_raw[5]"))
pairs(fit_all, pars = c("sigma_high[1]","stand_effect_high_raw[1]","sigma_high[2]","stand_effect_high_raw[2]","sigma_high[3]","stand_effect_high_raw[3]"))
pairs(fit_all, pars = c("sigma_high[4]","stand_effect_high_raw[4]","sigma_high[5]","stand_effect_high_raw[5]"))


mcmc_trace(
  draws,
  pars = c("sigma_low", "sigma_high",
           "theta1[1]", "theta2[1]")
)

mcmc_trace(
  draws,
  pars = c("log_means[1,1]","log_means[1,2]", "log_means[2,1]","log_means[2,2]",
           "log_means[3,1]","log_means[3,2]"))

mcmc_trace(
  draws,
  pars = c("log_means[4,1]","log_means[4,2]", "log_means[5,1]","log_means[5,2]"))


mcmc_pairs(
  draws,
  pars = c("sigma_high", "stand_effect_high_raw[1]", "log_means[1,2]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("sigma_high", "stand_effect_high_raw[1]", "log_means[1,2]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("sigma_low", "stand_effect_low_raw[1]", "log_means[1,1]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("log_phi1[1]", "log_phi2[1]", "sigma_low", "sigma_high","theta1[1]","theta2[1]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("log_phi1[1]", "log_phi2[1]", "sigma_low", "sigma_high","theta1[1]","theta2[1]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("log_phi1[4]", "log_phi2[3]", "log_phi2[4]", "log_phi2[5]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("sigma_low", "sigma_high","theta1[1]","theta2[1]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("log_phi1[1]", "log_phi2[1]", "theta1[1]","theta2[1]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("sigma_low", "sigma_high","log_phi1[1]","log_phi2[1]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("theta1[1]","theta2[1]","theta1[2]","theta2[2]","theta1[3]","theta2[3]","theta1[4]","theta2[4]","theta1[5]","theta2[5]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("log_phi1[1]","log_phi2[1]","log_phi1[2]","log_phi2[2]","log_phi1[3]","log_phi2[3]","log_phi1[4]","log_phi2[4]","log_phi1[5]","log_phi2[5]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("log_means[5,1]", "log_means[5,2]",
           "log_phi1[5]", "log_phi2[5]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("theta1[5]", "theta2[5]"),
  np = nuts_params(fit_all)
)

mcmc_pairs(
  draws,
  pars = c("sigma_low","stand_effect_low_raw[1]",
           "sigma_high","stand_effect_high_raw[1]"),
  np = nuts_params(fit_all)
)
mcmc_pairs(
  draws,
  pars = c("sigma_low","stand_effect_low_raw[2]",
           "sigma_high","stand_effect_high_raw[2]"),
  np = nuts_params(fit_all)
)
mcmc_pairs(
  draws,
  pars = c("sigma_low","stand_effect_low_raw[3]",
           "sigma_high","stand_effect_high_raw[3]"),
  np = nuts_params(fit_all)
)
mcmc_pairs(
  draws,
  pars = c("sigma_low","stand_effect_low_raw[4]",
           "sigma_high","stand_effect_high_raw[4]"),
  np = nuts_params(fit_all)
)
mcmc_pairs(
  draws,
  pars = c("sigma_low","stand_effect_low_raw[5]",
           "sigma_high","stand_effect_high_raw[5]"),
  np = nuts_params(fit_all)
)

library(posterior)

summ <- summarise_draws(draws)

# Rhat distribution
hist(summ$rhat, breaks = 50, main = "R-hat distribution")

# ESS
par(mfrow=c(1,2))
hist(summ$ess_bulk, breaks = 50, main = "Bulk ESS distribution")
hist(summ$ess_tail, breaks = 50, main = "Tail ESS distribution")

np <- nuts_params(fit_all)

mcmc_pairs(
  draws,
  pars = c("sigma_low", "stand_effect_low_raw[1]"),
  np = np
)

mcmc_dens_overlay(
  draws,
  pars = c("log_means[1,1]", "log_means[1,2]")
)

mcmc_intervals(
  draws,
  pars = paste0("stand_effect_high_raw[", 1:10, "]")
)

y_rep <- as.matrix(fit_all, pars = "y_rep")

ppc_dens_overlay(
  y = stan_data_all$y,
  yrep = y_rep[1:100, ]
)


draws_df <- as_draws_df(fit_all)
draws_df$diff_1 <- draws_df$`log_means[1,2]` - draws_df$`log_means[1,1]`

hist(draws_df$diff_1, breaks = 50)


# Prior only
fit_prior <- sampling(
  model,
  data = c(stan_data_all, list(prior_only = 1)),
  iter = 1000,
  chains = 4
)

# Posterior
fit_post <- sampling(
  model,
  data = c(stan_data_all, list(prior_only = 0)),
  iter = 2000, chains = 4
)

library(bayesplot)

draws_prior <- as_draws_df(fit_prior)
draws_post  <- as_draws_df(fit_post)

mcmc_dens_overlay(
  list(
    prior = draws_prior,
    posterior = draws_post
  ),
  pars = c("sigma_low", "sigma_high")
)

library(ggplot2)

df <- rbind(
  data.frame(value = draws_prior$sigma_low, type = "prior"),
  data.frame(value = draws_post$sigma_low,  type = "posterior")
)

ggplot(df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.4)


mcmc_dens_overlay(
  draws,
  pars = c("log_means[5,1]","log_means[5,2]")
)
