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
library(purrr)
library(patchwork) 

options(mc.cores = parallel::detectCores())




#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

list.files()
list.files("Data/")
list.files("Rscript_AnalysisEleonore/")


util <- new.env()
source("Rscript_AnalysisEleonore/mcmc_analysis_tools_rstan.R", local=util)
source("Rscript_AnalysisEleonore/mcmc_visualization_tools.R", local=util)

# Data import -------------------------------------------------------------


seed_data<-read.csv("Data/SeedData_all.csv")

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
#Current model used after discussion in git
fit_75 <- stan(
  file    = "Stan_code_Eleonore/Species_Stan_Model/MultispeciesStanCodeLast.stan",
  data    = stan_data_all,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

# PPC ---------------------------------------------------------------------

#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")
dev.off()

samples <- util$extract_expectand_vals(fit_75)

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

#Single observation

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


# Prior vs posterior ------------------------------------------------

samples <- rstan::extract(fit_75)
n <- length(samples$mu_log_low)

priors <- data.frame(
  mu_log_low           = rnorm(n, 2, 1.5),
  mu_log_delta         = rnorm(n, 1.5, 1.5),
  grand_logit_theta1   = rnorm(n, 0.5, 1),
  grand_logit_theta2   = rnorm(n, -1, 0.5),
  sigma_low_species    = abs(rnorm(n, 0.5, 1)),
  sigma_stand          = abs(rnorm(n,1,0.5)),
  sigma_log_delta_species   = abs(rnorm(n, 1.5, 0.7)),
  sigma_delta_stand    = abs(rnorm(n,1,0.5)),
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
  sigma_log_delta_species   = samples$sigma_log_delta_species,
  sigma_delta_stand    = samples$sigma_delta_stand,
  sigma_theta1_species = samples$sigma_theta1_species,
  sigma_theta2_species = samples$sigma_theta2_species,
  sigma_theta1_stand   = samples$sigma_theta1_stand,
  sigma_theta2_stand   = samples$sigma_theta2_stand,
  phi_species          = exp(rowMeans(samples$log_phi_species)),
  phi_low              = exp(samples$log_phi_state[, 1]),
  phi_high             = exp(samples$log_phi_state[, 2])
)

#labels 
param_labels <- c(
  mu_log_low              = "grand mean seeds, low state",
  mu_log_delta            = "grand mean fold-change",
  grand_logit_theta1      = "P stay low",
  grand_logit_theta2      = "P stay high",
  sigma_low_species       = "σ low species",
  sigma_stand             = "σ low stand",
  sigma_log_delta_species = "σ log delta species",
  sigma_delta_stand       = "σ log delta stand",
  sigma_theta1_species    = "σ θ₁ species",
  sigma_theta2_species    = "σ θ₂ species",
  sigma_theta1_stand      = "σ θ₁ stand",
  sigma_theta2_stand      = "σ θ₂ stand",
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


# Masting synchrony analysis-------------------------------------------------------------------------
# Species x Stand x Year latent state (1/2) synchrony analysis, based on posterior draws of a latent "state" variable from a fitted Stan model (fit_75).
# This script loads the data once, defines each helper function once, and then runs through every analysis/plot in a single pass.

#Data loading 
# Elevation lookup table, one row per stand, used later to plot synchrony against elevation
stand_elevation <- readRDS("Data/stand_elevation_table.rds")

# Posterior draws of the latent state for every species x stand x year unit.
# state_draws is a matrix: rows = posterior draws, columns = species-stand-year
# units (in the same order as stand_year_all)
state_draws <- rstan::extract(fit_75, pars = "state")$state
n_draws     <- nrow(state_draws)

# Which draws to actually loop over. Using every draw can be slow, so I can subsample instead, for example draw_idx_use <- seq(1, n_draws, by = 5)
draw_idx_use <- seq_len(n_draws)

# Lookup table linking each column of state_draws back to its species, stand, and year. row_id is the column index into state_draws for that unit.
meta <- stand_year_all %>%
  select(species = spp, stand = stand, year = year) %>%
  mutate(row_id = row_number())

# This pulls out one posterior draw and attaches the corresponding latent state to every species-stand-year observation in meta
get_unit_year_state <- function(draw_idx) {
  meta %>%
    mutate(state = state_draws[draw_idx, row_id])
}

# Same as above, but also flags whether each species-stand unit switched state from the previous year (used for the "mast switching" transition analyses)
get_transition_data <- function(draw_idx) {
  get_unit_year_state(draw_idx) %>%
    arrange(species, stand, year) %>%
    group_by(species, stand) %>%
    mutate(transition = state != lag(state)) %>%
    ungroup()
}

#Helper functions that I use across all analyses below
# Generic pairwise synchrony: for any group of rows (each with a "state"),returns the proportion of all pairwise combinations that share the same state. Used for interspecific (within stand), intraspecific (within species)
# NOTE ON INTERPRETATION (why synchrony often sits near/below 50%): state is binary, so under INDEPENDENT random assignment the expected synchrony is p^2 + (1-p)^2, which is minimized at p = 0.5, giving a floor of exactly 50%. Synchrony can't drop below 50% from chance pairing alone when state frequencies are roughly balanced, so values BELOW 50% indicate real ANTI-synchrony (pairs tend to be in opposite states more than chance would predict), not just "weak" synchrony. I could consider a permutation-based null (shuffle state labels within stand-year) if you want a formal chance line instead of a flat 50% reference.
pairwise_synchrony_generic <- function(df_one_group) {
  n <- nrow(df_one_group)
  if (n < 2) return(tibble(synchrony = NA_real_))
  idx <- combn(n, 2)
  same_state <- df_one_group$state[idx[1, ]] == df_one_group$state[idx[2, ]]
  tibble(synchrony = mean(same_state))
}

# Same idea as above, but computes within-stand and between-stand synchrony together for a single year's worth of data, so both numbers come from the same set of pairwise comparisons
pairwise_synchrony_within_between <- function(df_one_year) {
  n <- nrow(df_one_year)
  if (n < 2) return(tibble(within = NA_real_, between = NA_real_))
  idx        <- combn(n, 2)
  same_state <- df_one_year$state[idx[1, ]] == df_one_year$state[idx[2, ]]
  same_stand <- df_one_year$stand[idx[1, ]] == df_one_year$stand[idx[2, ]]
  tibble(
    within  = if (any(same_stand))  mean(same_state[same_stand])  else NA_real_,
    between = if (any(!same_stand)) mean(same_state[!same_stand]) else NA_real_
  )
}

# Proportion of pairs that BOTH switched state in the same year ( I could here maybe select the state?)
pairwise_transition_synchrony <- function(df_one_year) {
  n <- nrow(df_one_year)
  if (n < 2) return(tibble(sync = NA_real_))
  idx         <- combn(n, 2)
  both_switch <- df_one_year$transition[idx[1, ]] & df_one_year$transition[idx[2, ]]
  tibble(sync = mean(both_switch, na.rm = TRUE))
}

# Collapses a data frame of per-draw values into a posterior median plus 50% and 90% credible intervals, grouped by whatever variables you give it ( for example year, or species + year, or stand)
summarise_draws <- function(df, value_col, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      median = median(.data[[value_col]], na.rm = TRUE),
      lo50   = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      hi50   = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      lo90   = quantile(.data[[value_col]], 0.05, na.rm = TRUE),
      hi90   = quantile(.data[[value_col]], 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

# Shared 50%/90% credible-interval ribbon layers, so every plot below uses the same visual style
ribbon_layers <- function(fill = "#4C7C9B") {
  list(
    geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), fill = fill, alpha = 0.15),
    geom_ribbon(aes(ymin = lo50 * 100, ymax = hi50 * 100), fill = fill, alpha = 0.35)
  )
}

#Overall synchrony over time, within stand vs. between stands
# For each posterior draw and year, I compute the proportion of species pairs sharing a state, separately for pairs in the same stand and pairs in different stands
synchrony_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    group_by(year) %>%
    group_modify(~ pairwise_synchrony_within_between(.x)) %>%
    mutate(draw = d)
})

synchrony_summary <- synchrony_draws %>%
  pivot_longer(c(within, between), names_to = "type", values_to = "synchrony") %>%
  summarise_draws("synchrony", c("year", "type")) %>%
  mutate(type = recode(type, within = "Within stand", between = "Between stands"))

# Plot: park-wide synchrony over time, comparing within-stand vs between-stand
ggplot(synchrony_summary, aes(x = year, y = median * 100, color = type, fill = type)) +
  geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), alpha = 0.15, color = NA) +
  geom_ribbon(aes(ymin = lo50 * 100, ymax = hi50 * 100), alpha = 0.35, color = NA) +
  geom_line(aes(group = type), linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Within stand" = "#7A3B69", "Between stands" = "#4C7C9B")) +
  scale_fill_manual(values  = c("Within stand" = "#7A3B69", "Between stands" = "#4C7C9B")) +
  labs(x = NULL, y = "Synchrony (% of species pairs)", color = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# same but just within stand
main_synchrony_summary <- synchrony_summary %>%
  filter(type == "Within stand")

ggplot(main_synchrony_summary, aes(x = year, y = median * 100)) +
  ribbon_layers(fill = "#7A3B69") +
  geom_line(color = "#5A2B4D", linetype = "dashed", linewidth = 0.6) +
  geom_point(color = "#5A2B4D", size = 2.2) +
  labs(x = NULL, y = "Interspecific synchrony (% of species pairs)",
       title = "Park-wide interspecific synchrony in seed production",
       subtitle = "Species pairs at the same stand, same year") +
  theme_minimal(base_size = 13)

#Interspecific synchrony along the elevational gradient : Within-stand synchrony (species pairs), averaged across all years, then related to each stand's elevation
stand_synchrony_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    group_by(stand, year) %>%
    filter(n() >= 2) %>%                        # need >= 2 species to compare
    group_modify(~ pairwise_synchrony_generic(.x)) %>%
    group_by(stand) %>%
    summarise(synchrony = mean(synchrony), .groups = "drop") %>%
    mutate(draw = d)
})

# Number of species ever recorded in each stand, used to flag low-confidence stands (few species => noisier synchrony estimate)
n_species_per_stand <- meta %>%
  group_by(stand) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

stand_synchrony_summary <- stand_synchrony_draws %>%
  summarise_draws("synchrony", "stand") %>%
  left_join(stand_elevation %>% distinct(stand, elevation), by = "stand") %>%
  left_join(n_species_per_stand, by = "stand") %>%
  mutate(
    n_species_bin = case_when(
      n_species <= 2 ~ "1-2",
      n_species <= 4 ~ "3-4",
      n_species <= 6 ~ "5-6",
      TRUE           ~ "7+"
    ),
    n_species_bin = factor(n_species_bin, levels = c("1-2", "3-4", "5-6", "7+"))
  ) %>%
  arrange(elevation)

# Plot: interspecific synchrony vs. elevation, point color shows how many species that stand's estimate is based on
ggplot(stand_synchrony_summary, aes(x = elevation, y = median * 100)) +
  ribbon_layers() +
  geom_line(color = "#1B4C6B", linewidth = 0.8) +
  geom_point(aes(color = n_species_bin), size = 2.5) +
  scale_color_manual(
    values = c("1-2" = "#D9695F", "3-4" = "#E8A54B", "5-6" = "#8FBF7F", "7+" = "#2E6B4F"),
    name = "# species"
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
  labs(x = "Elevation (m)", y = "Synchrony index (within stand, % of pairs)",
       title = "Interspecific Synchrony Along the Elevational Gradient",
       subtitle = "Dark ribbon = 50% CI | Light ribbon = 90% CI | Line = posterior median\nDashed grey line = 50% chance floor under independent, balanced states") +
  theme_minimal(base_size = 13)

#Within-stand synchrony over time, broken out by stand 
# Same computation as above (interspecific synchrony within a stand), but kept at the stand x year level instead of averaged over years, so you can see the trajectory of each stand separately
# Each point below = % of species pairs *in that stand and year* sharing a state. Because state is binary, chance/independent pairing gives an EXPECTED synchrony of exactly 50% when state frequencies are balanced (p ~ 0.5), and this is a FLOOR, not a midpoint so imbalance only pushes the chance expectation UP, never down. Points sitting below the dashed 50% line therefore reflect species actively alternating states (anti-synchrony) among co-occurring species in that stand.
stand_synchrony_time_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    group_by(stand, year) %>%
    filter(n() >= 2) %>%  # need >= 2 species present that year
    group_modify(~ pairwise_synchrony_generic(.x)) %>%
    ungroup() %>%
    mutate(draw = d)
})

n_species_per_stand_year <- meta %>%
  group_by(stand, year) %>%
  summarise(n_species = n_distinct(species), .groups = "drop")

stand_synchrony_time_summary <- stand_synchrony_time_draws %>%
  summarise_draws("synchrony", c("stand", "year")) %>%
  left_join(n_species_per_stand_year, by = c("stand", "year")) %>%
  mutate(
    n_species_bin = case_when(
      n_species <= 2 ~ "1-2",
      n_species <= 4 ~ "3-4",
      n_species <= 6 ~ "5-6",
      TRUE           ~ "7+"
    ),
    n_species_bin = factor(n_species_bin, levels = c("1-2", "3-4", "5-6", "7+"))
  )

year_range <- range(meta$year)

# Plot: one facet per stand, showing within-stand synchrony over time, with the 50% chance floor drawn in for reference
ggplot(stand_synchrony_time_summary, aes(x = year, y = median * 100)) +
  ribbon_layers(fill = "#7A3B69") +
  geom_line(color = "#5A2B4D", linewidth = 0.6) +
  geom_point(aes(color = n_species_bin), size = 2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
  scale_color_manual(
    values = c("1-2" = "#D9695F", "3-4" = "#E8A54B", "5-6" = "#8FBF7F", "7+" = "#2E6B4F"),
    name = "# species"
  ) +
  scale_x_continuous(limits = year_range, breaks = seq(year_range[1], year_range[2], 2)) +
  facet_wrap(~ stand) +
  labs(x = NULL, y = "Synchrony across species (% of pairs)",
       title = "Within-stand synchrony over time, by stand",
       subtitle = "Dashed line = 50% chance floor; values below it indicate anti-synchrony") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Synchrony of state transitions ("mast switching") over time 
# Park-wide: for species pairs in the same stand, the proportion that BOTH switched state (rather than just sharing a state) in the same year, averaged across stands within each year
transition_sync_draws <- map_dfr(draw_idx_use, function(d) {
  get_transition_data(d) %>%
    group_by(stand, year) %>%              # pairs restricted to the same stand = interspecific
    filter(n() >= 2) %>%
    group_modify(~ pairwise_transition_synchrony(.x)) %>%
    ungroup() %>%
    group_by(year) %>%
    summarise(sync = mean(sync, na.rm = TRUE), .groups = "drop") %>%
    mutate(draw = d)
})

transition_sync_summary <- transition_sync_draws %>%
  summarise_draws("sync", "year")

ggplot(transition_sync_summary, aes(x = year, y = median * 100)) +
  ribbon_layers() +
  geom_line(color = "#1B4C6B", linewidth = 0.8) +
  geom_point(color = "#1B4C6B", size = 2) +
  labs(x = "Year", y = "Transition synchrony (% of species pairs)",
       title = "Synchrony of state transitions, park-wide",
       subtitle = "Proportion of species pairs (same stand) switching state together") +
  theme_minimal(base_size = 13)

# By stand: same transition synchrony metric, kept separate per stand instead of averaged into one park-wide number
stand_transition_draws <- map_dfr(draw_idx_use, function(d) {
  get_transition_data(d) %>%
    group_by(stand, year) %>%
    filter(n() >= 2) %>%
    group_modify(~ pairwise_transition_synchrony(.x)) %>%
    ungroup() %>%
    mutate(draw = d)
})

stand_transition_summary <- stand_transition_draws %>%
  summarise_draws("sync", c("stand", "year")) %>%
  left_join(n_species_per_stand, by = "stand") %>%
  mutate(
    n_species_bin = case_when(
      n_species <= 2 ~ "1-2",
      n_species <= 4 ~ "3-4",
      n_species <= 6 ~ "5-6",
      TRUE           ~ "7+"
    ),
    n_species_bin = factor(n_species_bin, levels = c("1-2", "3-4", "5-6", "7+"))
  )

ggplot(stand_transition_summary, aes(x = year, y = median * 100)) +
  ribbon_layers(fill = "#7A3B69") +
  geom_line(color = "#5A2B4D", linewidth = 0.7) +
  geom_point(aes(color = n_species_bin), size = 2.3) +
  scale_color_manual(
    values = c("1-2" = "#D9695F", "3-4" = "#E8A54B", "5-6" = "#8FBF7F", "7+" = "#2E6B4F"),
    name = "# species"
  ) +
  facet_wrap(~ stand, scales = "free_x") +
  labs(x = NULL, y = "Transition synchrony across species (% of pairs)",
       title = "Interspecific synchrony of switching events, by stand") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

#Intraspecific synchrony: same species, compared across stands 
# For each species and year, proportion of stand pairs sharing the same state
species_synchrony_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    group_by(species, year) %>%
    filter(n() >= 2) %>% # need >= 2 stands present that year
    group_modify(~ pairwise_synchrony_generic(.x)) %>%
    ungroup() %>%
    mutate(draw = d)
})

# Number of stands each species is recorded in per year, used to flag low-confidence species-year estimates (few stands => noisier)
n_stands_per_species_year <- meta %>%
  group_by(species, year) %>%
  summarise(n_stands = n_distinct(stand), .groups = "drop")

species_synchrony_summary <- species_synchrony_draws %>%
  summarise_draws("synchrony", c("species", "year")) %>%
  left_join(n_stands_per_species_year, by = c("species", "year")) %>%
  mutate(
    n_stands_bin = case_when(
      n_stands %in% 1:4   ~ "1-4",
      n_stands %in% 5:10  ~ "5-9",
      n_stands %in% 11:16 ~ "10-16",
      n_stands >= 17      ~ "17"
    ),
    n_stands_bin = factor(n_stands_bin, levels = c("1-4", "5-9", "10-16", "17"))
  )

# Plot: one facet per species, showing how synchronized that species is across the stands it occurs in, over time
ggplot(species_synchrony_summary, aes(x = year, y = median * 100)) +
  ribbon_layers() +
  geom_line(color = "#1B4C6B", linewidth = 0.7) +
  geom_point(aes(color = n_stands_bin), size = 2.3) +
  scale_color_manual(
    values = c("1-4" = "#D9695F", "5-9" = "#E8A54B", "10-16" = "#8FBF7F", "17" = "#1B4C6B"),
    name = "# stands"
  ) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = NULL, y = "Synchrony across stands (% of pairs)",
       title = "Intraspecific synchrony across sites, by species") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

# Same idea, but for state TRANSITIONS: does a species switch state at the same time across the different stands it's found in?
species_transition_draws <- map_dfr(draw_idx_use, function(d) {
  get_transition_data(d) %>%
    group_by(species, year) %>%             # same species, across stands = intraspecific
    filter(n() >= 2) %>%
    group_modify(~ pairwise_transition_synchrony(.x)) %>%
    ungroup() %>%
    mutate(draw = d)
})

species_transition_summary <- species_transition_draws %>%
  summarise_draws("sync", c("species", "year")) %>%
  left_join(n_stands_per_species_year, by = c("species", "year")) %>%
  mutate(
    n_stands_bin = case_when(
      n_stands %in% 1:4   ~ "1-4",
      n_stands %in% 5:10  ~ "5-9",
      n_stands %in% 11:16 ~ "10-16",
      n_stands >= 17      ~ "17"
    ),
    n_stands_bin = factor(n_stands_bin, levels = c("1-4", "5-9", "10-16", "17"))
  )

ggplot(species_transition_summary, aes(x = year, y = median * 100)) +
  ribbon_layers(fill = "#7A3B69") +
  geom_line(color = "#5A2B4D", linewidth = 0.7) +
  geom_point(aes(color = n_stands_bin), size = 2.3) +
  scale_color_manual(
    values = c("1-4" = "#D9695F", "5-9" = "#E8A54B", "10-16" = "#8FBF7F", "17" = "#1B4C6B"),
    name = "# stands"
  ) +
  facet_wrap(~ species, scales = "free_x") +
  labs(x = NULL, y = "Transition synchrony across stands (% of pairs)",
       title = "Intraspecific synchrony of switching events, by species") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))





# Comparing to Victor's code ----------------------------------------------
# Alternative synchrony metric: "majority-state" synchrony
# Adapted from Victor's script written for individual trees nested in sites, to my species-within-stand-within-year 
# Instead of the pairwise "do two species/stands share the same state" metric used in my synchrony analysis, this asks a slightly different question: within a stand-year, what fraction of species sit in whichever state is more common that year? That's a measure of how consolidated (vs. split) the stand is, and the two metrics are related but not identical - for a binary state with proportion p in the majority category, pairwise same-state synchrony = p^2 + (1-p)^2, while majority-fraction synchrony = p. 

masting_state    <- 2
nonmasting_state <- 1

# For one posterior draw, count how many species are masting/non-masting in each stand-year, and I turn those counts into fractions
compute_state_counts <- function(draw_idx) {
  get_unit_year_state(draw_idx) %>%
    group_by(stand, year) %>%
    filter(n() >= 2) %>%# need >= 2 species to define a "majority"
    summarise(
      n_total      = n(),
      n_masting    = sum(state == masting_state),
      n_nonmasting = sum(state == nonmasting_state),
      .groups = "drop"
    ) %>%
    mutate(
      masting_frac    = n_masting / n_total,
      nonmasting_frac = n_nonmasting / n_total,
      # fraction of species in the stand-year that are in the MAJORITY state, regardless of which state that happens to be
      majority_frac   = pmax(n_masting, n_nonmasting) / n_total
    )
}

# For each draw, average the stand-level fractions into two yearly numbers:
#  within stand synchrony: mean, across stands, of how consolidated each stand is into its own majority state  
#  across stands synchrony: mean, across stands, of the overall masting fraction and non-masting fraction, then the max of those two means 
yearly_metrics_draws <- map_dfr(draw_idx_use, function(d) {
  compute_state_counts(d) %>%
    group_by(year) %>%
    summarise(
      withinstand_synchrony       = mean(majority_frac),
      acrossstands_masting_frac   = mean(masting_frac),
      acrossstands_nonmasting_frac = mean(nonmasting_frac),
      .groups = "drop"
    ) %>%
    mutate(
      acrossstands_synchrony = pmax(acrossstands_masting_frac, acrossstands_nonmasting_frac),
      draw = d
    )
})

# Posterior median + 50%/90% credible intervals, by year, for both metrics
yearly_metrics_summary <- yearly_metrics_draws %>%
  select(year, draw, withinstand_synchrony, acrossstands_synchrony) %>%
  pivot_longer(c(withinstand_synchrony, acrossstands_synchrony),
               names_to = "type", values_to = "synchrony") %>%
  summarise_draws("synchrony", c("year", "type")) %>%
  mutate(type = recode(type,
                       withinstand_synchrony  = "Within stand",
                       acrossstands_synchrony = "Across stands"))

#same visual style as the rest 
ggplot(yearly_metrics_summary, aes(x = year, y = median * 100, color = type, fill = type)) +
  geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), alpha = 0.15, color = NA) +
  geom_ribbon(aes(ymin = lo50 * 100, ymax = hi50 * 100), alpha = 0.35, color = NA) +
  geom_line(aes(group = type), linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Within stand" = "#4f1d4b", "Across stands" = "#1b3f5a")) +
  scale_fill_manual(values  = c("Within stand" = "#bba6bb", "Across stands" = "#b7cddb")) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
  labs(x = NULL, y = "Majority-state synchrony (%)", color = NULL, fill = NULL,
       title = "Majority-state synchrony over time (yearly resolution)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")


# Other plots that I need to do : -----------------------------------------

# Analysis of Synchrony : ONE species -------------------------------------
# ABAM 
target_species <- "ABAM"


# STEP 1: pull the latent state draws out of the fit
# `state` in generated quantities is an array of length N (one entry per observation), sampled once per posterior draw. rstan::extract() merges all chains and returns a (n_draws x N) matrix - draws in rows, observations in columns - in the same observation order as the original `y` vector.
state_draws <- rstan::extract(fit_75, pars = "state")$state

# STEP 2: rebuild the observation-level metadata directly from stand_year_all
# Because `y`, `area`, and now `state` all share the same row order as stand_year_all, I don't need a separate row_id lookup as rhe row NUMBER in stand_year_all is the column index into state_draws. So here I re-attach the species/stand/year label to each column from the state draws
meta_all <- stand_year_all %>%
  mutate(row_id = row_number()) %>%
  rename(species = spp)   # renamed only for readability below

# STEP 3: small helper functions 
# Here I reconstruct a long-format (species, stand, year, state) table for one posterior draw, by combining the metadata with that draw's row of state_draws.
get_unit_year_state <- function(draw_idx) {
  tibble(
    species = meta_all$species,
    stand   = meta_all$stand,
    year    = meta_all$year,
    state   = state_draws[draw_idx, ]
  )
}

# Pairwise synchrony within one year, for one draw: the fraction of all stand PAIRS that are in the same state.
pairwise_synchrony_generic <- function(df) {
  n_units <- nrow(df)
  if (n_units < 2) return(tibble(synchrony = NA_real_))
  pairs <- combn(n_units, 2)                      # all unique stand pairs
  matches <- df$state[pairs[1, ]] == df$state[pairs[2, ]]
  tibble(synchrony = mean(matches))
}

# Victor's synchrony. Within one year, for one draw: the fraction of stands sitting in whichever state is more common that year (majority share).
majority_fraction_generic <- function(df) {
  state_counts <- table(df$state)
  tibble(consolidation = max(state_counts) / sum(state_counts))
}

# Here a function to summarise a per-draw statistic across draws: median + 90% credible interval. It takes a column of per-draw values (for example: synchrony for every draw × year) and, per year, computes the median and a 90% interval (5th/95th percentile) across draws. This is the step that turns "8000 numbers per year" into "one median line + one uncertainty band per year."
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

# Which draws to use. Using all post-warmup draws is the most correct but can be slow for the pairwise-combination step below as I have many draws/stands.
n_draws_total <- nrow(state_draws)
#draw_idx_use  <- seq_len(n_draws_total)
draw_idx_use <- seq(1, n_draws_total, by = 4)  # used this for the moment to thin and make it faster :)


# PART 1: posterior probability of the mast state, per stand and year
abam_meta <- meta_all %>%
  filter(species == target_species)

abam_state_draws <- state_draws[, abam_meta$row_id]

abam_posterior <- map_dfr(seq_len(nrow(abam_state_draws)), function(d) {
  tibble(
    draw  = d,
    stand = abam_meta$stand,
    year  = abam_meta$year,
    state = abam_state_draws[d, ]
  )
}) %>%
  group_by(stand, year) %>%
  summarise(prob_mast = mean(state == 2), .groups = "drop")

ggplot(abam_posterior, aes(x = year, y = prob_mast * 100, group = stand)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~stand) +
  labs(
    x = NULL, y = "Posterior probability of mast state (%)",
    title = "Predicted masting dynamics of ABAM across stands",
    subtitle = "Posterior probability of latent high-production state"
  ) +
  theme_minimal(base_size = 13)

ggplot(abam_posterior, aes(year, prob_mast * 100, colour = stand)) +
  geom_line(linewidth = 0.8) +
  labs(y = "Probability of mast state (%)", colour = "Stand",
       title = "ABAM predicted masting states across stands") +
  theme_minimal()


# PART 2: pairwise synchrony computed WITHIN each draw first, then summarized across draws 
abam_sync_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~pairwise_synchrony_generic(.x)) %>%
    summarise(synchrony = mean(synchrony, na.rm = TRUE), .groups = "drop") %>%
    mutate(draw = d)
})

abam_sync_summary <- summarise_draws(abam_sync_draws, "synchrony", "year")

ggplot(abam_sync_summary, aes(x = year, y = median * 100)) +
  geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), alpha = 0.2) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 50, linetype = "dashed") +
  labs(x = NULL, y = "ABAM stand synchrony (% of stand pairs)",
       title = "Spatial synchrony of ABAM masting") +
  theme_minimal()


# PART 3: Victor's synchrony. Fraction of stands in the majority state per year  (same per-draw-first logic as synchrony, for direct comparability)
abam_consolidation_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_fraction_generic(.x)) %>%
    mutate(draw = d)
})

abam_consolidation_summary <- summarise_draws(abam_consolidation_draws, "consolidation", "year")

comparison_df <- bind_rows(
  abam_sync_summary %>% mutate(metric = "Pairwise synchrony"),
  abam_consolidation_summary %>% mutate(metric = "majority share")
)

ggplot(comparison_df, aes(x = year, y = median * 100, colour = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  labs(x = NULL, y = "%", colour = NULL, fill = NULL,
       title = "ABAM: pairwise synchrony vs. majority state",
       subtitle = "Both computed per posterior draw, then summarised across draws") +
  theme_minimal()


# PART 4: comparison to the raw observed data
abam_raw <- meta_all %>%
  filter(species == target_species) %>%
  mutate(seed_density = y / area) %>%
  select(stand, year, y, area, seed_density)

posterior_vs_raw <- abam_posterior %>%
  left_join(abam_raw, by = c("stand", "year"))

p_prob <- ggplot(posterior_vs_raw, aes(x = year, y = prob_mast * 100, group = stand)) +
  geom_line(colour = "steelblue") +
  geom_point(size = 1.2, colour = "steelblue") +
  facet_wrap(~stand, nrow = 1) +
  labs(x = NULL, y = "Predicted P(mast) (%)") +
  theme_minimal(base_size = 11)

p_raw <- ggplot(posterior_vs_raw, aes(x = year, y = seed_density, group = stand)) +
  geom_col(fill = "grey40") +
  #scale_y_log10() +
  facet_wrap(~stand, nrow = 1) +
  labs(x = NULL, y = "Observed seed density\n(viable seeds / area)") +
  theme_minimal(base_size = 11)

ggplot(posterior_vs_raw, aes(x = year, y = seed_density, colour = stand)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  labs(x = "Year", y = "Observed seed density\n(viable seeds / area, log scale)",
       colour = "Stand") +
  theme_minimal(base_size = 11)


(p_prob / p_raw) +
  plot_annotation(
    title = "ABAM: posterior predicted mast probability vs. raw observed seed density",
    subtitle = "Top: model prediction  |  Bottom: raw data, no classification applied"
  )

ggplot(posterior_vs_raw, aes(x = seed_density + 1, y = prob_mast * 100)) +
  geom_point(aes(colour = stand), alpha = 0.7) +
  scale_x_log10() +
  labs(x = "Observed seed density + 1 (log scale)", y = "Predicted P(mast) (%)",
       colour = "Stand", title = "Predicted mast probability vs. raw seed density") +
  theme_minimal()

raw_vs_pred_corr <- posterior_vs_raw %>%
  filter(!is.na(seed_density)) %>%
  group_by(stand) %>%
  summarise(
    spearman_r = cor(prob_mast, seed_density, method = "spearman", use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )
print(raw_vs_pred_corr)


n_stands_per_year <- abam_meta %>%
  distinct(stand, year) %>%
  count(year, name = "n_stands")

print(n_stands_per_year)



# PART 5: sensitivity check:  does excluding low-signal stands change the story?
# Stands with a low Spearman correlation between predicted P(mast) and their own raw seed density are ones where the model's state lean heavily on the species-level partial pooling (shared transition rates) rather than on that stand's own counts. This re-runs synchrony and majority state excluding those stands, so I can see whether the overall pattern (dips/peaks over time) depends on them or holds up without them.

# Picked a correlation threshold below which a stand is considered "low signal".
low_r_threshold <- 0.7

low_signal_stands <- raw_vs_pred_corr %>%
  filter(spearman_r < low_r_threshold) %>%
  pull(stand)

message("Excluding these low-signal stands from the sensitivity check: ",
        paste(low_signal_stands, collapse = ", "))

# Re-running the same per-draw synchrony/majority state calculations, just adding a filter to drop the low-signal stands before the pairwise/majority steps.
abam_sync_draws_excl <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species, !stand %in% low_signal_stands) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~pairwise_synchrony_generic(.x)) %>%
    summarise(synchrony = mean(synchrony, na.rm = TRUE), .groups = "drop") %>%
    mutate(draw = d)
})
abam_sync_summary_excl <- summarise_draws(abam_sync_draws_excl, "synchrony", "year")

abam_consolidation_draws_excl <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species, !stand %in% low_signal_stands) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_fraction_generic(.x)) %>%
    mutate(draw = d)
})
abam_consolidation_summary_excl <- summarise_draws(abam_consolidation_draws_excl, "consolidation", "year")

# Overlay full-data vs. excluded-low-signal-stands versions of synchrony, so I can see directly whether dropping shaky stands changes the trend shape.
sensitivity_df <- bind_rows(
  abam_sync_summary       %>% mutate(version = "All stands"),
  abam_sync_summary_excl  %>% mutate(version = "Excluding low-signal stands")
)

ggplot(sensitivity_df, aes(x = year, y = median * 100, colour = version, fill = version)) +
  geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  labs(
    x = NULL, y = "ABAM stand synchrony (% of stand pairs)",
    colour = NULL, fill = NULL,
    title = "Sensitivity check: synchrony with and without low-signal stands",
    subtitle = paste0("Low-signal = Spearman r < ", low_r_threshold,
                      " between predicted P(mast) and raw seed density")
  ) +
  theme_minimal()

# Same comparison for Victor's synchrony
sensitivity_df_consol <- bind_rows(
  abam_consolidation_summary       %>% mutate(version = "All stands"),
  abam_consolidation_summary_excl  %>% mutate(version = "Excluding low-signal stands")
)

ggplot(sensitivity_df_consol, aes(x = year, y = median * 100, colour = version, fill = version)) +
  geom_ribbon(aes(ymin = lo90 * 100, ymax = hi90 * 100), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  labs(
    x = NULL, y = "ABAM (% majority share)",
    colour = NULL, fill = NULL,
    title = "Sensitivity check: majority with and without low-signal stands",
    subtitle = paste0("Low-signal = Spearman r < ", low_r_threshold,
                      " between predicted P(mast) and raw seed density")
  ) +
  theme_minimal()


# IMPORTANT CHECK: excluding stands can shrink an already-small year down to very few remaining stands, which can reintroduce the same issue as the overall low-n problem (Part 2). So I compare stand counts per year for the full data vs. the excluded-stands version before trusting any divergence between the two sensitivity lines.
n_stands_per_year_excl <- abam_meta %>%
  filter(!stand %in% low_signal_stands) %>%
  distinct(stand, year) %>%
  count(year, name = "n_stands_excl")

stand_count_check <- n_stands_per_year %>%
  left_join(n_stands_per_year_excl, by = "year")
print(stand_count_check)
# If n_stands_excl is much smaller than n_stands for a given year (especially in already-low-n years), than I should treat any divergence between "All stands" and"Excluding low-signal stands" for that year with extra caution as it may reflect a small-sample artifact rather than a real difference driven by excluding stands. But it is still ok so all good :)

# Which of the low-signal stands are actually still being sampled in the reduced post-2020 period (2021+)? These are the specific stands responsible for any divergence seen there between the two sensitivity lines.
low_signal_still_sampled_2021plus <- abam_meta %>%
  filter(year >= 2021, stand %in% low_signal_stands) %>%
  distinct(stand) %>%
  pull(stand)

message("Low-signal stands still sampled from 2021 onward: ",
        paste(low_signal_still_sampled_2021plus, collapse = ", "))
#SUNR and TA01

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
draws_idx_use<- seq(1,n_draws_total, by=4)

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
