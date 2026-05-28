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

str(seed_filtered)

# Count number of species per stand per year
species_count_check <- stand_year_all %>%
  group_by(stand, year) %>%
  summarise(n_species = n_distinct(spp), .groups = "drop")

# Look at raw data for one species x stand
stand_year_all %>%
  filter(spp == "ABLA", stand == "AE10") %>%
  ggplot(aes(x = year, y = y / area)) +
  geom_line() + geom_point() +
  labs(title = "Raw seed density - ABLA AE10",
       y = "Seeds per m²")

stand_year_all %>%
  mutate(density = y / area) %>%
  group_by(spp) %>%
  summarise(
    low_approx  = median(density),           # rough low state
    high_approx = quantile(density, 0.9),    # rough high state
    log_low     = log(median(density) + 0.5),  # +0.5 to avoid log(0)
    log_high    = log(quantile(density, 0.9))
  )


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

str(stand_year_all)
unique(stand_year_all$area)

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

unique(years_per_series$spp)

# Mock dataset ------------------------------------------------------------

library(tidyverse)
set.seed(42)

# Known mast years per species — ground truth to check against model output
mast_years <- list(
  ABLA = c(2012, 2015, 2018, 2021),
  PIEN = c(2012, 2015, 2018, 2021),
  PICO = c(2012, 2015, 2018, 2021),
  PIAL = c(2011, 2013, 2016, 2019, 2022),
  LALY = c(2011, 2013, 2016, 2019, 2022),
  ABAM = c(2010, 2014, 2017, 2020, 2023)
)

stand_info <- tibble(
  stand     = c("AB08","AE10","AG05","AM16","AO03","AV02","AV06","AV14","AX15","PP17","TB13"),
  area      = c(1.0204875, 0.8504062, 2.4477830, 1.4272955, 0.9811313,
                0.8110500, 0.5102437, 0.7716937, 2.0518938, 2.2383455, 0.9417750),
  last_year = c(2023, 2023, 2023, 2023, 2018,
                2023, 2016, 2023, 2023, 2023, 2019)
)

stand_year_all <- expand_grid(spp = names(mast_years), stand = stand_info$stand) %>%
  left_join(stand_info, by = "stand") %>%
  rowwise() %>%
  reframe(spp = spp, stand = stand, area = area, year = 2010:last_year) %>%
  rowwise() %>%
  mutate(
    state = as.integer(year %in% mast_years[[spp]]),
    # high state: large counts, low state: small counts
    y     = as.integer(rnbinom(1,
                               mu   = ifelse(state == 1, 150 * area, 5 * area),
                               size = ifelse(state == 1, 5, 3)))
  ) %>%
  ungroup() %>%
  select(spp, stand, year, y, area) %>%
  arrange(spp, stand, year)

str(stand_year_all)

# Fitting Model -----------------------------------------------------------

fit_all <- stan(
    file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_delta.stan",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_allPhi2 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_deltaPhi5.stan",
  data    = stan_data_all,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

fit_allPhi3 <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesPartialPooling_deltaPhi6",
  data    = stan_data_all,
  iter    = 2000, #change based on how much iterations you need
  warmup  = 1000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)

## Checks -------------------------------------------------------------------------
install.packages("tidyverse")
install.packages("tidybayes")
install.packages("ggdist")
library(tidyverse)
library(tidybayes)
library(ggdist)

# Extract log_delta_high — one value per series f, per draw
draws <- fit$draws(variables = "log_delta_high", format = "df")

# Pivot to long format
draws_long <- draws |>
  pivot_longer(
    cols = starts_with("log_delta_high"),
    names_to = "f",
    values_to = "log_delta_high"
  ) |>
  mutate(f = as.integer(str_extract(f, "\\d+")))

# Plot
ggplot(draws_long, aes(x = log_delta_high, y = factor(f))) +
  stat_halfeye(.width = c(0.5, 0.89)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(
    x = "log_delta_high (log difference: high state minus low state)",
    y = "Series (f)",
    title = "Posterior of log_delta_high per series",
    subtitle = "Red line = 0. If intervals overlap 0, states are not separated."
  ) +
  theme_minimal()
# Diagnostics for the state "stickiness" ----------------------------------

# hmm_latent_rng, which samples a single most-likely path (Viterbi-like), it can be prone to "collapsing" to one dominant state if the emission distributions overlap a lot or if transition probabilities are extreme.

#theta
draws <- as.data.frame(fit_all)

# Or visually
mcmc_hist(draws, pars = c("grand_logit_theta1", "grand_logit_theta2"))
inv_logit <- function(x) 1 / (1 + exp(-x))
# Summarise grand-level transition probabilities
draws |>
  mutate(
    theta1_grand = inv_logit(grand_logit_theta1),
    theta2_grand = inv_logit(grand_logit_theta2)
  ) |>
  summarise(
    theta1_mean = mean(theta1_grand),
    theta2_mean = mean(theta2_grand),
    theta1_q5   = quantile(theta1_grand, 0.05),
    theta2_q5   = quantile(theta2_grand, 0.05)
  )

#is the log_delta_high close to 0?
# Extract log_delta_high for all series
delta_cols <- grep("^log_delta_high\\[", names(draws), value = TRUE)

delta_summary <- data.frame(
  series = delta_cols,
  mean   = sapply(delta_cols, \(col) mean(draws[[col]])),
  q5     = sapply(delta_cols, \(col) quantile(draws[[col]], 0.05)),
  q95    = sapply(delta_cols, \(col) quantile(draws[[col]], 0.95))
)

print(delta_summary)

#State proportions
state_cols <- grep("^state\\[", names(draws), value = TRUE)

# Posterior mode per time point
state_mode <- sapply(state_cols, \(col) {
  as.integer(names(which.max(table(draws[[col]]))))
})

# Compare proportion in high state vs rho
cat("Prop in state 2:", mean(state_mode == 2), "\n")
cat("rho[2] mean:",     mean(draws$`rho[2]`), "\n")


# Check series lengths
series_lengths <- end_idxs - start_idxs + 1
summary(series_lengths)

# Then cross with state proportions per series
state_df <- data.frame(t = 1:stan_data_all$N, state = state_mode)

props_per_series <- sapply(1:stan_data_all$F, function(f) {
  idx <- start_idxs[f]:end_idxs[f]
  mean(state_mode[idx] == 2)
})

data.frame(
  series = 1:F,
  length = series_lengths,
  prop_high = props_per_series
) |> arrange(prop_high)


set.seed(42)
n_sim <- 2000
T <- 13  # your median series length

theta1_prior <- plogis(rnorm(n_sim, 0.4, 0.75))
theta2_prior <- plogis(rnorm(n_sim, 0.4, 0.75))

results <- sapply(1:n_sim, function(i) {
  state <- sample(c(1,2), 1)  # random start
  switches <- 0
  time_in_high <- 0
  
  for (t in 2:T) {
    if (state == 1) {
      new_state <- sample(c(1,2), 1, prob = c(theta1_prior[i], 1-theta1_prior[i]))
    } else {
      new_state <- sample(c(1,2), 1, prob = c(1-theta2_prior[i], theta2_prior[i]))
    }
    if (new_state != state) switches <- switches + 1
    if (new_state == 2) time_in_high <- time_in_high + 1
    state <- new_state
  }
  c(switches = switches, prop_high = time_in_high / T)
})

results <- as.data.frame(t(results))

cat("Mean switches per series:", mean(results$switches), "\n")
cat("Prop series with 0 switches:", mean(results$switches == 0), "\n")
cat("Prop series with 1+ switches:", mean(results$switches >= 1), "\n")
cat("Prop series with 3+ switches:", mean(results$switches >= 3), "\n")
cat("Mean prop time in high state:", mean(results$prop_high), "\n")

par(mfrow = c(1,2))
hist(results$switches, main = "Switches per series", xlab = "N switches", col = "steelblue")
hist(results$prop_high, main = "Prop time in high state", xlab = "Proportion", col = "steelblue")


# If label switching is happening, the grand means should show bimodality
# or the low/high assignment should be inconsistent across chains
chain1 <- as.data.frame(fit_all) |> filter(.chain == 1)
chain2 <- as.data.frame(fit_all) |> filter(.chain == 1)

# Compare grand_mean_low across chains
cat("Chain 1 grand_mean_low:", mean(chain1$grand_mean_low), "\n")
cat("Chain 2 grand_mean_low:", mean(chain2$grand_mean_low), "\n")

# If these are very different, you have label switching

# Plot grand_mean_low vs log_delta_high_grand_mean
# Label switching shows as two clusters
plot(draws$grand_mean_low, draws$log_delta_high_grand_mean,
     pch = 20, cex = 0.3, col = rgb(0,0,0,0.2),
     xlab = "grand_mean_low", ylab = "log_delta_high_grand_mean")

# Check if any series has negative log_delta_high (high < low)
delta_cols <- grep("^log_delta_high\\[", names(draws), value = TRUE)
any_negative <- sapply(delta_cols, \(col) mean(draws[[col]] < 0))
print(any_negative)  # should all be ~0

# Check convergence
s <- summary(fit_all)$summary
s[order(s[,"Rhat"], decreasing = TRUE)[1:20], c("mean","Rhat","n_eff")]

# Compare log_alpha_low and log_alpha_high for series 2 vs 15
params_to_check <- c(
  "log_alpha_low[2]",  "log_alpha_high[2]",  "log_delta_high[2]",
  "log_alpha_low[15]", "log_alpha_high[15]", "log_delta_high[15]"
)

s <- summary(fit_all)$summary
s[params_to_check, c("mean", "sd", "2.5%", "97.5%")]

s <- summary(fit_all)$summary
s[c("phi1", "phi2"), c("mean", "sd", "2.5%", "97.5%")]

# Proportion of zeros per series
zero_props <- sapply(1:stan_data_all$F, function(f) {
  idx <- start_idxs[f]:end_idxs[f]
  mean(stan_data_all$y[idx] == 0)
})

par(mfrow,c(1,1))
hist(zero_props, main = "Proportion of zeros per series", 
     xlab = "Prop zeros", col = "steelblue")
cat("Series with >50% zeros:", sum(zero_props > 0.5), "\n")
cat("Series with >30% zeros:", sum(zero_props > 0.3), "\n")


# Which species and stands have the high-zero series?
problem_df <- data.frame(
  series    = 1:stan_data_all$F,
  zero_prop = zero_props,
  species   = sapply(1:stan_data_all$F, function(f) stan_data_all$sp[start_idxs[f]]),
  stand     = stan_data_all$stand_id,
  mean_count = sapply(1:stan_data_all$F, function(f) mean(stan_data_all$y[start_idxs[f]:end_idxs[f]])),
  max_count  = sapply(1:stan_data_all$F, function(f) max(stan_data_all$y[start_idxs[f]:end_idxs[f]]))
) |> arrange(desc(zero_prop))

print(problem_df)

# Are zeros clustered by species or stand?
tapply(zero_props, sapply(1:F, function(f) sp[start_idxs[f]]), mean)
tapply(zero_props, stand_id, mean)

# Diagnostic plots --------------------------------------------------------

check_hmc_diagnostics(fit_allPhi2)

#Checking the divergence issues
# Extract posterior draws (excluding warmup)
draws2 <- as_draws_df(fit_allPhi2)
draws$
# Extract sampler params from rstan
sampler_params <- get_sampler_params(fit_allPhi2, inc_warmup = FALSE)

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


mcmc_pairs(
  draws2,
  pars = c("mu_log_phi", "sigma_log_phi", "log_phi_species[1]", "log_phi_species[2]")
)
draws2$

# Diagnostics -------------------------------------------------------------

source("mcmc_analysis_tools_rstan.R")
source("mcmc_visualization_tools.R")
file.exists("mcmc_analysis_tools_rstan.R")
file.exists("mcmc_visualization_tools.R")
util <- new.env()
#creating a source file to make all the plots. 
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit_allPhi2)
util$check_all_hmc_diagnostics(diagnostics)


# extraire les posterior values
samples <- util$extract_expectand_vals(fit_all_hudle)

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


samples <- util$extract_expectand_vals(fit_allPhi3)
par(mfrow = c(2,1))
#1,5,15,25,30,35,45,57
years_per_series[9, ]
f <- 9
start_id <- start_idxs[f];
end_id <- end_idxs[f];
util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                     display_ylim = c(1,2))
util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
points(x = 1:length(start_id:end_id), y = stan_data_all$y[start_id:end_id], pch = 20)


#new approach for plotting the states
# Extract what you need
y_rep_cols  <- paste0("y_rep[",  1:stan_data_all$N, "]")
state_cols  <- paste0("state[",  1:stan_data_all$N, "]")

draws_mat <- as.matrix(fit_all)

# Pick one series to plot
f        <- 15   # change this
idx      <- start_idxs[f]:end_idxs[f]
T_f      <- length(idx)

y_obs    <- stan_data_all$y[idx]
yrep_f   <- draws_mat[, paste0("y_rep[", idx, "]")]   # n_draws x T
state_f  <- draws_mat[, paste0("state[", idx, "]")]   # n_draws x T

n_draws  <- nrow(yrep_f)

# Compute quantiles separately for draws where each time point is in state 1 vs 2
quant_probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)

q_low  <- matrix(NA, length(quant_probs), T_f)
q_high <- matrix(NA, length(quant_probs), T_f)
q_marg <- matrix(NA, length(quant_probs), T_f)

for (t in 1:T_f) {
  in_low  <- state_f[, t] == 1
  in_high <- state_f[, t] == 2
  
  if (sum(in_low)  > 10) q_low[,  t] <- quantile(yrep_f[in_low,  t], quant_probs)
  if (sum(in_high) > 10) q_high[, t] <- quantile(yrep_f[in_high, t], quant_probs)
  q_marg[, t] <- quantile(yrep_f[, t], quant_probs)
}

# Plot
par(mfrow = c(1,1))
ylim_max <- max(y_obs, q_high[5,], na.rm = TRUE) * 1.2

plot(1:T_f, y_obs, pch = 20, ylim = c(0, ylim_max),
     xlab = "Time", ylab = "y", main = paste("Series", f))

# Marginal envelope (grey)
polygon(c(1:T_f, T_f:1),
        c(q_marg[1,], rev(q_marg[5,])),
        col = rgb(0.5, 0.5, 0.5, 0.2), border = NA)

# Low state envelope (blue)
polygon(c(1:T_f, T_f:1),
        c(q_low[1,], rev(q_low[5,])),
        col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
lines(1:T_f, q_low[3,], col = "blue", lwd = 1.5)

# High state envelope (red)  
polygon(c(1:T_f, T_f:1),
        c(q_high[1,], rev(q_high[5,])),
        col = rgb(0.8, 0.2, 0.2, 0.3), border = NA)
lines(1:T_f, q_high[3,], col = "red", lwd = 1.5)

# Observed data
points(1:T_f, y_obs, pch = 20, cex = 1.2)

# Add state probability as bar at bottom
prob_high_t <- colMeans(state_f == 2)
rect(1:T_f - 0.4, -ylim_max*0.05, 1:T_f + 0.4, -ylim_max*0.05 + prob_high_t * ylim_max*0.08,
     col = rgb(0.8, 0.2, 0.2, 0.6), border = NA)

legend("topright",
       legend = c("Observed", "Low state (90% CI)", "High state (90% CI)", "P(high state)"),
       col    = c("black", "blue", "red", "red"),
       lty    = c(NA, 1, 1, NA),
       pch    = c(20, NA, NA, 15),
       bty    = "n")

# Series 5 - ABAM AO03
f <- 5
idx <- start_idxs[f]:end_idxs[f]
cat("ABAM AO03 counts:\n")
print(stan_data_all$y[idx])

# Series 52 - TSHE AV14  
f <- 52
idx <- start_idxs[f]:end_idxs[f]
cat("TSHE AV14 counts:\n")
print(stan_data_all$y[idx])

plot_series_states <- function(fit, f, stan_data_all, start_idxs, end_idxs, sp_vec) {
  draws_mat <- as.matrix(fit)
  idx <- start_idxs[f]:end_idxs[f]
  T_f <- length(idx)
  y_obs <- stan_data_all$y[idx]
  
  # Soft state probs if you added hmm_hidden_state_prob
  # Otherwise use hard states
  state_cols <- paste0("state[", idx, "]")
  state_mat  <- draws_mat[, state_cols]
  prob_high  <- colMeans(state_mat == 2)
  
  par(mfrow = c(2,1), mar = c(3,4,2,1))
  
  # Top: raw data
  plot(1:T_f, y_obs, type = "b", pch = 20,
       xlab = "", ylab = "Seed count",
       main = paste("Series", f))
  
  # Bottom: P(high state) over time
  plot(1:T_f, prob_high, type = "l", col = "red", lwd = 2,
       ylim = c(0,1), xlab = "Time", ylab = "P(high state)")
  abline(h = 0.5, lty = 2, col = "grey")
  points(1:T_f, prob_high, pch = 20, col = "red")
}

# Compare all three models for series 52
par(mfrow = c(1,3))
plot_series_states(fit_all, 52, stan_data_all, start_idxs, end_idxs, sp)
plot_series_states(fit_all2, 52, stan_data_all, start_idxs, end_idxs, sp)
plot_series_states(fit_all_hudle, 52, stan_data_all, start_idxs, end_idxs, sp)
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




# Prior vs posterior 2 ----------------------------------------------------
n_prior <- 1e4

half_norm <- function(n, sd) abs(rnorm(n, 0, sd))


# 1. SAMPLE TOP-LEVEL PRIORS


prior <- list()

prior$grand_logit_theta1 <- rnorm(n_prior, 1.4, 1)
prior$grand_logit_theta2 <- rnorm(n_prior, 1.4, 1)

alpha_theta1_species <- rnorm(n_prior, 0, prior$sigma_theta1_species)
alpha_theta2_species <- rnorm(n_prior, 0, prior$sigma_theta2_species)
prior$sigma_theta1_species <- half_norm(n_prior, 0.5)
prior$sigma_theta2_species <- half_norm(n_prior, 0.5)

alpha_theta1_stand   <- rnorm(n_prior, 0, prior$sigma_theta1_stand)
alpha_theta2_stand   <- rnorm(n_prior, 0, prior$sigma_theta2_stand)
prior$sigma_theta1_stand   <- half_norm(n_prior, 0.5)
prior$sigma_theta2_stand   <- half_norm(n_prior, 0.5)


# Emissions
prior$grand_mean_low  <- rnorm(n_prior, 2.6, 1.0)

alpha_low_species <- rnorm(n_prior, 0, prior$sigma_low_species)
prior$sigma_low_species  <- half_norm(n_prior, 0.5)

alpha_low_stand   <- rnorm(n_prior, 0, prior$sigma_low_stand)
prior$sigma_low_stand    <- half_norm(n_prior, 0.5)

prior$log_delta_high_grand_mean <- rnorm(n_prior, 3, 1.0)

log_delta_high_species <- rnorm(n_prior, 0, prior$sigma_log_delta_high_species)
prior$sigma_log_delta_high_species <- half_norm(n_prior, 1.0)

log_delta_high_stand   <- rnorm(n_prior, 0, prior$sigma_log_delta_high_stand)
prior$sigma_log_delta_high_stand   <- half_norm(n_prior, 1.0)

# Dispersion (MATCHES STAN)
prior$log_phi_species <- rnorm(n_prior,log(4),0.6)




# 3. BUILD REALIZED PARAMETERS (THIS IS THE KEY)


# Transition probabilities (realized)
theta1_realized <- plogis(
  prior$grand_logit_theta1 +
    alpha_theta1_species +
    alpha_theta1_stand
)

theta2_realized <- plogis(
  prior$grand_logit_theta2 +
    alpha_theta2_species +
    alpha_theta2_stand
)

# Emission means (log scale)
log_alpha_low_realized <-
  prior$grand_mean_low +
  alpha_low_species +
  alpha_low_stand

log_delta_high_realized <-
  prior$log_delta_high_grand_mean +
  log_delta_high_species +
  log_delta_high_stand

# High state mean (important!)
log_alpha_high_realized <-
  log(exp(log_alpha_low_realized) + exp(log_delta_high_realized))


# 4. CREATE FINAL DATAFRAME

prior_df <- data.frame(
  theta1 = theta1_realized,
  theta2 = theta2_realized,
  log_alpha_low = log_alpha_low_realized,
  log_alpha_high = log_alpha_high_realized,
  phi = prior$log_phi_species
)

post <- rstan::extract(fit_allPhi2)

posterior_df <- data.frame(
  theta1 = as.vector(post$theta1),
  theta2 = as.vector(post$theta2),
  log_alpha_low = as.vector(post$log_alpha_low),
  log_alpha_high = as.vector(post$log_alpha_high),
  phi = post$log_phi_species
)

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
    title  = "Prior vs Posterior — Realized Quantities (HMM)",
    x      = "Value",
    y      = "Density",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text      = element_text(size = 8),
    legend.position = "top"
  )

# Prior vs posterior ------------------------------------------------------


n_prior <- 1e4

half_norm <- function(n, sd) abs(rnorm(n, 0, sd))


# 2. SAMPLE FROM PRIORS

prior_list <- list()

# Grand means (emissions)
prior_list$grand_mean_low  <- rnorm(n_prior, 2.6, 1.0)
prior_list$log_delta_high_grand_mean <- rnorm(n_prior, 5.7, 1.0)

# Transition logits
prior_list$grand_logit_theta1 <- rnorm(n_prior, 0, 0.5)
prior_list$grand_logit_theta2 <- rnorm(n_prior, 0, 0.5)

# Variance parameters (half-normal)
prior_list$sigma_theta1_species <- half_norm(n_prior, 0.3)
prior_list$sigma_theta2_species <- half_norm(n_prior, 0.3)
prior_list$sigma_theta1_stand   <- half_norm(n_prior, 0.3)
prior_list$sigma_theta2_stand   <- half_norm(n_prior, 0.3)

prior_list$sigma_low_species  <- half_norm(n_prior, 0.5)
prior_list$sigma_log_delta_high_species <- half_norm(n_prior, 0.5)
prior_list$sigma_low_stand    <- half_norm(n_prior, 0.5)
prior_list$sigma_log_delta_high_stand   <- half_norm(n_prior, 0.5)

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
post$
post_list$grand_mean_low  <- post$grand_mean_low
post_list$log_delta_high_grand_mean <- post$log_delta_high_grand_mean

# Transition logits
post_list$grand_logit_theta1 <- post$grand_logit_theta1
post_list$grand_logit_theta2 <- post$grand_logit_theta2

# Variance parameters
post_list$sigma_theta1_species <- post$sigma_theta1_species
post_list$sigma_theta2_species <- post$sigma_theta2_species
post_list$sigma_theta1_stand   <- post$sigma_theta1_stand
post_list$sigma_theta2_stand   <- post$sigma_theta2_stand

post_list$sigma_low_species  <- post$sigma_low_species
post_list$sigma_log_delta_high_species <- post$sigma_log_delta_high_species
post_list$sigma_low_stand    <- post$sigma_low_stand
post_list$sigma_log_delta_high_stand   <- post$sigma_log_delta_high_stand

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



# Synchrony ---------------------------------------------------------------------

library(tidyverse)
library(rstan)
library(patchwork)

stand_elev <- data.frame(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13",
            "AO03","AG05","AV06","AX15","AB08","PP17",
            "AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600,668,700,850,1450,850,
                900,950,1060,1090,1100,1150,
                1150,1200,1450,1600,1700,1800)
)

sync_draws <- rstan::extract(fit_allPhi2, pars = "stand_sync")$stand_sync

stand_lookup <- years_per_series %>%
  distinct(stand, stand_id) %>%
  arrange(stand_id) %>%
  left_join(stand_elev, by = "stand")

sync_summary <- as.data.frame(sync_draws) %>%
  set_names(paste0("stand_", 1:N_stands)) %>%
  pivot_longer(everything(), names_to = "stand_id_label", values_to = "sync") %>%
  mutate(stand_id = as.integer(str_extract(stand_id_label, "\\d+"))) %>%
  left_join(stand_lookup, by = "stand_id") %>%
  group_by(stand, stand_id, elevation) %>%
  summarise(
    mean_sync  = mean(sync),
    lower_95   = quantile(sync, 0.025),
    upper_95   = quantile(sync, 0.975),
    lower_50   = quantile(sync, 0.25),
    upper_50   = quantile(sync, 0.75),
    .groups    = "drop"
  ) %>%
  arrange(elevation)

print(sync_summary)

dev.off()
p_dot <- ggplot(sync_summary,
                aes(x = reorder(stand, mean_sync), y = mean_sync)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 colour = "blue", linewidth = 0.6) +
  geom_linerange(aes(ymin = lower_50, ymax = upper_50),
                 colour = "lightblue", linewidth = 1.5) +
  geom_point(colour = "red", size = 2.5) +
  coord_flip() +
  scale_y_continuous(limits = c(0.5, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Stand-level Masting Synchrony",
    subtitle = "Proportion of species in majority state, averaged across years\nThick = 50% CI | Thin = 95% CI",
    x        = "Stand",
    y        = "Synchrony index"
  ) +
  theme_bw(base_size = 11)

p_density <- as.data.frame(sync_draws) %>%
  set_names(stand_lookup$stand) %>%
  pivot_longer(everything(), names_to = "stand", values_to = "sync") %>%
  ggplot(aes(x = sync, fill = stand, colour = stand)) +
  geom_density(alpha = 0.15, linewidth = 0.6) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(limits = c(0.5, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Posterior Distributions of Synchrony by Stand",
    x        = "Synchrony index",
    y        = "Density",
    fill     = NULL, colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")


combined <- p_dot | p_density

p_elev_rank <- ggplot(sync_summary,
                      aes(x = reorder(stand, elevation), y = mean_sync)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 colour = "#4e9af1", linewidth = 0.6) +
  geom_linerange(aes(ymin = lower_50, ymax = upper_50),
                 colour = "#1a5fa8", linewidth = 1.5) +
  geom_point(aes(colour = elevation), size = 3) +
  scale_colour_gradient(low = "#a8dadc", high = "#e63946",
                        name = "Elevation (m)") +
  coord_flip() +
  scale_y_continuous(limits = c(0.5, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Stand-level Synchrony Ranked by Elevation",
    subtitle = "Thick = 50% CI | Thin = 95% CI",
    x        = "Stand (low → high elevation)",
    y        = "Synchrony index"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

p_scatter <- ggplot(sync_summary,
                    aes(x = elevation, y = mean_sync)) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#e63946", fill = "#e63946", alpha = 0.15,
              linewidth = 0.8) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 colour = "grey60", linewidth = 0.5) +
  geom_linerange(aes(ymin = lower_50, ymax = upper_50),
                 colour = "#1a5fa8", linewidth = 1.2) +
  geom_point(aes(fill = elevation), shape = 21,
             size = 3, colour = "white", stroke = 0.5) +
  geom_text(aes(label = stand), size = 2.8,
            vjust = -1, hjust = 0.5, colour = "grey30") +
  scale_fill_gradient(low = "#a8dadc", high = "#e63946",
                      name = "Elevation (m)") +
  scale_y_continuous(limits = c(0.5, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Synchrony vs Elevation",
    subtitle = "Each point = one stand | Error bars = 50% and 95% posterior CI\nLine = linear trend",
    x        = "Elevation (m)",
    y        = "Synchrony index"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

draws_long <- as.data.frame(sync_draws) %>%
  set_names(stand_lookup$stand) %>%
  pivot_longer(everything(), names_to = "stand", values_to = "sync") %>%
  left_join(stand_elev, by = "stand")

p_ribbon <- draws_long %>%
  group_by(stand, elevation) %>%
  summarise(
    q05  = quantile(sync, 0.05),
    q25  = quantile(sync, 0.25),
    q50  = quantile(sync, 0.50),
    q75  = quantile(sync, 0.75),
    q95  = quantile(sync, 0.95),
    .groups = "drop"
  ) %>%
  arrange(elevation) %>%
  ggplot(aes(x = elevation)) +
  geom_ribbon(aes(ymin = q05, ymax = q95),
              fill = "#4e9af1", alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill = "#1a5fa8", alpha = 0.35) +
  geom_line(aes(y = q50), colour = "#1a5fa8", linewidth = 0.9) +
  geom_point(aes(y = q50, fill = elevation), shape = 21,
             size = 2.5, colour = "white") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  scale_fill_gradient(low = "#a8dadc", high = "#e63946",
                      name = "Elevation (m)") +
  scale_y_continuous(limits = c(0.5, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Synchrony Along the Elevational Gradient",
    subtitle = "Dark ribbon = 50% CI | Light ribbon = 90% CI | Line = posterior median",
    x        = "Elevation (m)",
    y        = "Synchrony index"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

combined <- (p_elev_rank | p_scatter) / p_ribbon +
  plot_annotation(
    title = "Masting Synchrony Across the Elevational Gradient",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )



# Synchrony like Victor ---------------------------------------------------

library(tidyverse)
library(patchwork)


sync_draws <- rstan::extract(fit_allPhi2, pars = "stand_sync")$stand_sync
# [n_draws x N_stands]

stand_lookup <- years_per_series %>%
  distinct(stand, stand_id) %>%
  arrange(stand_id)

# We also need to know which years each stand covers
stand_years_lookup <- stand_year_all %>%
  group_by(stand) %>%
  summarise(years = list(sort(unique(year))), .groups = "drop")


# 2. WITHIN-STAND synchrony per year
# For each draw, we need stand_sync broken down by year, not just averaged.
# This requires re-extracting from the state draws.

state_draws <- rstan::extract(fit_all, pars = "state")$state

# Compute per-stand per-year synchrony for each draw
base_df <- stand_year_all %>%
  select(spp, stand, year) %>%
  mutate(row_idx = row_number())

set.seed(42)
n_draws     <- min(500, nrow(state_draws))
draw_sample <- sample(1:nrow(state_draws), n_draws)

within_per_draw <- map_dfr(seq_along(draw_sample), function(i) {
  base_df %>%
    mutate(state = state_draws[draw_sample[i], ]) %>%
    pivot_wider(id_cols = c(stand, year),
                names_from = spp, values_from = state) %>%
    rowwise() %>%
    mutate(
      vals       = list(c_across(-c(stand, year))),
      count_high = sum(unlist(vals) == 2, na.rm = TRUE),
      count_low  = sum(unlist(vals) == 1, na.rm = TRUE),
      total      = count_high + count_low,
      sync       = max(count_high, count_low) / total
    ) %>%
    ungroup() %>%
    select(stand, year, sync) %>%
    mutate(draw = i)
}, .progress = TRUE)


# 3. BETWEEN-STAND synchrony per year
# For each year and draw: proportion of stands where majority species agree
# i.e. landscape-level — are stands in the same state as each other?
# Computed as: for each year, what is the mean stand_sync across stands?
# (landscape sync = mean of within-stand syncs across stands that year)

between_per_draw <- within_per_draw %>%
  group_by(year, draw) %>%
  summarise(sync = mean(sync), .groups = "drop")  # landscape mean

# 4. POSTERIOR SUMMARIES

within_summary <- within_per_draw %>%
  group_by(stand, year) %>%
  summarise(
    median = median(sync),
    q25    = quantile(sync, 0.25),
    q75    = quantile(sync, 0.75),
    q05    = quantile(sync, 0.05),
    q95    = quantile(sync, 0.95),
    .groups = "drop"
  ) %>%
  # average within-stand across stands per year for the "within" line
  group_by(year) %>%
  summarise(
    median = median(median),
    q25    = median(q25),
    q75    = median(q75),
    q05    = median(q05),
    q95    = median(q95),
    .groups = "drop"
  ) %>%
  mutate(type = "Within stand")

between_summary <- between_per_draw %>%
  group_by(year) %>%
  summarise(
    median = median(sync),
    q25    = quantile(sync, 0.25),
    q75    = quantile(sync, 0.75),
    q05    = quantile(sync, 0.05),
    q95    = quantile(sync, 0.95),
    .groups = "drop"
  ) %>%
  mutate(type = "Between stands")

plot_df <- bind_rows(within_summary, between_summary)


# 5. PLOT

pal <- c("Within stand"   = "#7b4f8e",   # purple like reference
         "Between stands" = "#4a7fa5")    # blue like reference

p <- ggplot(plot_df, aes(x = year, colour = type, fill = type)) +
  # 90% CI box (light)
  geom_rect(aes(xmin = year - 0.3, xmax = year + 0.3,
                ymin = q05, ymax = q95),
            alpha = 0.2, colour = NA) +
  # 50% CI box (dark)
  geom_rect(aes(xmin = year - 0.3, xmax = year + 0.3,
                ymin = q25, ymax = q75),
            alpha = 0.4, colour = NA) +
  # median line inside box
  geom_segment(aes(x = year - 0.3, xend = year + 0.3,
                   y = median, yend = median),
               linewidth = 1.2) +
  # dashed line connecting medians
  geom_line(aes(y = median), linetype = "dashed",
            linewidth = 0.6, alpha = 0.8) +
  # reference line
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  annotate("text", x = min(plot_df$year), y = 0.505,
           label = "Minimum synchrony", colour = "grey60",
           hjust = 0, vjust = 0, size = 3) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal,   name = NULL) +
  scale_y_continuous(
    limits = c(0.5, 1),
    labels = scales::percent_format(accuracy = 1),
    name   = "Synchrony (% species in majority state)"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
  labs(
    title    = "Masting Synchrony Through Time",
    subtitle = "Boxes = 50% CI | whiskers = 90% CI | dashed = median trend",
    x        = "Year"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "top",
    legend.key.width = unit(1.5, "cm"),
    panel.grid.minor = element_blank()
  )

ggsave("synchrony_through_time.png", p,
       width = 12, height = 6, dpi = 150, bg = "white")
