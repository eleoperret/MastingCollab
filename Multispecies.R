##Code for multispecies HMM model 
# Stan code: 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

####Doesn't include PARA and SUNR as I do not have a map of the species. 
####Issues with ABLA so I will remove it. 

#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)


#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

seed_data<-read.csv("SeedData_all.csv")


str(seed_data)

# 1) Remove species I'm not interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","CANO","PSME","TSHE","THPL"))

# 2) Removing invalid stands 
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)

lol <- seed_filtered %>%
  group_by(stand, year) %>%
  summarise(n_traps = n_distinct(trapno), .groups = "drop")
str(lol)

unique(seed_filtered$spp)
unique(seed_filtered$stand)

str(seed_filtered)

# # Load stand elevations from saved RDS
# stand_elev <- readRDS("data/stand_elevation_table.rds")#from General.Data.R
# species_data <- left_join(seed_filtered, stand_elev, by = "stand")


#Creating the data for stan

##Add the removal of year 2009
stand_year_df <- seed_filtered %>%
  group_by(spp, stand, year) %>%
  filter(!year=="2009")%>% #removing the year 2009
  summarise(
    y = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  )

stand_year_df <- stand_year_df %>%
  arrange(spp, stand, year)

stand_year_df$species_id <- as.numeric(as.factor(stand_year_df$spp))

S <- length(unique(stand_year_df$spp))

years_per_series <- stand_year_df %>%
  group_by(spp, stand) %>%
  summarise(T_i = n(), .groups = "drop")

G <- nrow(years_per_series)
T_i <- years_per_series$T_i

start_idxs<-c()
end_idxs<-c()

id<-1

for (g in 1:G) {
  start_idxs <- c(start_idxs, id)
  id <- id + T_i[g] - 1
  end_idxs <- c(end_idxs, id)
  id <- id + 1
}


y <- stand_year_df$y
area <- stand_year_df$area
sp <- stand_year_df$species_id

N <- length(y)

stan_data <- list(
  N = N,
  F = G,     # number of species x stand combination
  S = S,     # number of species
  sp = sp,
  start_idxs = start_idxs,
  end_idxs = end_idxs,
  y = y,
  area = area
)


init_fn <- function() {
  list(
    rho              = c(0.5, 0.5),
    theta1           = rep(0.75, stan_data$S),
    theta2           = rep(0.50, stan_data$S),
    log_means        = matrix(
      c(rep(1.5, stan_data$S),   # low  state
        rep(5.0, stan_data$S)),  # high state
      nrow = stan_data$S,
      ncol = 2
    ),
    log_phi1         = rep(log(6.5), stan_data$S),
    log_phi2         = rep(log(6.5), stan_data$S),
    sigma            = 0.3,
    stand_effect_raw = rep(0, stan_data$F)
  )
}

options(mc.cores = parallel::detectCores())

fit <- stan(
  file = "Stan_code/Species_Stan_Model/MultispeciesTransitionSpecies.stan",  # your Stan file
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123, 
  init = init_fn, 
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)



# 1. Find which parameters have NA or high Rhat
rhat_vals <- rhat(fit)
rhat_vals[is.na(rhat_vals)]          # which are NA
rhat_vals[rhat_vals > 1.05 & !is.na(rhat_vals)]  # which are bad

# 2. Find which parameters have low ESS
neff_vals <- neff_ratio(fit)
neff_vals[neff_vals < 0.1 & !is.na(neff_vals)]   # worst offenders

# 3. Full parameter summary — focus on log_means
print(fit, pars = c("theta1", "theta2", "log_means",
                    "log_phi1", "log_phi2", "sigma", "rho"))

# 4. Where are the divergences occurring?
library(bayesplot)
posterior <- as.array(fit)

# Pairs plot for the most suspicious parameters
mcmc_pairs(posterior, 
           pars = c("log_means[3,1]", "log_means[3,2]", 
                    "theta2[3]", "sigma"),
           np = nuts_params(fit))   # highlights divergences in red

# 5. Trace plots for previously problematic parameters
mcmc_trace(posterior, 
           pars = c("log_means[1,1]", "log_means[1,2]",
                    "log_means[3,1]", "log_means[3,2]",
                    "log_means[5,1]", "log_means[5,2]"))

# 6. Check sigma specifically — collapsing to zero causes NA Rhat
mcmc_trace(posterior, pars = c("sigma"))
print(fit, pars = "sigma")


print(fit, pars = c("sigma", "stand_effect_raw[1]", 
                    "stand_effect_raw[2]", "stand_effect_raw[3]"))


# Are THPL means (species 4) close together?
print(fit, pars = c("log_means[4,1]", "log_means[4,2]"))
# If 97.5% of [4,1] overlaps with 2.5% of [4,2], the states are not well separated



# Check which stand-species combinations have year gaps
gap_check <- stand_year_df %>%
  arrange(spp, stand, year) %>%
  group_by(spp, stand) %>%
  summarise(
    years      = list(year),
    has_gaps   = any(diff(year) > 1),
    gap_years  = paste(year[which(diff(year) > 1)], 
                       year[which(diff(year) > 1) + 1], 
                       sep = "->", collapse = ", "),
    .groups    = "drop"
  ) %>%
  filter(has_gaps)

print(gap_check)


stand_year_df %>%
  group_by(spp) %>%
  summarise(
    n_stands    = n_distinct(stand),
    year_min    = min(year),
    year_max    = max(year),
    n_obs       = n(),
    n_zeros     = sum(y == 0),
    pct_zeros   = round(100 * mean(y == 0), 1),
    mean_y_nonzero = round(mean(y[y > 0]), 1)
  )


library(ggplot2)

stand_year_df %>%
  filter(spp %in% c("PSME", "TSHE")) %>%
  ggplot(aes(x = log1p(y))) +
  geom_histogram(bins = 30) +
  facet_wrap(~spp, scales = "free") +
  labs(x = "log(seeds + 1)", title = "Are distributions bimodal?")

# Compare all species
stand_year_df %>%
  ggplot(aes(x = log1p(y))) +
  geom_histogram(bins = 30) +
  facet_wrap(~spp, scales = "free") +
  labs(x = "log(seeds + 1)")





# new part: ---------------------------------------------------------------

library(dplyr)
library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())

# -------------------------------------------------------
# 1. Data preparation — remove PSME
# -------------------------------------------------------

stand_year_4sp <- seed_filtered %>%
  filter(year != 2009) %>%
  filter(spp != "PSME") %>%                # remove PSME
  group_by(spp, stand, year) %>%
  summarise(
    y    = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(spp, stand, year)

# Species index — should now be 1:4 for ABAM, CANO, THPL, TSHE
stand_year_4sp$species_id <- as.numeric(as.factor(stand_year_4sp$spp))
S <- length(unique(stand_year_4sp$spp))

# Confirm species order — important for matching priors to species
cat("Species order (check this matches your Stan priors):\n")
print(levels(as.factor(stand_year_4sp$spp)))
# Should print: ABAM CANO THPL TSHE
# s=1: ABAM, s=2: CANO, s=3: THPL, s=4: TSHE

# Build segment indices
years_per_series <- stand_year_4sp %>%
  group_by(spp, stand) %>%
  summarise(T_i = n(), .groups = "drop")

G          <- nrow(years_per_series)
T_i        <- years_per_series$T_i
start_idxs <- cumsum(c(1, T_i[-G]))
end_idxs   <- cumsum(T_i)

# Sanity checks
cat("Total observations N       :", nrow(stand_year_4sp), "\n")
cat("Number of segments F       :", G, "\n")
cat("Number of species S        :", S, "\n")
cat("Indices correctly computed :", all(end_idxs - start_idxs + 1 == T_i), "\n")
cat("Final end_idx matches N    :", tail(end_idxs, 1) == nrow(stand_year_4sp), "\n")

stan_data_4sp <- list(
  N          = nrow(stand_year_4sp),
  F          = G,
  S          = S,
  sp         = stand_year_4sp$species_id,
  start_idxs = start_idxs,
  end_idxs   = end_idxs,
  y          = stand_year_4sp$y,
  area       = stand_year_4sp$area
)

# -------------------------------------------------------
# 2. Initial values
# -------------------------------------------------------

init_fn_4sp <- function() {
  list(
    rho              = c(0.5, 0.5),
    theta1           = rep(0.75, stan_data_4sp$S),
    theta2           = rep(0.50, stan_data_4sp$S),
    # Rows = species, columns = [low state, high state]
    # s=1 ABAM, s=2 CANO, s=3 THPL, s=4 TSHE
    log_means        = matrix(
      c(0.5,  0.5,  2.0,  4.5,   # low  state per species
        3.0,  3.5,  5.0,  5.5),  # high state per species
      nrow = stan_data_4sp$S,
      ncol = 2
    ),
    log_phi1         = rep(log(6.5), stan_data_4sp$S),
    log_phi2         = rep(log(6.5), stan_data_4sp$S),
    sigma            = 0.5,
    stand_effect_raw = rep(0, stan_data_4sp$F)
  )
}

# -------------------------------------------------------
# 3. Fit the model
# -------------------------------------------------------

fit_4sp <- stan(
  file    = "Stan_code/Species_Stan_Model/2ndMultispeciesTransitionSpecies.stan",
  data    = stan_data_4sp,
  iter    = 4000,
  warmup  = 2000,
  chains  = 4,
  seed    = 123,
  init    = init_fn_4sp,
  control = list(
    adapt_delta   = 0.95,
    max_treedepth = 12
  )
)

# -------------------------------------------------------
# 4. Convergence diagnostics
# -------------------------------------------------------

# Key parameters summary
print(fit_4sp, pars = c("theta1", "theta2", "log_means",
                        "log_phi1", "log_phi2", "sigma", "rho"))

# Check Rhat — all should be < 1.05
rhat_vals <- rhat(fit_4sp)
cat("\nParameters with Rhat > 1.05:\n")
print(rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.05])

cat("\nParameters with NA Rhat (check if only state[]):\n")
print(names(rhat_vals[is.na(rhat_vals)]))

# Check ESS — all should be > 400 ideally
neff_vals <- neff_ratio(fit_4sp)
cat("\nParameters with ESS ratio < 0.1:\n")
print(neff_vals[!is.na(neff_vals) & neff_vals < 0.1])

# Trace plots for all key parameters
posterior_4sp <- as.array(fit_4sp)

mcmc_trace(posterior_4sp,
           pars = c("theta1[1]", "theta1[2]", "theta1[3]", "theta1[4]"))

mcmc_trace(posterior_4sp,
           pars = c("log_means[1,1]", "log_means[1,2]",
                    "log_means[2,1]", "log_means[2,2]",
                    "log_means[3,1]", "log_means[3,2]",
                    "log_means[4,1]", "log_means[4,2]"))

mcmc_trace(posterior_4sp, pars = c("sigma"))

# Rhat and ESS overview plots
mcmc_rhat(rhat_vals[!is.na(rhat_vals)])
mcmc_neff(neff_vals[!is.na(neff_vals)])

# -------------------------------------------------------
# 5. Posterior predictive check (if convergence looks good)
# -------------------------------------------------------

y_rep <- extract(fit_4sp, "y_rep")$y_rep

# Overall distribution
ppc_dens_overlay(
  y     = stan_data_4sp$y,
  yrep  = y_rep[1:100, ],
  trans = "log1p"           # log scale for readability
)

# Per species — need species index for each observation
ppc_dens_overlay_grouped(
  y      = stan_data_4sp$y,
  yrep   = y_rep[1:100, ],
  group  = stand_year_4sp$spp,
  trans  = "log1p"
)
