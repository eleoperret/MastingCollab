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


#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)

options(mc.cores = parallel::detectCores())


#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

seed_data<-read.csv("SeedData_all.csv")

# 1) Keep species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","CANO","PSME","TSHE","THPL"))
#I did not include Abla because it is present only at two stands and the other ones are the PARA and SUNR which I still do not have the data for


# 2) Removing invalid stands 
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

# # Overlook at the data
# stand_year_all<- seed_filtered %>%
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


# Fitting Model -----------------------------------------------------------

fit_all <- stan(
  file    = "Stan_code/Species_Stan_Model/MultispeciesB5.stan",
  data    = stan_data_all,
  iter    = 3000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)





###Today's work
# Step 1 — find the worst parameters
all_summary <- summary(fit_all)$summary
all_summary[order(all_summary[,"Rhat"], decreasing=TRUE)[1:20], 
            c("mean","sd","n_eff","Rhat")]

# Step 2 — check the key structural parameters
print(fit_all, pars = c("sigma_low", "sigma_high",
                        "mu_log_phi1", "mu_log_phi2",
                        "sigma_log_phi1", "sigma_log_phi2",
                        "log_phi1", "log_phi2"))

# Step 3 — pairs plot to see where the geometry is bad
pairs(fit_all, pars = c("sigma_low", "sigma_high",
                        "mu_log_phi1", "sigma_log_phi1",
                        "mu_log_phi2", "sigma_log_phi2"))

# Are any phi values near zero or negative?
print(fit_all, pars = c("log_phi1", "log_phi2"))
# If log_phi1 is very negative (e.g. -5), exp(-5) = 0.007 
# which can cause numerical issues in neg_binomial_2_log_rng


library(bayesplot)
library(bayesplot)
library(patchwork)

posterior_array <- as.array(fit_all)
divergences <- nuts_params(fit_all)

# Fixed style — remove unsupported arguments
np_style <- pairs_style_np(
  div_color = "red",
  div_shape = 16,
  div_size  = 1.5
)

# ── 1. Hyperparameters ────────────────────────────────────────────────────────
p_hyper <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = c("sigma_low", "sigma_high",
               "mu_log_phi1", "sigma_log_phi1",
               "mu_log_phi2", "sigma_log_phi2"),
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

# ── 2. Species-level means (log scale) ───────────────────────────────────────
p_means <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = c("log_means[1,1]", "log_means[1,2]",
               "log_means[2,1]", "log_means[2,2]",
               "log_means[3,1]", "log_means[3,2]",
               "log_means[4,1]", "log_means[4,2]",
               "log_means[5,1]", "log_means[5,2]"),
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

# ── 3. Overdispersion per species ─────────────────────────────────────────────
p_phi <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = c("log_phi1[1]", "log_phi1[2]", "log_phi1[3]",
               "log_phi1[4]", "log_phi1[5]",
               "log_phi2[1]", "log_phi2[2]", "log_phi2[3]",
               "log_phi2[4]", "log_phi2[5]"),
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

# ── 4. Transition probabilities — sample of stands ───────────────────────────
# Select a few representative stands per species
# Adjust indices based on your years_per_series ordering
sample_stands_theta <- c(
  paste0("theta1[", c(1, 2, 3, 14, 15, 25, 26, 37, 38, 50, 51, 56), "]")
)
sample_stands_theta2 <- c(
  paste0("theta2[", c(1, 2, 3, 14, 15, 25, 26, 37, 38, 50, 51, 56), "]")
)

p_theta1 <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = sample_stands_theta,
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

p_theta2 <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = sample_stands_theta2,
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

# ── 5. Stand-level log_alpha — sample across species ─────────────────────────
# Pick ~3 stands per species
sample_low  <- paste0("log_alpha_low[",  c(1,2,3, 14,15,16, 25,26,27, 37,38,39, 50,51,52), "]")
sample_high <- paste0("log_alpha_high[", c(1,2,3, 14,15,16, 25,26,27, 37,38,39, 50,51,52), "]")

p_alpha_low <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = sample_low,
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

p_alpha_high <- mcmc_pairs(
  posterior_array,
  np       = divergences,
  np_style = np_style,
  pars     = sample_high,
  diag_fun     = "hist",
  off_diag_fun = "scatter"
)

# ── 6. Traceplot — all key parameters ────────────────────────────────────────
p_trace_hyper <- mcmc_trace(
  posterior_array,
  pars     = c("sigma_low", "sigma_high",
               "mu_log_phi1", "mu_log_phi2",
               "sigma_log_phi1", "sigma_log_phi2"),
  np       = divergences,
  np_style = trace_style_np(div_color = "red", div_size = 0.5)
) + labs(title = "Traceplots — hyperparameters")

p_trace_means <- mcmc_trace(
  posterior_array,
  pars     = c("log_means[1,1]", "log_means[1,2]",
               "log_means[2,1]", "log_means[2,2]",
               "log_means[3,1]", "log_means[3,2]",
               "log_means[4,1]", "log_means[4,2]",
               "log_means[5,1]", "log_means[5,2]"),
  np       = divergences,
  np_style = trace_style_np(div_color = "red", div_size = 0.5)
) + labs(title = "Traceplots — species means")

p_trace_phi <- mcmc_trace(
  posterior_array,
  pars     = c("log_phi1[1]", "log_phi1[2]", "log_phi1[3]",
               "log_phi1[4]", "log_phi1[5]",
               "log_phi2[1]", "log_phi2[2]", "log_phi2[3]",
               "log_phi2[4]", "log_phi2[5]"),
  np       = divergences,
  np_style = trace_style_np(div_color = "red", div_size = 0.5)
) + labs(title = "Traceplots — overdispersion")

# ── 7. Save everything to PDF ─────────────────────────────────────────────────
pdf("HMM_diagnostics.pdf", width = 14, height = 12)

print(p_hyper);       title("Pairs — hyperparameters (red = divergent)")
print(p_means);       title("Pairs — species log_means")
print(p_phi);         title("Pairs — overdispersion log_phi")
print(p_theta1);      title("Pairs — theta1 sample of stands")
print(p_theta2);      title("Pairs — theta2 sample of stands")
print(p_alpha_low);   title("Pairs — log_alpha_low sample of stands")
print(p_alpha_high);  title("Pairs — log_alpha_high sample of stands")
print(p_trace_hyper)
print(p_trace_means)
print(p_trace_phi)

dev.off()

cat("Saved to HMM_diagnostics.pdf\n")





###














stand_year_all %>%
  group_by(spp) %>%
  summarise(n_stands = n_distinct(stand), 
            avg_years = mean(table(stand)))

# 1. Which parameters have the worst R-hat / ESS?
print(fit_all, pars = c("theta1","theta2","log_means","mu_log_phi1", "mu_log_phi2", 
                        "sigma_log_phi1", "sigma_log_phi2",
                        "rho_conc", "sigma_low", "sigma_high"))

# 2. Where are the divergences happening?
pairs(fit_all, pars = c("mu_log_phi1", "sigma_log_phi1", 
                        "mu_log_phi2", "sigma_log_phi2"))

pairs(fit_all, pars = c("rho_conc", "sigma_low", "sigma_high"))

# 3. Check the chain that got stuck
traceplot(fit_all, pars = c("mu_log_phi1", "sigma_log_phi1", 
                            "rho_conc"), inc_warmup = TRUE)


stand_year_all %>%
  filter(species_id == 1) %>%
  summarise(
    n_stands    = n_distinct(stand),
    n_years     = n_distinct(year),
    n_obs       = n(),
    total_seeds = sum(y),
    zeros       = mean(y == 0)
  )

stand_year_all %>%
  filter(species_id == 1) %>%
  group_by(stand) %>%
  summarise(max_y = max(y), mean_y = mean(y), zeros = mean(y==0)) %>%
  arrange(zeros)

all_summary <- summary(fit_all)$summary
all_summary[order(all_summary[,"Rhat"], decreasing=TRUE)[1:20], 
            c("mean","sd","n_eff","Rhat")]

# Find the worst parameters across ALL parameters, not just these ones
all_summary <- summary(fit_all)$summary
all_summary[order(all_summary[,"Rhat"], decreasing=TRUE)[1:20], 
            c("mean","sd","n_eff","Rhat")]

# Which stand is series 45?
years_per_series[45, ]

# How many observations does it span?
cat("Observations:", start_idxs[45], "to", end_idxs[45], "\n")
# This should show ~589:599 — confirming it's stand 45

# Look at its data
stand_45_info <- years_per_series[45, ]
stand_year_all %>%
  filter(spp == stand_45_info$spp, 
         stand == stand_45_info$stand)

# 4. Diagnostics --------------------------------------------------------

#New style

#Fit of my parameters
print(fit_all, pars = c("theta1", "theta2", "log_means",
                        "log_phi1", "log_phi2", "sigma_low","sigma_high", "rho"))

#Caterpillar plots for some of my parameters. 
traceplot(fit_all, pars = c("log_means[1,1]", "log_means[1,2]", "sigma"))

#Checking Rhat values 
rhat_vals <- rhat(fit_all)
cat("\nRhat > 1.05 (excluding state and log_omega):\n")
bad_rhat <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.05]
print(bad_rhat[!grepl("state|log_omega|y_rep", names(bad_rhat))])

#Checking nESS ratio (efficiency of sampler -- values <1% bad as it means that draws were redundant or hiigly correlated)
cat("\nESS ratio < 0.1 (excluding state and log_omega):\n")
neff_vals <- neff_ratio(fit_all)
bad_neff  <- neff_vals[!is.na(neff_vals) & neff_vals < 0.1]
print(bad_neff[!grepl("state|log_omega|y_rep", names(bad_neff))])

# Trace plots for seed production
posterior_all <- as.array(fit_all)
mcmc_trace(posterior_all,
           pars = c("log_means[1,1]", "log_means[1,2]",
                    "log_means[2,1]", "log_means[2,2]",
                    "log_means[3,1]", "log_means[3,2]"))
mcmc_trace(posterior_all,
           pars = c("log_means[4,1]", "log_means[4,2]",
                    "log_means[5,1]", "log_means[5,2]"))
mcmc_trace(posterior_all,
           pars = c("sigma_log_phi1", "sigma_log_phi2","rho_conc[1]",
                    "rho_conc[2]"))


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

# -------------------------------------------------------
# 7. Stan's util plots if you have Mike's tools loaded
# -------------------------------------------------------

if (exists("util")) {
  
  samples <- util$extract_expectand_vals(fit_all)
  
  # HMC diagnostics
  diagnostics <- util$extract_hmc_diagnostics(fit_all)
  util$check_all_hmc_diagnostics(diagnostics)
  
  # Parameter diagnostics — using correct parameter names
  base_samples <- util$filter_expectands(
    samples,
    c("rho", "theta1", "theta2",
      "log_means", "log_phi1", "log_phi2",
      "sigma", "stand_effect_raw"),
    check_arrays = TRUE
  )
  util$check_all_expectand_diagnostics(base_samples)
  
  # PPC histogram — overall
  par(mfrow = c(1, 1))
  util$plot_hist_quantiles(
    samples, "y_rep",
    bin_min        = 0,
    bin_max        = 600,
    bin_delta      = 20,
    baseline_values = stan_data_all$y
  )
  
  # PPC per stand
  par(mfrow = c(4, 3), mar = c(5, 5, 2, 1))
  for (f in 1:stan_data_all$F) {
    idxs <- stan_data_all$start_idxs[f]:stan_data_all$end_idxs[f]
    subsamples <- util$filter_expectands(
      samples,
      paste0("y_rep[", idxs, "]")
    )
    spp_name   <- stand_year_all$spp[idxs[1]]
    stand_name <- stand_year_all$stand[idxs[1]]
    util$plot_hist_quantiles(
      subsamples, "y_rep",
      bin_min         = 0,
      bin_max         = 500,
      bin_delta       = 10,
      baseline_values = stan_data_all$y[idxs],
      main            = paste0(stand_name, " (", spp_name, ")")
    )
  }
  
  # rho
  par(mfrow = c(1, 2), mar = c(5, 5, 1, 1))
  util$plot_expectand_pushforward(
    samples[["rho[1]"]], 50, flim = c(0, 1),
    display_name = bquote(rho * " (Low)")
  )
  util$plot_expectand_pushforward(
    samples[["rho[2]"]], 50, flim = c(0, 1),
    display_name = bquote(rho * " (High)")
  )
  
  # theta per species
  par(mfrow = c(2, 5), mar = c(5, 5, 2, 1))
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      samples[[paste0("theta1[", s, "]")]], 50, flim = c(0, 1),
      display_name = paste0("theta1 ", species_names[s])
    )
  }
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      samples[[paste0("theta2[", s, "]")]], 50, flim = c(0, 1),
      display_name = paste0("theta2 ", species_names[s])
    )
  }
  
  # log_means per species on natural scale
  par(mfrow = c(2, 5), mar = c(5, 5, 2, 1))
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      exp(samples[[paste0("log_means[", s, ",1]")]]), 50,
      display_name = paste0("lambda ", species_names[s])
    )
  }
  for (s in seq_along(species_names)) {
    util$plot_expectand_pushforward(
      exp(samples[[paste0("log_means[", s, ",2]")]]), 50,
      display_name = paste0("mu ", species_names[s])
    )
  }
  
  # sigma
  par(mfrow = c(1, 1))
  util$plot_expectand_pushforward(
    samples[["sigma_low"]], 50,
    display_name = expression(sigma)
  )
  
  # log_alpha per stand
  par(mfrow = c(4, 3), mar = c(5, 5, 2, 1))
  for (f in 1:stan_data_all$F) {
    spp_name   <- stand_year_all$spp[stan_data_all$start_idxs[f]]
    stand_name <- stand_year_all$stand[stan_data_all$start_idxs[f]]
    util$plot_expectand_pushforward(
      exp(samples[[paste0("log_alpha[", f, "]")]]), 50,
      display_name = paste0(stand_name, " (", spp_name, ")")
    )
  }
}



#Probability of being in mast or not
par(mfrow = c(4,3), mar = c(5,5,2,1))
for(i in 1:stan_data_all$F){
  idxs <- stan_data_all$start_idxs[s]:stan_data_all$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples,names,
                                       baseline_values = stan_data_all$y[idxs],
                                       main = paste0(stand_order[i], " (",
                                                     stand_elev$elevation[stand_elev$stand == stand_order[i]],
                                                     " m)")
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
  ABAM = list(low = list(mu = 0,   sd = 0.8),  high = list(mu = 5.2, sd = 0.6)),
  CANO = list(low = list(mu = 2.3, sd = 0.8),  high = list(mu = 5.8, sd = 0.8)),
  PSME = list(low = list(mu = 0,   sd = log(30)/2.57), high = list(mu = 5.2, sd = 0.5)),
  THPL = list(low = list(mu = 4,   sd = 0.5),  high = list(mu = 7.0, sd = 0.5)),
  TSHE = list(low = list(mu = 4,   sd = 0.5),  high = list(mu = 7.0, sd = 0.5))
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
  prior_list[[paste0("theta1_", sp)]] <- rbeta(n_prior, 5, 1)
  prior_list[[paste0("theta2_", sp)]] <- rbeta(n_prior, 1, 4)
}

# sigma
prior_list[["sigma"]] <- abs(rnorm(n_prior, 0, 0.5))   # half-normal (sigma > 0)

# # stand effects  (marginal: stand_effect_raw * sigma)
# prior_sigma <- abs(rnorm(n_prior, 0, 0.5))
# for (f in seq_len(stan_data_all$F)) {
#   prior_list[[paste0("stand_effect_", f)]] <- rnorm(n_prior, 0, 1) * prior_sigma
# }

# log_phi1, log_phi2 — stored on log scale in Stan, plot phi = exp(log_phi)
prior_list[["phi1"]] <- exp(rnorm(n_prior, log(2.0), 0.5))
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

post_list[["sigma"]] <- post_samples$sigma

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


y_rep <- rstan::extract(fit_all)$y_rep
bayesplot::ppc_dens_overlay(
  y = stan_data_all$y,
  yrep = y_rep[1:200, ]
)
# Also useful per species:
bayesplot::ppc_dens_overlay_grouped(
  y = stan_data_all$y,
  yrep = y_rep[1:200, ],
  group = stan_data_all$sp
)



# New diagnostics  --------------------------------------------------------


# See R-hat and ESS for every parameter
print(fit_all, pars = c("log_means", "theta1", "theta2", "sigma", "log_phi1", "log_phi2"))

# Which parameters have the worst R-hat?
summary_fit <- summary(fit_all)$summary
summary_fit[is.na(summary_fit[,"Rhat"]) | summary_fit[,"Rhat"] > 1.05, ]

# Traceplot — visually shows if chains are mixing
traceplot(fit_all, pars = c("log_means[2,1]", "log_means[2,2]", 
                            "theta1[1]", "theta2[1]", "sigma"))


# ESS per chain
sapply(1:4, function(i) {
  chain <- rstan::extract(fit_all, permuted = FALSE)[, i, ]
  apply(chain, 2, function(x) coda::effectiveSize(x))
})

table(stand_year_all$species_id)  # how many obs per species?


# How many observations per species?
table(stand_year_all$species_id)

# For CANO specifically — what do the raw counts look like?
hist(stand_year_all$y[stand_year_all$species_id == 2], 
     main = "CANO seed counts", xlab = "Count")

# Are there actually low-production years for CANO?
quantile(stand_year_all$y[stand_year_all$species_id == 2])


cano_data <- stand_year_all[stand_year_all$species_id == 2, ]

summary(cano_data$y)

cano_data %>%
  group_by(stand) %>%
  summarise(n_years = n(), mean_y = mean(y), max_y = max(y))


stand_year_all %>%
  filter(species_id == 2, stand %in% c("AE10", "AM16", "AR07")) %>%
  ggplot(aes(x = year, y = y, colour = stand)) +
  geom_line() +
  geom_point() +
  labs(title = "CANO high-producing stands over time",
       y = "Seed count") +
  theme_minimal()


stand_year_all %>%
  filter(species_id == 2, stand %in% c("AB08", "AV06", "AV14", "AX15", "PP17")) %>%
  ggplot(aes(x = year, y = y, colour = stand)) +
  geom_line() +
  geom_point() +
  labs(title = "CANO low-producing stands over time") +
  theme_minimal()



# PLOTS -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming your data frame is called stand_year_all and has columns:
# species_id, stand_id (or a stand identifier), year, y (seed count), area

# First, create species names if you only have IDs
species_names <- c("ABAM", "CANO", "PSME", "THPL", "TSHE")

# Add species name to your data
stand_year_all <- stand_year_all %>%
  mutate(species = species_names[species_id])

# ============================================================
# Plot 1: Distribution of seed counts by species (histograms)
# ============================================================

ggplot(stand_year_all, aes(x = y + 1)) +
  geom_histogram(aes(fill = species), bins = 50, alpha = 0.7) +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~ species, scales = "free_y", ncol = 2) +
  labs(
    title = "Seed Count Distribution by Species",
    x = "Seed count (log scale, +1)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# ============================================================
# Plot 2: Seed density (seeds/area) by species and stand
# ============================================================

# Calculate seed density
stand_year_all <- stand_year_all %>%
  mutate(seed_density = y / area)

# Create a stand identifier if you don't have one
# (assuming rows are grouped by stand based on start_idxs/end_idxs)
stand_year_all <- stand_year_all %>%
  group_by(species_id) %>%
  mutate(stand_id = cumsum(row_number() %in% start_idxs)) %>%
  ungroup()

ggplot(stand_year_all, aes(x = factor(stand_id), y = seed_density + 1)) +
  geom_boxplot(aes(fill = species), outlier.size = 0.5, alpha = 0.7) +
  scale_y_log10() +
  facet_wrap(~ species, scales = "free", ncol = 2) +
  labs(
    title = "Seed Density by Stand within Species",
    x = "Stand",
    y = "Seed density (seeds/area, log scale)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )

# ============================================================
# Plot 3: Time series by species and stand
# ============================================================

ggplot(stand_year_all, aes(x = year, y = y + 1, group = stand_id)) +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~ species, scales = "free_y", ncol = 2) +
  labs(
    title = "Seed Production Time Series by Stand",
    subtitle = "Each line = one stand",
    x = "Year",
    y = "Seed count (log scale)"
  ) +
  theme_minimal()

# ============================================================
# Plot 4: Heatmap of seed production (species × year)
# ============================================================

# Aggregate by species and year
species_year_summary <- stand_year_all %>%
  group_by(species, year) %>%
  summarise(
    mean_seeds = mean(y),
    median_seeds = median(y),
    .groups = "drop"
  )

ggplot(species_year_summary, aes(x = year, y = species, fill = log10(mean_seeds + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "log10(mean + 1)") +
  labs(
    title = "Mean Seed Production by Species and Year",
    x = "Year",
    y = "Species"
  ) +
  theme_minimal()

# ============================================================
# Plot 5: Density plots comparing species distributions
# ============================================================

ggplot(stand_year_all, aes(x = y + 1, fill = species, color = species)) +
  geom_density(alpha = 0.3) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Seed Count Density by Species",
    x = "Seed count (log scale)",
    y = "Density"
  ) +
  theme_minimal()

# ============================================================
# Plot 6: Summary statistics table
# ============================================================

summary_stats <- stand_year_all %>%
  group_by(species) %>%
  summarise(
    n_obs = n(),
    n_stands = n_distinct(stand_id),
    mean = round(mean(y), 1),
    median = median(y),
    sd = round(sd(y), 1),
    q10 = quantile(y, 0.10),
    q90 = quantile(y, 0.90),
    pct_zero = round(100 * mean(y == 0), 1),
    .groups = "drop"
  )

print(summary_stats)

# ============================================================
# Plot 7: Faceted by species AND stand (if not too many stands)
# ============================================================

# Only practical if you have a manageable number of stands
n_stands <- length(unique(stand_year_all$stand_id))

if (n_stands <= 30) {
  ggplot(stand_year_all, aes(x = year, y = y + 1)) +
    geom_col(aes(fill = species), width = 0.8) +
    scale_y_log10() +
    facet_grid(stand_id ~ species, scales = "free_y") +
    labs(
      title = "Seed Production by Stand and Species",
      x = "Year",
      y = "Seed count (log scale)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(size = 6, angle = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
    )
}


stand_year_all %>%
  group_by(spp, stand) %>%
  summarise(n_species = n_distinct(species_id)) %>%
  filter(n_species > 1)
