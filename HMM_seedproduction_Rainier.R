#HMM model
#Analysis on seed production for one specie over elevation (stand used as proxy)
#Eléonore Perret - February 2026


#TO DO: 
#Re do with same distribution and then check with a partial pooling for the difference between stands. Check the parameters differnece between model. Create git issue and with the summary adn the description of the model (likelihood) adn thoughts also the posterior check visuals (using mickeals diagnostics, or density together with histogramms). What I see and what to do next. WHcih parameters to change based on the fact that wach stand has its own distribution. 

library(dplyr)
library(ggplot2)
library(rstan)
library(posterior)
library(bayesplot)
library(tidyr)

getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")
list.files()

seed_data<-read.csv("SeedData_all.csv")

# PSME --------------------------------------------------------------------
#selecting one species 
#Changing my dataset to PSME where only stands with seeds -- removing AE10, AR07, PARA, SUNR, SPRY, AV14, AM16
psme_data<-seed_data%>%
  filter(spp=="PSME")%>%
  #Added this constraint
  filter(!stand %in% c("AE10", "AR07", "PARA", "SUNR", "SPRY", "AV14", "AM16"))

total_stand_psme<-psme_data%>%
  group_by(stand)%>%
  summarise(total_seeds=sum(total_viable_sds),na.rm=TRUE)

ggplot(total_stand_psme, aes(x = stand, y = total_seeds)) +
  geom_col() +
  theme_bw() +
  labs(x = "Stand",
       y = "Total viable seeds",
       title = "Total viable seeds per stand") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

summary(psme_data$total_viable_sds)
mean(psme_data$total_viable_sds) #6.48657
var(psme_data$total_viable_sds) #139.4104
#This is overdispersed Var>> mean

#I can see here that some stand have 0 seed production in 15 years (very unlikely that it will change) except for SPRY and AR07 (but to be honest I'm a bit suprised; in total it is 2 seeds found for the AR07 and 3 for the SPRY over the last years; what to do with that. I honestly cannot refute the misidentification because there are no adult tree anywhere close to this area so something to think about)

# HMM: One specie x One stand ------------------------------------------------------------
psme_data_simple<-psme_data%>%
  filter (stand=="TO04")%>%
  group_by(year) %>%
  summarise(y = sum(total_viable_sds), .groups = "drop")

stan_data <- list(
  T = length(psme_data_simple$y),
  y = as.integer(psme_data_simple$y)
)

mod <- cmdstan_model("hmm_simple.stan")

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4
)

fit$summary()

summary(fit$draws(c("phi","log_mu", "Gamma","rho")))
#where phi is the overdispersion for each state
#where log_mu is the mean seed production for each state
#where gamma is the propability of transition from one state to another
#where rho is the probability to be in a state which is 0,665 to be in low and 0,345 in high

# Preparing the data for plotting the posterior distribution
yrep <- fit$draws("y_rep")
yrep_mat <- as_draws_matrix(yrep)
yrep_quants <- apply(yrep_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))

years <- psme_data_simple$year
y_obs <- psme_data_simple$y

df_plot <- data.frame(
  year = years,
  y_obs = y_obs,
  lower = yrep_quants[1, ],
  median = yrep_quants[2, ],
  upper = yrep_quants[3, ]
)

#Plotting the posterior predictive check vs real data
ggplot(df_plot, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = median), color = "blue", size = 1) +
  geom_point(aes(y = y_obs), color = "red", size = 2) +
  labs(
    y = "Seed counts",
    x = "Year",
    title = "Posterior Predictive Check"
  ) +
  theme_minimal()


#Preparing to plot the posterior predictive check for each state
posterior <- fit$draws(variables = c("mu", "phi"), format = "df")
set.seed(123)
posterior_sub <- posterior[sample(nrow(posterior), 500), ]# Ssubset of posterior draws to reduce plotting size
# Simulating distributions
sim_data <- data.frame()

for (i in 1:nrow(posterior_sub)) {
  # State 1 (low)
  sim_low <- rnbinom(100, size = posterior_sub$`phi[1]`[i],
                     mu = posterior_sub$`mu[1]`[i])
  sim_data <- rbind(sim_data,
                    data.frame(seed = sim_low, state = "low"))

  # State 2 (high)
  sim_high <- rnbinom(100, size = posterior_sub$`phi[2]`[i],
                      mu = posterior_sub$`mu[2]`[i])
  sim_data <- rbind(sim_data,
                    data.frame(seed = sim_high, state = "high"))
}
#Plotting the posterior predictive check for each state
ggplot(sim_data, aes(x = seed, fill = state)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(name = "Seed Production") +
  scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
  theme_minimal() +
  ggtitle("Posterior Predictive Distributions of Seed Production per State")

#Adding the means seed production for each state
mu_means <- colMeans(posterior_sub[, c("mu[1]", "mu[2]")])
ggplot(sim_data, aes(x = seed, fill = state)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = mu_means[1], color = "blue", linetype = "dashed") +
  geom_vline(xintercept = mu_means[2], color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
  theme_minimal() +
  xlim(c(0,300))

#Density plot of seed production from real data
ggplot(psme_data_simple, aes(x = y)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  ggtitle("Density")

#Plotting everything togeter: Posterior distribution with real distribution
ggplot() +
  # Posterior predictive densities plot (by state)
  geom_density(data = sim_data, aes(x = seed, fill = state), alpha = 0.5) +
  geom_vline(xintercept = mu_means[1], color = "blue", linetype = "dashed") +
  geom_vline(xintercept = mu_means[2], color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
  # Actual observed data density
  geom_density(data = psme_data_simple, aes(x = y),
               fill = "black", alpha = 0.3) +
  theme_minimal() +
  xlim(c(0, 300)) +
  ggtitle("Posterior Predictive Densities (by state) with Observed Data")

# HMM: One specie x Almost all stands (see selection in PSME section) ------------------------------------------------------------
psme_data_simple2<-psme_data%>%
  group_by(year,stand) %>%
  summarise(y = sum(total_viable_sds), .groups = "drop")

missing_years <- psme_data_simple2 %>%
  group_by(stand) %>%
  summarise(
    start_year = min(year),
    end_year   = max(year),
    observed_years = list(sort(unique(year))),
    full_seq       = list(seq(start_year, end_year)),
    missing_years  = list(setdiff(full_seq[[1]], observed_years[[1]])),
    .groups = "drop"
  )

#As there is some discrepency with start and end year for each stand. 
stand_ranges <- psme_data_simple2 %>%
  group_by(stand) %>%
  summarise(
    start_year = min(year),
    end_year   = max(year),
    .groups = "drop"
  )

stand_ranges <- stand_ranges %>%
  mutate(
    year_seq = purrr::map2(start_year, end_year, seq)
  )

# Making sure the data is sorted properly
data_sorted <- psme_data_simple2 %>%
  arrange(stand, year)

#Defining some stuff
y <- data_sorted$y  
#Amount of years per stand
T_i <- data_sorted %>%
  count(stand) %>%
  pull(n)       

#Creating a ragged vector to make sure years, stands and observations match
start_idxs <- c()
end_idxs <- c()
id <- 1
for(s in 1:F){
  start_idxs <- c(start_idxs, id)
  id <- id + T_i[s]-1
  end_idxs <- c(end_idxs, id)
  id <- id + 1
}

#For stan code
data_list <- list(
  F   = length(T_i),#number of stands
  N   = length(y),#total number of observations
  start_idxs = start_idxs, 
  end_idxs = end_idxs, 
  y   = y
)

mod <- cmdstan_model("hmm_stand.stan")

fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4
)

fit$summary()
summary(fit$draws(c("phi","log_mu", "Gamma","rho")))

# Preparing the data for plotting the posterior distribution
yrep <- fit$draws("y_rep")
yrep_mat <- as_draws_matrix(yrep)

#Checking the parameters per stands
# Posterior mean per observation
yrep_mean <- colMeans(yrep_mat)
df_pred <- psme_data_simple2 %>%
  mutate(pred_mean = yrep_mean)

# Mean predicted per stand
df_pred %>%
  group_by(stand) %>%
  summarise(
    mean_observed = mean(y),
    mean_predicted = mean(pred_mean)
  )

#DENSITIES
# yrep_quants <- apply(yrep_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
# 
# df_plot <- psme_data_simple2 %>%
#   mutate(
#     lower  = yrep_quants[1, ],
#     median = yrep_quants[2, ],
#     upper  = yrep_quants[3, ]
#   )
# 
# # Plot per stand
# ggplot(df_plot, aes(x = year)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper),fill = "lightblue", alpha = 0.4) +
#   geom_line(aes(y = median),linewidth = 1, color = "blue") +
#   geom_point(aes(y = y),color = "red", size = 1.5) +
#   facet_wrap(~ stand, scales = "free_y") +
#   labs(y = "Seed counts",x = "Year",title = "Posterior Predictive Check by Stand") +
#   theme_minimal()
# 
# #Plot per state :
# posterior <- fit$draws(variables = c("mu", "phi"), format = "df")
# set.seed(123)
# posterior_sub <- posterior[sample(nrow(posterior), 500), ]
# sim_data <- data.frame()
# 
# for (i in 1:nrow(posterior_sub)) {
#   sim_low <- rnbinom(
#     100,
#     size = posterior_sub$`phi[1]`[i],
#     mu   = posterior_sub$`mu[1]`[i]
#   )
#   sim_high <- rnbinom(
#     100,
#     size = posterior_sub$`phi[2]`[i],
#     mu   = posterior_sub$`mu[2]`[i]
#   )
#   sim_data <- rbind(
#     sim_data,
#     data.frame(seed = sim_low,  state = "low"),
#     data.frame(seed = sim_high, state = "high")
#   )
# }
# 
# ggplot(sim_data, aes(x = seed, fill = state)) +
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(name = "Seed Production") +
#   scale_fill_manual(values = c("low" = "skyblue", "high" = "orange")) +
#   theme_minimal() +
#   ggtitle("Posterior Predictive Distributions of Seed Production per State")
# 
# #Densities with observed data
# mu_means <- colMeans(
#   posterior_sub[, c("mu[1]", "mu[2]")]
# )
# 
# ggplot() +
#   geom_density(data = sim_data,aes(x = seed, fill = state),alpha = 0.5) +
#   geom_density(data = psme_data_simple2,aes(x = y),fill = "black",alpha = 0.3) +
#   geom_vline(xintercept = mu_means[1],linetype = "dashed") +
#   geom_vline(xintercept = mu_means[2],linetype = "dashed") +
#   scale_fill_manual(values = c("low" = "skyblue","high" = "orange")) +
#   theme_minimal() +
#   xlim(c(0, 300)) +
#   ggtitle("Posterior Predictive Densities with Observed Data")
# 
# 
# ggplot() +
#   geom_histogram(data = df_hist,aes(x = observed, y = after_stat(density)),bins = 30,fill = "black",alpha = 0.4) +
#   geom_histogram(data = df_hist,aes(x = predicted, y = after_stat(density)), bins = 30,fill = "skyblue",alpha = 0.5) +
#   theme_minimal() +
#   xlim(c(0, 300)) +
#   labs(title = "Observed vs Posterior Predictive Histogram",x = "Seed Production",y = "Density" )
# 
# #Density per stand
# ggplot(psme_data_simple2, aes(x = y)) +
#   geom_density(fill = "grey70", alpha = 0.5) +
#   facet_wrap(~ stand, scales = "free") +
#   theme_minimal() +
#   ggtitle("Observed Seed Distribution per Stand")
# 
# 
# #states per stand
# z_draws <- as_draws_matrix(fit$draws("state"))
# prob_high <- colMeans(z_draws == 2)
# df_states <- psme_data_simple2 %>%
#   mutate(prob_high = prob_high)
# 
# ggplot(df_states, aes(x = year, y = prob_high)) +
#   geom_line() +
#   facet_wrap(~ stand) +
#   ylim(0,1) +
#   theme_minimal() +
#   labs( y = "Pr(High production state)",x = "Year",title = "State Probabilities by Stand")


#HISTROGRAMMS

set.seed(123)
draw_id <- sample(1:nrow(yrep_mat), 1)
yrep_one <- yrep_mat[draw_id, ]

df_hist <- data.frame(
  observed  = as.numeric(psme_data_simple2$y),
  predicted = as.numeric(yrep_one)
)
colnames(df_hist)

#Observed vs predicted
ggplot(df_hist) +
  geom_histogram(aes(x = observed, y = after_stat(density)),bins = 30,fill = "black",alpha = 0.4) +
  geom_histogram(aes(x = predicted, y = after_stat(density)),bins = 30,fill = "skyblue",alpha = 0.5) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 300)) +
  labs(title = "Observed vs Posterior Predictive Histogram",x = "Seed Production",y = "Density"
  )


df_hist2 <- psme_data_simple2 %>%
  mutate(predicted = yrep_one)

#observed vs predicted per stand
ggplot(df_hist2, aes(x = y)) +
  geom_histogram(bins = 20,fill = "black",alpha = 0.4) +
  geom_histogram(aes(x = predicted),bins = 20,fill = "skyblue",alpha = 0.5) +
  facet_wrap(~ stand, scales = "free") +
  theme_minimal() +
  labs(title = "Observed vs Posterior Predictive per Stand",x = "Seed Production")

#per states
posterior <- fit$draws(variables = c("mu", "phi"), format = "df")

set.seed(123)
posterior_sub <- posterior[sample(nrow(posterior), 300), ]

sim_data <- data.frame()

for (i in 1:nrow(posterior_sub)) {
  sim_low <- rnbinom(
    200,
    size = posterior_sub$`phi[1]`[i],
    mu   = posterior_sub$`mu[1]`[i]
  )
  sim_high <- rnbinom(
    200,
    size = posterior_sub$`phi[2]`[i],
    mu   = posterior_sub$`mu[2]`[i]
  )
  sim_data <- rbind(
    sim_data,
    data.frame(seed = sim_low,  state = "Low"),
    data.frame(seed = sim_high, state = "High")
  )
}

ggplot(sim_data, aes(x = seed, fill = state)) +
  geom_histogram(position = "identity",alpha = 0.5,bins = 40) +
  scale_fill_manual(values = c("Low" = "skyblue","High" = "orange")) +
  theme_minimal() +
  xlim(c(0, 300)) +
  labs(title = "Posterior Predictive Histograms by State",x = "Seed Production")

# HMM: One specie x Almost all stands (PARTIAL POOLING) ------------------------------------------------------------

#before my model assumed that both means from high and low states are identidical across stands, but now I will a stand-leved random effect. 

#For stan code
data_list <- list(
  F   = length(T_i),#number of stands
  N   = length(y),#total number of observations
  start_idxs = start_idxs, 
  end_idxs = end_idxs, 
  y   = y
)

mod <- stan_model("Stan_code/feb18/hmm_lowpoisson_highnegbin_standpooling.stan")

fit <- sampling(
  mod,
  data = data_list,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)

fit$summary()
summary(fit$draws(c("phi","log_mu", "Gamma","rho")))


####
##Plots
#Based on Mike's code 

util<- new.env()
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

# extraire les posterior values
samples <- util$extract_expectand_vals(fit)

# diagnostics parametre par parametre
base_samples <- util$filter_expectands(samples,
                                       c("rho", "Gamma_raw",
                                         "log_lambda", "log_mu",
                                         "stand_effect_raw", "phi",
                                         "sigma"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

# ppc for all 
par(mfrow = c(1,1))
names_yrep <- paste0("y_rep[", 1:data_list$N, "]")
util$plot_hist_quantiles(samples, "y_rep")
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 1000, bin_delta = 100)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = data_list$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = data_list$y)

par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:data_list$F){
  idxs <- data_list$start_idxs[s]:data_list$end_idxs[s]
  subsamples <-  util$filter_expectands(samples,
                                        paste0("y_rep[", idxs, "]"))
  util$plot_hist_quantiles(subsamples, "y_rep",
                           bin_min = 0, bin_max = 500, bin_delta = 10,
                           baseline_values = data_list$y[idxs])
}

#Probability of being in a mast vs non-mast 
par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:data_list$F){
  idxs <- data_list$start_idxs[s]:data_list$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples, names,
                                       baseline_values = data_list$y[idxs])
}

#Probability of being in a mast vs non-mast scaled
par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:data_list$F){
  idxs <- data_list$start_idxs[s]:data_list$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples, names,
                                       baseline_values = data_list$y[idxs],
                                       display_ylim = c(0,800))
}

unique(psme_data$stand)

#Probability of switching
par(mfrow = c(2,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["Gamma[1,1]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["Gamma[1,2]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["Gamma[2,1]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["Gamma[2,2]"]], 50, flim = c(0,1))

#Probability in starting in the low state
par(mfrow=c(1,2))
util$plot_expectand_pushforward(samples[["rho[1]"]], 50, flim=c(0,1))
util$plot_expectand_pushforward(samples[["rho[2]"]], 50, flim=c(0,1))


#Seed production for high state
for(s in 1:data_list$F){
  newname <- paste0("alpha[", s, "]")
  samples[[newname]] <- exp(samples[[paste0("log_alpha[", s, "]")]])
}
#Mast means vary by stand
par(mfrow = c(1,1), mar = c(5,5,1,1))
util$plot_disc_pushforward_quantiles(samples, paste0("alpha[", 1:data_list$F, "]"))
#Mast means vary by stand [1] or [2]
util$plot_expectand_pushforward(samples[["log_alpha[1]"]], 50)
util$plot_expectand_pushforward(samples[["log_alpha[2]"]], 50)
#Mast means vary by all stands
par(mfrow=c(4,3))
for(f in 1:F){
  util$plot_expectand_pushforward(samples[[paste0("log_alpha[", f, "]")]], 50)
}
#Seed production under mast year per stands.
par(mfrow=c(4,3))
for(f in 1:F){
  util$plot_expectand_pushforward(
    exp(samples[[paste0("log_alpha[", f, "]")]]), 50,display_name = bquote(alpha * "(high)"))
}

forplot <- psme_data %>%
  group_by(stand) %>%
  summarise(
    mean_seed = mean(total_viable_sds, na.rm = TRUE),
    max_seed  = max(total_viable_sds, na.rm = TRUE)
  )



#mean seed production of the low state
par(mfrow = c(1,1), mar = c(5,5,1,1))
util$plot_expectand_pushforward(exp(samples[["log_lambda"]]), 20,
                                display_name = bquote(lambda * " (low)"))



####


# Preparing the data for plotting the posterior distribution
yrep <- fit$draws("y_rep")
yrep_mat <- as_draws_matrix(yrep)

# Posterior mean per observation
yrep_mean <- colMeans(yrep_mat)

df_pred <- psme_data_simple2 %>%
  mutate(pred_mean = yrep_mean)

# Mean predicted per stand
df_pred %>%
  group_by(stand) %>%
  summarise(
    mean_observed  = mean(y),
    mean_predicted = mean(pred_mean)
  )


###First plot
##
#OVERALL
set.seed(123)
draw_id <- sample(1:nrow(yrep_mat), 1)
yrep_one <- yrep_mat[draw_id, ]

df_hist <- data.frame(
  observed  = as.numeric(psme_data_simple2$y),
  predicted = as.numeric(yrep_one)
)
#Plot for the posterior distribution vs the observed for all model
ggplot(df_hist) +
  geom_histogram(aes(x = observed, y = after_stat(density)),
                 bins = 30, fill = "black", alpha = 0.4) +
  geom_histogram(aes(x = predicted, y = after_stat(density)),
                 bins = 30, fill = "skyblue", alpha = 0.5) +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 300)) +
  labs(title = "Observed vs Posterior Predictive Histogram",
       x = "Seed Production",
       y = "Density")
#PER STAND
df_hist2 <- psme_data_simple2 %>%
  mutate(predicted = yrep_one)
#Plot for the posterior distribution vs the observed for each stand
ggplot(df_hist2, aes(x = y)) +
  geom_histogram(bins = 20, fill = "black", alpha = 0.4) +
  geom_histogram(aes(x = predicted), bins = 20, fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ stand, scales = "free") +
  theme_minimal() +
  labs(title = "Observed vs Posterior Predictive per Stand", x = "Seed Production")
##

##
#MEAN PER STATE and STAND
posterior <- fit$draws(variables = c("alpha","sigma_stand","stand_effect_raw","phi"), format = "df")
set.seed(123)
posterior_sub <- posterior[sample(nrow(posterior), 300), ]
#Simulation
sim_data <- data.frame()
F <- length(unique(psme_data_simple2$stand))
for (i in 1:nrow(posterior_sub)) {
  for (f in 1:F) {
    # Compute stand-specific mu for low and high
    mu_low  <- exp(posterior_sub[[paste0("alpha[1]")]][i] +
                     posterior_sub[[paste0("stand_effect_raw[",f,",1]")]][i] *
                     posterior_sub[[paste0("sigma_stand[1]")]][i])
    
    mu_high <- exp(posterior_sub[[paste0("alpha[2]")]][i] +
                     posterior_sub[[paste0("stand_effect_raw[",f,",2]")]][i] *
                     posterior_sub[[paste0("sigma_stand[2]")]][i])
    
    # Simulate for this stand
    sim_low  <- rnbinom(200, size = posterior_sub$`phi[1]`[i], mu = mu_low)
    sim_high <- rnbinom(200, size = posterior_sub$`phi[2]`[i], mu = mu_high)
    
    sim_data <- rbind(
      sim_data,
      data.frame(seed = sim_low, state = "Low", stand = f),
      data.frame(seed = sim_high, state = "High", stand = f)
    )
  }
}
# Compute mean for each stand & state
sim_means <- sim_data %>%
  group_by(stand, state) %>%
  summarise(mean_seed = mean(seed), .groups = "drop")
#Plot of the posterior mean per state and stand. 
ggplot(sim_data, aes(x = seed, fill = state)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
  # Add vertical lines for posterior predictive mean
  geom_vline(data = sim_means, aes(xintercept = mean_seed, color = state),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("Low" = "skyblue","High" = "orange")) +
  scale_color_manual(values = c("Low" = "blue", "High" = "red")) +
  facet_wrap(~ stand, scales = "free") +
  theme_minimal() +
  xlim(c(0, 300)) +
  labs(
    title = "Posterior Predictive Histograms by State and Stand",
    x = "Seed Production",
    subtitle = "Dashed lines = posterior predictive mean per state"
  )
##

##
#STATE switching
# Extract latent states from posterior draws
z_draws <- fit$draws("state")  # or "z" if that’s the variable in your Stan model
z_mat <- as_draws_matrix(z_draws)
# Take posterior mode per observation
z_mode <- apply(z_mat, 2, function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
})
# Add to dataset
df_pred <- df_pred %>%
  mutate(state_mode = z_mode)
# Count number of transitions per stand
state_transitions <- df_pred %>%
  group_by(stand) %>%
  summarise(
    n_obs = n(),
    n_transitions = sum(diff(state_mode) != 0),
    first_state = state_mode[1],
    last_state  = state_mode[n()]
  ) %>%
  arrange(n_transitions)
state_transitions
##


##
##STAND SPECIFIC MEAN VS OBSERVED DATA
post <- as_draws_matrix(
  fit$draws(c("alpha","sigma_stand","stand_effect_raw"))
)
stand_summary <- data.frame()
for (f in 1:F) {
  mu_draws <- exp(
    post[, paste0("alpha[1]")] +
      post[, paste0("stand_effect_raw[",f,",1]")] *
      post[, paste0("sigma_stand[1]")]
  )
  stand_summary <- rbind(
    stand_summary,
    data.frame(
      stand = f,
      mean  = mean(mu_draws),
      lower = quantile(mu_draws, 0.025),
      upper = quantile(mu_draws, 0.975)
    )
  )
}
# Compute mean observed seed production per stand
observed_summary <- psme_data_simple2 %>%
  group_by(stand) %>%
  summarise(
    mean_observed = mean(y),
    .groups = "drop"
  )
# Extract stand names from observed_summary
stand_names <- observed_summary$stand
# Add them to stand_summary
stand_summary$stand_name <- stand_names
observed_summary$stand_name <- stand_names
#Plot of the stand specific mean for posterior distribution vs the data
ggplot(stand_summary, aes(x = stand_name, y = mean)) +
  # Posterior mean ± CI
  geom_point(size = 3, color = "skyblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "skyblue") +
  # Observed mean per stand
  geom_point(data = observed_summary, aes(x = stand_name, y = mean_observed),
             size = 3, shape = 17, color = "black") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Stand",
    y = "Seed production",
    title = "Stand-specific posterior mean μ vs Observed data",
    subtitle = "Skyblue = posterior mean ± 95% CI, Black triangle = observed mean"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
##

##PLOTTING PARAMETERS
alpha_post  <- as_draws_matrix(fit$draws("alpha"))
sigma_post  <- as_draws_matrix(fit$draws("sigma_stand"))
phi_post    <- as_draws_matrix(fit$draws("phi"))
rho_post    <- as_draws_matrix(fit$draws("rho"))
gamma_post  <- as_draws_matrix(fit$draws("Gamma"))
summarise_param <- function(mat, param_name) {
  data.frame(
    parameter = colnames(mat),
    mean  = colMeans(mat),
    lower = apply(mat, 2, quantile, 0.025),
    upper = apply(mat, 2, quantile, 0.975),
    group = param_name
  )
}
summary_all <- bind_rows(
  summarise_param(alpha_post,  "Alpha (log mean)"),
  summarise_param(sigma_post,  "Sigma_stand"),
  summarise_param(phi_post,    "Phi (dispersion)"),
  summarise_param(rho_post,    "Rho (initial prob)"),
  summarise_param(gamma_post,  "Gamma (transition)")
)
summary_all$parameter <- gsub("\\[|\\]", "", summary_all$parameter)
ggplot(summary_all,
       aes(x = parameter, y = mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~ group, scales = "free_x") +   # << free x-axis
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # labels upright
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Posterior mean (95% CI)",
    title = "Posterior summaries of global HMM parameters"
  )
###






# HMM: One specie x Stands (Partial Pooled)  x Traps (as offset) x same start year and also with new transition matrix per stand ---------------------------------------

head(psme_data)

min(psme_data$size); max(psme_data$size);var(psme_data$size)
min(psme_data$total_viable_sds); max(psme_data$total_viable_sds);var(psme_data$total_viable_sds)

stand_year_df <- psme_data %>%
  group_by(stand, year) %>%
  filter(year!="2009") %>% #All starting in 2010 
  summarise(
    y = sum(total_viable_sds, na.rm = TRUE),
    area = sum(size, na.rm = TRUE), #Summarizing the total trap aerea
    .groups = "drop"
  ) %>%
  arrange(stand, year)


# Compute seed density (seeds per m²)
stand_year_df <- stand_year_df %>%
  mutate(seed_density = y / area)
# Seed density over time
ggplot(stand_year_df, aes(x = year, y = seed_density, color = stand, group = stand)) +
  geom_line(size = 1) +
  geom_point() +
  labs(
    title = "Seed Density (seeds/m²) per Stand over Time",
    x = "Year",
    y = "Seed Density"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

years_per_stand <- stand_year_df %>%
  group_by(stand) %>%
  summarise(T_i = n(), .groups = "drop")

F <- nrow(years_per_stand)
T_i <- years_per_stand$T_i

start_idxs <- c()
end_idxs <- c()
id <- 1

for (s in 1:F) {
  start_idxs <- c(start_idxs, id)
  id <- id + T_i[s] - 1
  end_idxs <- c(end_idxs, id)
  id <- id + 1
}

y <- stand_year_df$y
area <- stand_year_df$area
N <- length(y)

stan_data <- list(
  F = F,              
  N = N,              
  T_i = T_i,          
  start_idxs = start_idxs,
  end_idxs = end_idxs,
  y = y,              
  area = area        
)



mod1 <- stan_model("Stan_code/feb18/hmm_lowpoisson_highnegbin_standpooling.stan")

mod2 <- stan_model("Stan_code/feb18/hmm_lowpoisson_highnegbin_standpooling_TrapScaling.stan")

mod3 <- stan_model("Stan_code/feb18/hmm_lowpoisson_highnegbin_standpooling_TrapScaling_TransitionStand.stan")

mod4 <- stan_model("Stan_code/feb18/hmm_2highnegbin_standpooling_TrapScaling_TransitionStand.stan")

#Issues as some values of y rep in the initialization cannot be computed so I changed the log_mu_raw

fit1 <- sampling(
  mod1,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)

fit2 <- sampling(
  mod2,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)

fit3 <- sampling(
  mod3,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)

fit4 <- sampling(
  mod4,
  data = stan_data,
  chains = 4, cores = 4,
  iter = 2000, warmup = 1000
)

####
##Plots
#Based on Mike's code 

util<- new.env()
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)

# diagnostics generaux HMC (chain behavior)
diagnostics <- util$extract_hmc_diagnostics(fit1)
diagnostics <- util$extract_hmc_diagnostics(fit2)
diagnostics <- util$extract_hmc_diagnostics(fit3)
diagnostics <- util$extract_hmc_diagnostics(fit4)
util$check_all_hmc_diagnostics(diagnostics)

# extraire les posterior values
samples <- util$extract_expectand_vals(fit1)
samples <- util$extract_expectand_vals(fit2)
samples <- util$extract_expectand_vals(fit3)
samples <- util$extract_expectand_vals(fit4)

# diagnostics parametre par parametre
base_samples <- util$filter_expectands(samples,
                                       c("rho", "theta1","theta2",
                                         "log_lambda", "log_mu",
                                         "stand_effect_raw", "phi1", "phi2",                                         "sigma"), check_arrays = TRUE)
util$check_all_expectand_diagnostics(base_samples)

# ppc for all 
par(mfrow = c(1,1))

names_yrep <- paste0("y_rep[", 1:data_list$N, "]")
util$plot_hist_quantiles(samples, "y_rep")



# Plot PPC
par(mfrow = c(1,1))

util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 1000, bin_delta = 100)

util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 600, bin_delta = 20,
                         baseline_values = data_list$y)
util$plot_hist_quantiles(samples, "y_rep",
                         bin_min = 0, bin_max = 100, bin_delta = 2,
                         baseline_values = data_list$y)

par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:data_list$F){
  idxs <- data_list$start_idxs[s]:data_list$end_idxs[s]
  subsamples <-  util$filter_expectands(samples,
                                        paste0("y_rep[", idxs, "]"))
  util$plot_hist_quantiles(subsamples, "y_rep",
                           bin_min = 0, bin_max = 500, bin_delta = 10,
                           baseline_values = data_list$y[idxs])
}

#Probability of being in a mast vs non-mast 
par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:data_list$F){
  idxs <- data_list$start_idxs[s]:data_list$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples, names,
                                       baseline_values = data_list$y[idxs])
}

#Probability of being in a mast vs non-mast scaled
par(mfrow = c(4,3), mar = c(5,5,1,1))
for(s in 1:data_list$F){
  idxs <- data_list$start_idxs[s]:data_list$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  util$plot_disc_pushforward_quantiles(samples, names,
                                       baseline_values = data_list$y[idxs],
                                       display_ylim = c(0,800))
}

all(sort(unlist(mapply(seq,
                       data_list$start_idxs,
                       data_list$end_idxs))) == 1:data_list$N)


unique(psme_data$stand)

#Probability of switching
par(mfrow = c(2,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["Gamma[1,1]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["Gamma[1,2]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["Gamma[2,1]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["Gamma[2,2]"]], 50, flim = c(0,1))

par(mfrow = c(1,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["theta1[1]"]], 50, flim = c(0,1))
util$plot_expectand_pushforward(samples[["theta2[1]"]], 50, flim = c(0,1))
F <- length(unique(stand_year_df$stand))
par(mfrow = c(F,2), mar = c(3,3,1,1))  # adjust layout if you have many stands
for(f in 1:F){
  util$plot_expectand_pushforward(1 - samples[[paste0("theta1[",f,"]")]], 50, flim = c(0,1))
  title(paste0("Stand ", f, ": 1 → 2"))
  
  util$plot_expectand_pushforward(1 - samples[[paste0("theta2[",f,"]")]], 50, flim = c(0,1))
  title(paste0("Stand ", f, ": 2 → 1"))
}

#Probability in starting in the low state
par(mfrow=c(1,2))
util$plot_expectand_pushforward(samples[["rho[1]"]], 50, flim=c(0,1))
util$plot_expectand_pushforward(samples[["rho[2]"]], 50, flim=c(0,1))


#Seed production for high state
for(s in 1:data_list$F){
  newname <- paste0("alpha[", s, "]")
  samples[[newname]] <- exp(samples[[paste0("log_alpha[", s, "]")]])
}
#Mast means vary by stand
par(mfrow = c(1,1), mar = c(5,5,1,1))
util$plot_disc_pushforward_quantiles(samples, paste0("alpha[", 1:data_list$F, "]"))
# #Mast means vary by stand [1] or [2]
# util$plot_expectand_pushforward(samples[["log_alpha[1]"]], 50)
# util$plot_expectand_pushforward(samples[["log_alpha[2]"]], 50)
#Mast means vary by all stands
par(mfrow=c(4,3))
for(f in 1:F){
  util$plot_expectand_pushforward(samples[[paste0("log_alpha[", f, "]")]], 50)
}
#Seed production under mast year per stands.
par(mfrow=c(4,3))
for(f in 1:F){
  util$plot_expectand_pushforward(
    exp(samples[[paste0("log_alpha[", f, "]")]]), 50,display_name = bquote(alpha * "(high)"))
}

forplot <- psme_data %>%
  group_by(stand) %>%
  summarise(
    mean_seed = mean(total_viable_sds, na.rm = TRUE),
    max_seed  = max(total_viable_sds, na.rm = TRUE)
  )



#mean seed production of the low state
par(mfrow = c(1,1), mar = c(5,5,1,1))
util$plot_expectand_pushforward(exp(samples[["log_lambda"]]), 20,
                                display_name = bquote(lambda * " (low)"))



