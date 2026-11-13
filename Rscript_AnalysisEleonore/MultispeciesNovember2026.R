##Code for multispecies HMM model 
#HMM with two states each defined by a NB distribution
##Last updated : 12.11.2026
##Structure of the code 
##1. Libraries needed + setting working directory +Loading Mike's code
##2. Data import + modifications : Data needed to be added from the trees as not all the trees are present at each site and also data from the stand not surveyed by the PSP (SUNR, PARA and SPRY).Also adding the elevation. 
##3. Data preparation for Stan
##4. Model fitting
##5. Model diagnostic (model convergence, parameters, model summary and ppc)
##6. State assignments
##7. Posterior distribution of the parameters
##8. Prior vs posterior check
##9. Model results 
##10. Synchrony functions
##11. Intraspecific synchrony
##12. Interspecific synchrony


# 1. Library + Working directory + Mike's code-----------------------------------------------------------------
#Libraries
library(dplyr)
library(ggplot2)
library(cmdstanr)
library(rstan)
library(tidyr)
library(bayesplot)
library(tidyverse)
library(posterior)
library(purrr)


options(mc.cores = parallel::detectCores()) ##Important to run my code faster


#Setting working directory 
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

#Loading Mike's code
setwd("Rscript_AnalysisEleonore/")
#install.packages("colormap") #if not already installed
util <- new.env() #to import the R code from Mike that I have on the Rscript Analysis Eleonore folder
source("mcmc_visualization_tools.R", local=util)
source("mcmc_analysis_tools_rstan.R", local=util)
ls(util)


setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

# 2. Data import + modification -------------------------------------------------------------
seed_data<-read.csv("Data/SeedData_all.csv") #the merged seed data created by Janneke

# 1) Keeping species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","ABLA","CANO","PSME","TSHE","TSME","THPL"))

unique(seed_data$stand)

# 2) I want only possible combination of stands (ie. for CANO I don't want it at low elevation stands)
# a) Extra stand (to add manually as there is no map/records of the species present at SPRY, PARA and CANO as it was not installed by the PSP crew but Janneke)
# Add SPRY, SUNR, PARA for CANO, ABAM, ABLA
extra_stands <- seed_data %>%
  filter(spp %in% c("CANO", "ABAM", "ABLA"),
         stand %in% c("SPRY", "SUNR", "PARA")) %>%
  filter(year != 2009) %>%
  group_by(spp, stand, year) %>%
  summarise(y = sum(total_viable_sds, na.rm = TRUE),area = sum(size, na.rm = TRUE),.groups = "drop") %>%
  arrange(spp, stand, year)
# Checking
extra_stands %>%
  group_by(spp, stand) %>%
  summarise(n_years = n(), mean_density = mean(y/area), .groups = "drop")

# b) Getting the information from the PSP data
mapRainier<-read.csv("data/rawdata/Cleaned_mapping_2017 Rainier.csv", sep=";")

#Species per stand if present 
stands_per_species <- mapRainier %>%
  group_by(species) %>%  
  summarise(stands = list(unique(stand_id)),.groups = "drop") 

# Saving the object (I prefer this way but can be done also directly)
saveRDS(stands_per_species, file = "data/stands_per_species.rds")

#c) Adding only the possible combination of species x stands
stands_per_species <- readRDS("data/stands_per_species.rds")#from General.Data.R (orignal script but put there here for better readability)
stands_long <- stands_per_species %>%
  unnest(stands) %>%
  rename(spp = species,stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)
#Seeds_filtered contains now all the species x stand combination 

#d)Stand elevation
stand_elev <- data.frame(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13","AO03","AG05","AV06","AX15","AB08","PP17","AV14","AM16","AR07","PARA"  ,"SPRY","SUNR"),
  elevation = c(600,668,700,850,1450,850,900,950,1060,1090,1100,1150,1150,1200,1450,1600,1700,1800))
saveRDS(stand_elev, file= "data/stand_elevation_table.rds")


# 3. Data prep for Stan ---------------------------------------------------------------

#Removing year 2009, and merging all seed counts
stand_year_all <- seed_filtered %>%
  filter(year != 2009) %>% #because like this they start all the same year
  group_by(spp, stand, year) %>%
  summarise(y= sum(total_viable_sds, na.rm = TRUE),area = sum(size, na.rm = TRUE),.groups = "drop") %>%
  arrange(spp, stand, year)

#Adding the data from SUNR, SPRY and PARA
stand_year_all <- bind_rows(stand_year_all, extra_stands) %>%
  arrange(spp, stand, year)

#As I have a year missing for SUNR. 
missing_rows <- tibble(
  spp   = c("ABAM", "ABLA", "CANO"),
  stand = "SUNR",
  year  = 2023,
  y     = NA_integer_,
  area  = NA_real_
)

#Adding the missing year to the dataset
stand_year_all <- bind_rows(stand_year_all, missing_rows) %>%
  group_by(spp, stand) %>%
  mutate(area = ifelse(is.na(area), mean(area, na.rm = TRUE), area)) %>%
  ungroup() %>%
  arrange(spp, stand, year)

#Checking where I have missing "y" (sum of seed count)
obs_missing <- as.integer(is.na(stand_year_all$y))
#Replacing those NA values by 0 as I cannot have NA values in my dataset
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
         species_id = as.numeric(as.factor(spp))) 

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


#Final Stan data list 
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

# 4. Fitting Model -----------------------------------------------------------
#Current model used after discussion in git
fit_75 <- stan(
  file    = "Stan_code_Eleonore/Species_Stan_Model/MultispeciesStanCodeLast.stan",
  data    = stan_data_all,
  iter    = 4000, #change based on how much iterations you need
  warmup  = 2000, #make the warmup longer 
  chains  = 4,
  seed    = 123,
)


# 5. Model diagnostic --------------------------------------------------------

## a. Overall summary of the model fit. Here we can check for ESS (not bigger than the parameters themselves?) and Rhat (~1) values. 
summary(fit_75)
print(fit_75)

## b. Is the model good?
#Using Mike's function to check the fit of the model (if warning message, like divergence it would appear here)
diagnostics <- util$extract_hmc_diagnostics(fit_75) 
util$check_all_hmc_diagnostics(diagnostics)
#All is good. 

## c. We can check now the model parameters using Mike's function. 
samples <- util$extract_expectand_vals(fit_75)
names(samples)
base_samples1 <- util$filter_expectands(samples,c("rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2"))
base_samples2 <- util$filter_expectands(samples,c(paste0("alpha_theta1_species_nc[",1:S, "]"),paste0("alpha_theta2_species_nc[",   1:S, "]"), "sigma_theta1_species","sigma_theta2_species"))
base_samples3 <- util$filter_expectands(samples,c("mu_log_low",paste0('alpha_low_species[',1:S, ']'),paste0('alpha_stand[',rep(1:S, times= N_stands), ",", rep(1:N_stands, each = S), "]"), "sigma_low_species","sigma_stand"))
base_samples4<-util$filter_expectands(samples,c("mu_log_delta",paste0('log_delta_species_nc[',1:S, ']')),"sigma_log_delta_species","sigma_delta_stand")
base_samples5 <- util$filter_expectands(samples,c(paste0('log_phi_species[',1:S, ']'),"log_phi_state[1]","log_phi_state[2]"))
util$check_all_expectand_diagnostics(base_samples1);util$check_all_expectand_diagnostics(base_samples2); util$check_all_expectand_diagnostics(base_samples3) ; util$check_all_expectand_diagnostics(base_samples4) 
#It looks good no warning messages.  


#If the model diagnostics are well we can save the model summary: 

## d. Model summary: 
#Defining the parameters I'm interested in
model_parameters<- c("rho","grand_logit_theta1", "grand_logit_theta2","sigma_theta1_species", "sigma_theta2_species","mu_log_low","alpha_low_species", "sigma_low_species", "sigma_stand", "mu_log_delta","log_delta_species","sigma_log_delta_species","sigma_delta_stand","log_phi_species","log_phi_state")
summary_model<-summary(fit_75, pars= model_parameters, probs = c(0.025,0.25,0.5,0.75,0.975))$summary
model_table_summary<-as.data.frame(summary_model)%>% #model summary into a table
  tibble::rownames_to_column("parameter")%>%
  mutate(across(c(mean, se_mean, sd, '2.5%','25%','50%','75%', '97.5%'), ~round (.x, 3)),n_eff= round(n_eff),Rhat= round (Rhat,3))
model_table_summary #checking
write.table(model_table_summary, file = 'data/Tables/model_table_summary.txt', col.names = TRUE,row.names = FALSE, sep = "\t")#saving

#But is the model fitting our data? 

## e. Checking the PPC. 
#This allows us to know if the generated quantiles from the model match our actual data. 

# samples<-util$extract_expectand_vals(fit_75) #Not necessary if done in the previous step
#Overall fit: PPC
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep",bin_min = 0, bin_max = 600, bin_delta = 20,baseline_values = stan_data_all$y)

#Per series PPC
#We can also check per series : 
par(mfrow = c(4, 4))
for (f in 1:G) {
  idx    <- start_idxs[f]:end_idxs[f]
  years_f <- stand_year_all$year[idx] 
  pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
  util$plot_conn_pushforward_quantiles(samples,pred_names_f,years_f,xlab = paste("Series", f),ylab = "y")
  points(years_f, stan_data_all$y[idx], pch = 16, cex = 0.8)
}

#or for a single combination of F
f <- 26
idx     <- start_idxs[f]:end_idxs[f]
years_f <- stand_year_all$year[idx]
y_obs_f <- stan_data_all$y[idx]
pred_names_f <- sapply(idx, function(n) paste0("y_rep[", n, "]"))
util$plot_conn_pushforward_quantiles(samples, pred_names_f, years_f,xlab = paste("Series", f), ylab = "y")
points(years_f, pmin(y_obs_f, 400), pch = 16, cex = 0.8)# Adding point

#Everything looks good :) 

#Now we can go into the next step of checking the state identification. Does the model assign correctly the states (low and high?) and does it matches our data? 


# 6. State assignement -------------------------------------------------------

#Here we can see the model state assignements. I have made it into 3 different styles. a. you can see all the state assigment with the seed production for that same year and F combination (species x stand), you can ACrtl + Shift + C to uncomment this part if you want to look, b) seperated by species and stand ONLY state assignement c) seperated by species and stand ONLY seed production and d) specie seed production averaged. 
#I also have for the stands if we are interested in it :) in e) and finally, inf ) the seed production parameters and state transition parameters. 

#same as above no need to re-run if already done above. 
# samples<-util$extract_expectand_vals(fit_75) 

#a) all together 
# par(mfrow = c(4,6))
# for (f in 1:length(stan_data_all$start_idxs)) {
#   start_id <- stan_data_all$start_idxs[f]
#   end_id   <- stan_data_all$end_idxs[f]
#   species_id<-as.character(stand_year_all$spp[start_idxs[f]])
#   stand_id<-as.character(stand_year_all$stand[start_idxs[f]])
#   util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),display_ylim = c(1, 2))
#   title (main = paste(species_id, stand_id))
#   util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
#   points(x = 1:length(start_id:end_id), y = stan_data_all$y[start_id:end_id], pch = 20)
#   readline(prompt = paste0("Stand ", f, "/", length(stan_data_all$start_idxs), " -- Press [Enter] for next..."))
# }

#b) seperated state assignements
#the state assigment for each species at each stands. 
for (f in 1:length(stan_data_all$start_idxs)) {
  start_id <- stan_data_all$start_idxs[f]
  end_id   <- stan_data_all$end_idxs[f]
  idx<- start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", idx, "]"),display_ylim = c(1, 2), xticklabs = stand_year_all$year[idx],xlab= "Year", ylab= "State (1 = low, 2= high)",main = paste(stand_year_all$spp[start_id], stand_year_all$stand[start_id]))
}

#c) separated seed production per species and stands
#So here are the plots showing the model predictions for seed production (with greater uncertainty for the species and stands with lower values)
#I have a plot that can either show each species*stand combination with their own scale or if you want the same scale for all, we need to add this display_lim to the util$plot_disc_pushforward_quantiles(display_ylim = c(0,max(stan_data_all$y[stand_year_all$spp=="CANO"])))
#However this is sometimes not ideal because it crushes a lot the plots with lower seed production- but both can be displayed.

spp_f<-as.character(stand_year_all$spp[stan_data_all$start_idxs])
st_f<-as.character(stand_year_all$stand[stan_data_all$start_idxs])

#ABAM
for(f in which(spp_f=="ABAM")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]))
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}
#ABLA
par(mfrow=c(2,3))
for(f in which(spp_f=="ABLA")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]))
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}
#CANO
par(mfrow=c(2,2))
for(f in which(spp_f=="CANO")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]), display_ylim = c(0,max(stan_data_all$y[stand_year_all$spp=="CANO"])))
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}
#PSME
par(mfrow=c(4,3))
for(f in which(spp_f=="PSME")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]), #display_ylim = c(0,max(stan_data_all$y[stand_year_all$spp=="CANO"])))
  )
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}
#TSHE
par(mfrow=c(4,4))
for(f in which(spp_f=="TSHE")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]))
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}
#TSME
par(mfrow=c(2,3))
for(f in which(spp_f=="TSME")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]))
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}
#THPL
par(mfrow=c(4,3))
for(f in which(spp_f=="THPL")){
  start_id<-stan_data_all$start_idxs[f]
  end_id<-stan_data_all$end_idxs[f]
  idx<-start_id:end_id
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[",idx,"]"),xticklabs = stand_year_all$year[idx], xlab= "Year",ylab= "Seed count", main = paste(spp_f[f], st_f[f]))
  points(seq_along(idx),stan_data_all$y[idx], pch=18)
}


#d) Now a mean of all the stands per year per species. SO I have only one plot per species, all the stands seed production are together averaged.
par(mfrow=c(1,1))
plot_species_mean <- function(sp, samples, normalise_area = TRUE) {
  # Rows for this species, excluding missing observations
  sel <- which(stand_year_all$spp == sp & stan_data_all$obs_missing == 0)
  yrs <- sort(unique(stand_year_all$year[sel]))
  # so stands with different traps are comparable and because there are years with less or more traps than others
  scale_f <- if (normalise_area) stan_data_all$area / min(stan_data_all$area) else rep(1, length(stan_data_all$y))
  #mean y_rep across stands for each year
  pooled <- list()
  for (yr in yrs) {
    rows <- sel[stand_year_all$year[sel] == yr] #observations (stands) for this species and year
    arrs <- lapply(rows, function(r) samples[[paste0("y_rep[", r, "]")]] / scale_f[r]) #getting the posterior yrep for those stands
    pooled[[paste0("mean_yrep[", yr, "]")]] <- Reduce(`+`, arrs) / length(arrs) # each poaterior draw gets itw own stand average so it keeps the full uncertainity 
  }
  # Observed mean per year
  obs_mean <- tapply(stan_data_all$y[sel] / scale_f[sel], stand_year_all$year[sel], mean)
  #getting the plot
  util$plot_disc_pushforward_quantiles(
    pooled, names(pooled),
    xticklabs = yrs, xlab = "Year", ylab = "Mean seed count across stands",
    main = sp, display_ylim = c(0, max(obs_mean) * 1.15)) # * 1.15 to get a bit of margin above but not really necessary ifw we prefer not :) 
  points(seq_along(yrs), obs_mean, pch = 18)
}

plot_species_mean("ABAM", samples)
plot_species_mean("ABLA", samples)
plot_species_mean("CANO", samples)
plot_species_mean("PSME", samples)
plot_species_mean("TSHE", samples)
plot_species_mean("TSME", samples)
plot_species_mean("THPL", samples)


#State identification for the stand interspecific
#e) stands
par(mfrow=c(1,1))
plot_stand_mean <- function(st, samples, normalise_area = TRUE) {
  # Rows for this species, excluding missing observations
  sel <- which(stand_year_all$stand == st & stan_data_all$obs_missing == 0)
  yrs <- sort(unique(stand_year_all$year[sel]))
  # so stands with different traps are comparable and because there are years with less or more traps than others
  scale_f <- if (normalise_area) stan_data_all$area / min(stan_data_all$area) else rep(1, length(stan_data_all$y))
  #mean y_rep across stands for each year
  pooled <- list()
  for (yr in yrs) {
    rows <- sel[stand_year_all$year[sel] == yr] #observations (stands) for this species and year
    arrs <- lapply(rows, function(r) samples[[paste0("y_rep[", r, "]")]] / scale_f[r]) #getting the posterior yrep for those stands
    pooled[[paste0("mean_yrep[", yr, "]")]] <- Reduce(`+`, arrs) / length(arrs) # each poaterior draw gets itw own stand average so it keeps the full uncertainity 
  }
  # Observed mean per year
  obs_mean <- tapply(stan_data_all$y[sel] / scale_f[sel], stand_year_all$year[sel], mean)
  #getting the plot
  util$plot_disc_pushforward_quantiles(
    pooled, names(pooled),
    xticklabs = yrs, xlab = "Year", ylab = "Mean seed count across species",
    main = st, display_ylim = c(0, max(obs_mean) * 1.15)) # * 1.15 to get a bit of margin above but not really necessary ifw we prefer not :) 
  points(seq_along(yrs), obs_mean, pch = 18)
}

plot_stand_mean("TO04", samples)
plot_stand_mean("AV06", samples)
plot_stand_mean("AV06", samples)
#can do more... if this is something of interest :) 





# 7. Posterior distribution of the parameters -----------------------------

#distribution for each parameter. First the parameters for seed count , log_mu_low, log_delta, dispersion per state and per species, then for the state transition, rho, logit 1 et 2, high to low, low to low etc...

species <- c("ABAM","ABLA","CANO","PSME","THPL","TSHE","TSME")
B <- 25              # number of  bins

#Seed production parameters 
par(mfrow = c(2, 3))
util$plot_expectand_pushforward(samples[["mu_log_low"]], B, display_name = "Low seed production (log, global mean)")
util$plot_expectand_pushforward(samples[["mu_log_delta"]], B, display_name = "mu_log_delta (global gap, high vs low)")
util$plot_expectand_pushforward(samples[["sigma_low_species"]], B, display_name = "SD among species (low state)")
util$plot_expectand_pushforward(samples[["sigma_stand"]], B, display_name = "SD among stands (low state)")
util$plot_expectand_pushforward(samples[["sigma_log_delta_species"]], B, display_name = "SD among species (delta)")
util$plot_expectand_pushforward(samples[["sigma_delta_stand"]], B, display_name = "SD among stands (delta)")

# Species deviations in the low state
par(mfrow = c(3, 3))
for (s in 1:S)
  util$plot_expectand_pushforward(samples[[paste0("alpha_low_species[", s, "]")]], B,
                                  display_name = paste("Low-state deviation:", species[s]))

# Species deviations in delta (change between low and high)
par(mfrow = c(3, 3))
for (s in 1:S)
  util$plot_expectand_pushforward(samples[[paste0("log_delta_species[", s, "]")]], B,
                                  display_name = paste("Delta deviation:", species[s]))

# Dispersion per state 
par(mfrow = c(1, 2))
for (k in 1:2)
  util$plot_expectand_pushforward(samples[[paste0("log_phi_state[", k, "]")]], B,
                                  display_name = paste("log phi,", c("low", "high")[k], "state"))

# Dispersion per species (log scale)
par(mfrow = c(3, 3))
for (s in 1:S)
  util$plot_expectand_pushforward(samples[[paste0("log_phi_species[", s, "]")]], B,
                                  display_name = paste("log phi:", species[s]))

# Dispersion per species and state, on the natural scale --> phi = exp(species + state)
par(mfrow = c(4, 4))
for (s in 1:S) {
  phi_low  <- exp(samples[[paste0("log_phi_species[", s, "]")]] + samples[["log_phi_state[1]"]])
  phi_high <- exp(samples[[paste0("log_phi_species[", s, "]")]] + samples[["log_phi_state[2]"]])
  util$plot_expectand_pushforward(phi_low,  B, display_name = paste("phi low:",  species[s]))
  util$plot_expectand_pushforward(phi_high, B, display_name = paste("phi high:", species[s]))
}

#State transition parameters
# Initial state probabilities 
par(mfrow = c(1, 2))
util$plot_expectand_pushforward(samples[["rho[1]"]], B, display_name = "Start in low state")
util$plot_expectand_pushforward(samples[["rho[2]"]], B, display_name = "Start in high state")

# Global transition parameters and standard variation
par(mfrow = c(2, 2))
util$plot_expectand_pushforward(samples[["grand_logit_theta1"]], B, display_name = "Grand logit theta1 (stay low)")
util$plot_expectand_pushforward(samples[["grand_logit_theta2"]], B, display_name = "Grand logit theta2 (stay high)")
util$plot_expectand_pushforward(samples[["sigma_theta1_species"]], B, display_name = "SD among species (theta1)")
util$plot_expectand_pushforward(samples[["sigma_theta2_species"]], B, display_name = "SD among species (theta2)")

# Transition probabilities per species
for (s in 1:S) {
  theta1 <- plogis(samples[["grand_logit_theta1"]] + samples[[paste0("alpha_theta1_species[", s, "]")]])
  theta2 <- plogis(samples[["grand_logit_theta2"]] + samples[[paste0("alpha_theta2_species[", s, "]")]])
  
  par(mfrow = c(2, 2))
  util$plot_expectand_pushforward(theta1,     B, display_name = paste(species[s], ": low -> low"))
  util$plot_expectand_pushforward(1 - theta1, B, display_name = paste(species[s], ": low -> high"))
  util$plot_expectand_pushforward(1 - theta2, B, display_name = paste(species[s], ": high -> low"))
  util$plot_expectand_pushforward(theta2,     B, display_name = paste(species[s], ": high -> high"))
}

# 8. Prior vs posterior ------------------------------------------------

#Checking how much the data updated each priors. If a parameter barely moves, then the prior was weakly informed by the data. That is what we want. 

samples <- rstan::extract(fit_75) #extracting the posterior draws if not done already
n <- nrow(samples$mu_log_low) #number of posterior draws

#This function allows draws of random numbers with a minimum of 0. So only positive values.  first, it gives a probability of getting a value below 0 in that normal distribution (pnorm) and then it draws n random numbers unifromly beween that probability and 1 (runif). then qnorm converts those probabilities into values on the normal scale. This allows me to not have values lower then 0 for the sigma as using rnorn creates negative values. 
rtnorm_function<-function(n,mean, sd){ 
  qnorm(runif(n,pnorm(0,mean,sd),1), mean, sd)
}

#simulating draws from the priors (based on the .stan script)
set.seed(123)
priors <- data.frame(
  mu_log_low                = rnorm(n, 2, 1.5),
  mu_log_delta              = rnorm(n, 1.5, 1.5),
  grand_logit_theta1        = plogis(rnorm(n, 0.5, 1)),
  grand_logit_theta2        = plogis(rnorm(n, -1, 0.5)),
  sigma_low_species         = rtnorm_function(n, 0.5, 1),
  sigma_stand               = rtnorm_function(n,1,0.5),
  sigma_log_delta_species   = rtnorm_function(n, 1.5, 0.7),
  sigma_delta_stand         = rtnorm_function(n,1,0.5),
  sigma_theta1_species      = abs(rnorm(n, 0, 0.7)),
  sigma_theta2_species      = abs(rnorm(n, 0, 0.7))
)

#what are the actual draws. This is now taken from the model. and I back-transformed with plogis so everything (prior and posterior) are on the same scale. 
posteriors <- data.frame(
  mu_log_low                = samples$mu_log_low,
  mu_log_delta              = samples$mu_log_delta,
  grand_logit_theta1        = plogis(samples$grand_logit_theta1),
  grand_logit_theta2        = plogis(samples$grand_logit_theta2),
  sigma_low_species         = samples$sigma_low_species,
  sigma_stand               = samples$sigma_stand,
  sigma_log_delta_species   = samples$sigma_log_delta_species,
  sigma_delta_stand         = samples$sigma_delta_stand,
  sigma_theta1_species      = samples$sigma_theta1_species,
  sigma_theta2_species      = samples$sigma_theta2_species
)

#labels 
param_labels <- c(
  mu_log_low              = "low state mean",
  mu_log_delta            = "log delta (gap between states)",
  grand_logit_theta1      = "P stay low",
  grand_logit_theta2      = "P stay high",
  sigma_low_species       = "Species variation (low)",
  sigma_stand             = "Stands variation (low)",
  sigma_log_delta_species = "Species variation (delta)",
  sigma_delta_stand       = "Stands variation (delta)",
  sigma_theta1_species    = "Species transition variation (k=1)",
  sigma_theta2_species    = "Species transition variation (k=2)"
)

#Reshape to long for ggplot and add a label if prior or posterior
prior_long <- priors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>% #everything means all columnss. 
  mutate(type = "Prior")

post_long <- posteriors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(type = "Posterior")

#merging the two together. 
df <- bind_rows(prior_long, post_long) %>%
  mutate(
    label= param_labels [param],
    label = factor(label, levels= unname(param_labels)),
    type  = factor(type, levels = c("Prior", "Posterior"))
  )

#Plot
prior_posterior_plot<-ggplot(df, aes(x = value, fill = type, color = type)) +
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
prior_posterior_plot
#saving the plot
ggsave(filename= "prior_posterior_plot.png", plot= prior_posterior_plot,path = "Plots/Manuscript")


#Overdispersion parameter phi by species and state (because too difficult to implement above)
species<- c("ABAM" ,"ABLA" ,"CANO" ,"PSME", "THPL" ,"TSHE" ,"TSME")
state_name<-c("non-masting", "masting")

#the prior on the log scale is the sum of a species effect and a state effect so it is simulated as two draws and added together.Then exponential. 
set.seed(123)
prior_phi_sp<-data.frame(
  species= rep(rep(species, each=n),times=2), #each species repeated n times (and one per draw) and then the whole blovk repeated twice
  state= rep(state_name, each= n*length(species)),
  value= exp(rnorm(n*length(species) * 2, log(4), 0.6) + rnorm(n*length(species)*2, 0, 0.3)),
  type= "Prior"
)

#same here, giving a log phi for every draws  and species.
posterior_phi_sp<-data.frame(
  species= rep(rep(species, each=n),times=2),
  state= rep(state_name, each= n*length(species)),
  value= c(exp(as.vector(samples$log_phi_species + samples$log_phi_state[,1])),exp(as.vector(samples$log_phi_species + samples$log_phi_state[,2]))),
  type= "Posterior"
)

#combining them together
df_phi<-bind_rows(prior_phi_sp, posterior_phi_sp)%>%
  mutate(type= factor(type, levels = c("Prior", "Posterior")), state= factor(state, levels= state_name))

#plot
species_overdipersion_priorvspost<-ggplot(df_phi, aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(state~ species, scales = "free", ncol = 4) +
  coord_cartesian(xlim=c(0,20))+
  scale_fill_manual(values  = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
  scale_color_manual(values = c("Prior" = "#7B9FBF", "Posterior" = "#E07B54")) +
  labs(
    title = "Prior vs posterior distributions of the overdispersion parameter",
    x = NULL, y = "Density", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
species_overdipersion_priorvspost
ggsave(filename= "species_overdipersion_priorvspost.png", plot= species_overdipersion_priorvspost,path = "Plots/Manuscript")

# Model results -----------------------------------------------------------
samples <- rstan::extract(fit_75) #extracting the posterior draws if not done already
#Observed data summary and plot.
#I wanted to plot the mean seed count by species, by species and year 
sum_observed_data<- stand_year_all%>%
  group_by(year, spp)%>%
  summarise(mean_seed= mean(y), .groups = "drop")
sum_observed_data2<- stand_year_all%>%
  group_by(spp)%>%
  summarise(mean_seed= mean(y), .groups = "drop")
sum_seed_production<-ggplot(sum_observed_data, aes(x=year, y= log(mean_seed), colour =spp)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  scale_x_continuous( breaks = seq(min(stand_year_all$year), max(stand_year_all$year)))+
  scale_color_viridis_d(option= "viridis")+
  labs(x= "Year", y= "Mean Seed production across stands", color = "Species")+
  theme_classic()
sum_seed_production
ggsave(filename= "sum_seed_production.png", plot= sum_seed_production,path = "Plots/Manuscript")

#Stand table
#I wanted a stand*species table to see which species are at which stand (similar than what I had in #2. but here I have the overall data (with the not PSP stand))
species_stand<-stand_year_all%>%
  distinct(spp, stand)%>%
  group_by(spp)%>%
  summarise(stands= paste(stand, collapse = ","), .groups = "drop")
write.table(species_stand, file = 'data/Tables/species_stand.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#Species-level model results
#Here I'm taking the posterior draws from the model (fit 75) and for each species, summarize the seed production per state, the transition probabilities and the dispersion.  
#samples<-rstan::extract(fit_75) #not necessary if already done
species<- c("ABAM" ,"ABLA" ,"CANO" ,"PSME", "THPL" ,"TSHE" ,"TSME")
species_id<-stan_data_all$species_id
baseline_area<-min(stan_data_all$area)

#This is what I want in my table as measures- median, 90% credible interval, but I could also add the mean for example. 
summary_masting<-function(x,name){
  tibble(parameter= name, 
         median = median(x),
         lo90= unname(quantile(x,0.05)),
         hi90= unname(quantile(x,0.95)))
}

results<-list()
for(s in seq_along(species)) {
  f<-which (species_id == s) #the species and stand series 
  #seeds per state
  low <-rowMeans (exp(samples$log_alpha_low[,f, drop= FALSE]))/baseline_area
  high <-rowMeans (exp(samples$log_alpha_high[,f, drop= FALSE]))/baseline_area
  #transition probabilites
  theta1<-plogis(samples$grand_logit_theta1 + samples$alpha_theta1_species[,s])
  theta2<-plogis(samples$grand_logit_theta2 + samples$alpha_theta2_species[,s])
  #overdispersion
  phi_low<-exp(samples$log_phi_species[,s]+ samples$log_phi_state[,1])
  phi_high<-exp(samples$log_phi_species[,s]+ samples$log_phi_state[,2])
  #results
  results[[s]]<-bind_rows(
    summary_masting(low, "Non-masting"),
    summary_masting(high, "Masting"),
    summary_masting(high/low, "Fold change"),#maybe can remove? 
    summary_masting(theta1, "P stay low"), 
    summary_masting(theta2, "P stay high"),
    summary_masting((1-theta1)/((1-theta1)+(1-theta2)), "P time high"), #the longest strech of staying in mast, but could remove. 
    summary_masting(phi_low, "Dispersion low"),
    summary_masting(phi_high, "Dispersion high")
  ) %>%
    mutate(species= species[s])
}
#Putting this into a table
species_table<-bind_rows(results)%>% mutate (across(c(median, lo90, hi90), ~round (.x,2)))
species_table
write.table(species_table, file = 'data/Tables/species_table_model_results.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")



#here I checked the model results for the seed production vs the actual data, but not really sure it is necessary... 
# #checking the model seed production results vs real observed data
# sp_check<-"ABAM"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# abam_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,]))+ labs(title= sp_check, x= "Seed count (log scale)", y= "Density")  
# abam_observed_model
# ggsave(filename= "abam_observed_model.png", plot= abam_observed_model,path = "Plots/Manuscript")
# 
# sp_check<-"ABLA"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# abla_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
# abla_observed_model
# ggsave(filename= "abla_observed_model.png", plot= abla_observed_model,path = "Plots/Manuscript") 
# 
# 
# sp_check<-"CANO"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# cano_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
# cano_observed_model
# ggsave(filename= "cano_observed_model.png", plot= cano_observed_model,path = "Plots/Manuscript") 
# 
# 
# sp_check<-"PSME"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# psme_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
# psme_observed_model
# ggsave(filename= "psme_observed_model.png", plot= psme_observed_model,path = "Plots/Manuscript")
# 
# 
# sp_check<-"THPL"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# thpl_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
# thpl_observed_model
# ggsave(filename= "thpl_observed_model.png", plot= thpl_observed_model,path = "Plots/Manuscript")
# 
# 
# sp_check<-"TSHE"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# tshe_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,]))+ labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
# tshe_observed_model
# ggsave(filename= "tshe_observed_model.png", plot= tshe_observed_model,path = "Plots/Manuscript")
# 
# sp_check<-"TSME"
# selection<-which(stand_year_all$spp==sp_check)
# y_obs<-stan_data_all$y[selection]
# y_model<-samples$y_rep[,selection]
# set.seed(1)
# d<-(sample(nrow(y_model),200))
# tsme_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
# tsme_observed_model
# ggsave(filename= "tsme_observed_model.png", plot= tsme_observed_model,path = "Plots/Manuscript")




# Synchrony functions-------------------------------------------------


#pulling ths state for each draws (one row per posterior draw and one colum per stand/year/specie row)
state_draws <- rstan::extract(fit_75, pars = "state")$state

meta_all <- stand_year_all %>% #I named it meta as it is overall data but can be named as we want
  mutate(row_id = row_number()) %>%
  rename(species = spp)

get_unit_year_state <- function(draw_idx) { #function to get the species, stand, year and state per draw
  tibble(
    species = meta_all$species,
    stand   = meta_all$stand,
    year    = meta_all$year,
    state   = state_draws[draw_idx, ]
  )
}

# Synchrony function (based on Victor synchrony). Within one year, for one draw: the fraction of stands in whichever state is more common that year 
majority_state <- function(df) { 
  state_counts <- table(df$state) #how many rows are in each state
  max_count<-max(state_counts) #the max count in the majority state
  major_state_ID<-names(state_counts)[state_counts==max_count] # the most frequent state ID for that year and that draw
  #majority_proportion = max(state_counts) / sum(state_counts) # the proportion of rows in that majority state, removed for the moment
  tibble (major_state = paste(major_state_ID, collapse= "/"), proportion= max_count / sum(state_counts))
}

#takes one row per drawx group and returns one row per group. 
summarise_draws2 <- function(df, value_col, group_col) { 
  df %>%
    group_by(across(all_of(group_col))) %>% #will split the group by the variable I like with ABAM year for example.
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE), #the median of proportions across all draws that year
      lo90   = quantile(.data[[value_col]], 0.05, na.rm = TRUE), #the 5th percentile of proportion across draws
      hi90   = quantile(.data[[value_col]], 0.95, na.rm = TRUE), #the 95th percentile
      majority_value = names(table(.data[["major_state"]]))[which.max(table(.data[["major_state"]]))], #summarizes the major state across draws
      conf = max(table(.data[["major_state"]]))/n(), #the share of draws in which the winning state was the majority,
      p_tie= mean(grepl("/", .data[["major_state"]])),
      .groups = "drop"
    )
}


# Intraspecfic ------------------------------------------------------------
##Here I first explain how I calculate the intraspecific synchrony (for one specie across stands in which it is present). I show first for ABAM and then do the same for the other 6 species. I could do a loop (I did that for the stands ) but the draws take time to be sampled ("abamm draws") and it results in some computing time, which, if there is a mistake made me loose a lot of time :P so for the moment I have them like this also, so I can check them better individually. But of course I can also make a loop later on. 
##In general, I pull out the state assignements for each draws for a specie (so year x stand) and make a summary (see synchrony functions)that I then plot.Then I also check which stands I have for a specific year so I an overview. and then I have which stands are within the majority state. 

####ABAM
target_species <- "ABAM" #selecting the species
target_meta <- meta_all %>% #selecting from my loop
  filter(species == target_species)
n_draws_total <-nrow(state_draws) #total number of draws
draw_idx_use<- seq(1,n_draws_total, by=4)#so it is faster

#mapdfr from the purr package: It takes list and vector and applies a function to each element combining afterwards everything into a single dataframe. 
abam_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>% #latent state of every stand/species per year
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    ungroup()%>%
    mutate(draw = d)
})
abam_summary <- summarise_draws(abam_draws, "proportion", "year")
write.table(abam_summary, file = 'data/Tables/abam_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
abam_stand_state<-abam_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop= FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
#This is better because it takes into account if a year there was no majority state. 
year_prop <- abam_stand_state %>%
  group_by(year) %>%
  summarise(n_stands = n(),
            n_high   = sum(state == 2),
            prop_high = n_high / n_stands,
            split_year = prop_high == 0.5,            
            .groups = "drop")
#Stands within the majority state
abam_majority_standsID<-abam_stand_state%>%
  left_join(abam_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(abam_majority_standsID, file = 'data/Tables/abam_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")


####ABLA
target_species <- "ABLA"
target_meta <- meta_all %>%
  filter(species == target_species)
target_state_draws <- state_draws[, target_meta$row_id]
n_draws_total <-nrow(target_state_draws) 
abla_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})
abla_summary <- summarise_draws(abla_draws, "proportion", "year")
write.table(abla_summary, file = 'data/Tables/abla_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
abla_majority_standsID<-target_stand_state%>% #stands within the majority state
  left_join(abla_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(abla_majority_standsID, file = 'data/Tables/abla_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")



####CANO
target_species <- "CANO"
target_meta <- meta_all %>%
  filter(species == target_species)
target_state_draws <- state_draws[, target_meta$row_id]
n_draws_total <-nrow(target_state_draws) 
cano_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})
cano_summary <- summarise_draws(cano_draws, "proportion", "year")
write.table(cano_summary, file = 'data/Tables/cano_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
cano_majority_standsID<-target_stand_state%>%#stands within the majority state
  left_join(cano_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(cano_majority_standsID, file = 'data/Tables/cano_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

####PSME
target_species <- "PSME"
target_meta <- meta_all %>%
  filter(species == target_species)
target_state_draws <- state_draws[, target_meta$row_id]
n_draws_total <-nrow(target_state_draws) 
psme_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})
psme_summary <- summarise_draws(psme_draws, "proportion", "year")
write.table(psme_summary, file = 'data/Tables/psme_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
psme_majority_standsID<-target_stand_state%>%#stands within the majority state
  left_join(psme_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(psme_majority_standsID, file = 'data/Tables/psme_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

####TSHE
target_species <- "TSHE"
target_meta <- meta_all %>%
  filter(species == target_species)
target_state_draws <- state_draws[, target_meta$row_id]
n_draws_total <-nrow(target_state_draws) 
tshe_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})
tshe_summary <- summarise_draws(tshe_draws, "proportion", "year")
write.table(tshe_summary, file = 'data/Tables/tshe_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
tshe_majority_standsID<-target_stand_state%>%
  left_join(tshe_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(tshe_majority_standsID, file = 'data/Tables/tshe_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

####TSME
target_species <- "TSME"
target_meta <- meta_all %>%
  filter(species == target_species)
target_state_draws <- state_draws[, target_meta$row_id]
n_draws_total <-nrow(target_state_draws) 
tsme_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})
tsme_summary <- summarise_draws(tsme_draws, "proportion", "year")
write.table(tsme_summary, file = 'data/Tables/tsme_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
tsme_majority_standsID<-target_stand_state%>%
  left_join(tsme_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(tsme_majority_standsID, file = 'data/Tables/tsme_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

####THPL
target_species <- "THPL"
target_meta <- meta_all %>%
  filter(species == target_species)
target_state_draws <- state_draws[, target_meta$row_id]
n_draws_total <-nrow(target_state_draws) 
thpl_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    mutate(draw = d)
})
thpl_summary <- summarise_draws(thpl_draws, "proportion", "year")
write.table(thpl_summary, file = 'data/Tables/thpl_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
#the average of the latent state
target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)
thpl_majority_standsID<-target_stand_state%>%
  left_join(thpl_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")
write.table(thpl_majority_standsID, file = 'data/Tables/thpl_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")


#merging tables
abam_table<-abam_summary%>%
  left_join(abam_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean, lo90, hi90, majority_value, n_in_majority, stands)
abla_table<-abla_summary%>%
  left_join(abla_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean,lo90, hi90, majority_value, n_in_majority, stands)
cano_table<-cano_summary%>%
  left_join(cano_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean,lo90, hi90, majority_value, n_in_majority, stands)
psme_table<-psme_summary%>%
  left_join(psme_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean,lo90, hi90, majority_value, n_in_majority, stands)
tshe_table<-tshe_summary%>%
  left_join(tshe_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean,lo90, hi90, majority_value, n_in_majority, stands)
tsme_table<-tsme_summary%>%
  left_join(tsme_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean,lo90, hi90, majority_value, n_in_majority, stands)
thpl_table<-thpl_summary%>%
  left_join(thpl_majority_standsID, by= c("year", "majority_value"))%>%
  select(year, median, mean,lo90, hi90, majority_value, n_in_majority, stands)

all_species_intraspecific_synchrony_table<-bind_rows(ABAM=abam_table,ABLA= abla_table, CANO= cano_table, PSME= psme_table, THPL= thpl_table, TSHE= tshe_table,TSME= tsme_table, .id= "species")%>%
  select (species, year, median, mean,lo90, hi90, majority_value, n_in_majority, stands )

all_species_intraspecific_synchrony_table<- all_species_intraspecific_synchrony_table%>%
  mutate(across(where (is.numeric), ~round(.x, 4)))

write.table(all_species_intraspecific_synchrony_table, file = 'data/Tables/all_species_intraspecific_synchrony_table.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")


#Merging all the plots for the intraspecific synchrony together
species_order <- c("THPL", "TSHE", "PSME", "ABAM", "CANO", "TSME", "ABLA")
#binding all the summaries together
all_majority <- bind_rows(
  THPL = thpl_summary, TSHE = tshe_summary, PSME = psme_summary,
  ABAM = abam_summary, CANO = cano_summary, TSME = tsme_summary,
  ABLA = abla_summary,
  .id = "species"                               
) %>%
  mutate(species = factor(species, levels = species_order))
#final plot
all_majority_plot <- ggplot(all_majority, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), alpha = 0.4, fill = "lightblue") +
  coord_cartesian(ylim = c(0.5, 1)) +
  facet_wrap(~ species, ncol = 2,
             labeller = labeller(species = latin_names)) +   
  labs(x = "Year", y = "Proportion of stands in the majority state",
       title = "Intraspecific synchrony") +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "italic", hjust = 0))
all_majority_plot

ggsave("all_majority_plot.png", all_majority_plot, path = "Plots/Manuscript",
       width = 8, height = 10)






#Overall mean (checking when there is a difference if per stand per year and not)
overall_intra_mean_species<-function(draws, sp){
  w<-meta_all%>%
    filter(species==sp)%>%
    count(year, name= "n_units")
  
  draws %>% 
    left_join(w, by= "year") %>%
    group_by(draw) %>%
    summarise (mean_prop = mean(proportion), mean_prop_weighted= weighted.mean(proportion, n_units),.groups= "drop") %>%
    summarise(species= sp, mean= mean(mean_prop), median= median(mean_prop), lo90= quantile(mean_prop, 0.05), hi90= quantile(mean_prop, 0.95),   mean_wght= mean(mean_prop_weighted), median_wght= median(mean_prop_weighted), lo90_wght= quantile(mean_prop_weighted, 0.05), hi90_wght= quantile(mean_prop_weighted, 0.95))
}

all_species_overall_intra_std<-bind_rows(
  overall_intra_mean_species(abam_draws, "ABAM"),
  overall_intra_mean_species(abla_draws, "ABLA"),
  overall_intra_mean_species(cano_draws, "CANO"),
  overall_intra_mean_species(psme_draws, "PSME"),
  overall_intra_mean_species(thpl_draws, "THPL"),
  overall_intra_mean_species(tshe_draws, "TSHE"),
  overall_intra_mean_species(tsme_draws, "TSME"),
)

all_species_overall_intra_std<-all_species_overall_intra_std%>%
  mutate(across(where (is.numeric), ~round(.x, 4)))

write.table(all_species_overall_intra_std, file = 'data/Tables/all_species_overall_intra_std.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")



# Interspecific synchrony -------------------------------------------------
draw_idx_use<- seq(1,n_draws_total, by=4)#so it is faster
n_draws_total <-nrow(state_draws)



#A) Overall synchrony (all species and stands together)
overall_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    ungroup()%>%
    mutate(draw = d)
})
overall_summary <- summarise_draws2(overall_draws, "proportion", "year")
overall_summary<-overall_summary%>%
  mutate(across(where (is.numeric), ~round(.x, 4)))
write.table(overall_summary, file = 'data/Tables/overall_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")
overall_synchrony<-ggplot(overall_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = NULL, y = "Proportion in the majority state",
    title = "Synchrony of the overall dataset",
    subtitle = "Species and Stand combined"
  ) +
  theme_minimal(base_size = 13)
overall_synchrony
ggsave(filename= "overall_synchrony.png", plot= overall_synchrony,path = "Plots/Manuscript")

#B) Intraspecific 
#1. defining for each species and year the majority state across all stands
#this is the intraspecific synchrony, so for each species and year the majority state across its stands
species_year_draws<- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    group_by(species, year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    ungroup()%>%
    mutate(draw = d)
})
species_year_summary<- summarise_draws2(species_year_draws, "proportion", c("species", "year"))
species_year_summary<-species_year_summary%>%
  mutate(across(where (is.numeric), ~round(.x, 4)))
write.table(species_year_summary, file = 'data/Tables/species_year_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

species_order <- c("THPL", "TSHE", "PSME", "ABAM", "CANO", "TSME", "ABLA")
latin_names <- c(THPL = "Thuja plicata",TSHE = "Tsuga heterophylla",PSME = "Pseudotsuga menziesii",ABAM = "Abies amabilis",CANO = "Callitropsis nootkatensis",TSME = "Tsuga mertensiana",ABLA = "Abies lasiocarpa")
species_year_summary_plot <- species_year_summary %>%
  mutate(species = factor(species, levels = species_order)) %>%   # panel order
  ggplot(aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), alpha = 0.4, fill = "lightblue") +
  facet_wrap(~ species, labeller = labeller(species = latin_names)) +  # Latin panel titles
  labs(
    x = "Year", y = "Proportion in the majority state",
    title = "Intraspecific synchrony"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "italic", hjust = 0))
species_year_summary_plot
ggsave(filename= "species_year_summary_plot.png", plot= species_year_summary_plot,path = "Plots/Manuscript")

#per species year and stands. 
#Information about the stand and year (how many and which ones) 
year_stands_species<-meta_all%>%
  group_by(year, species)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))
write.table(year_stands_species, file = 'data/Tables/year_stands_species.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#stand level states and majority stand ID. 
#all species intraspecific synchrony table 
#boxplot
#overall mean table. 


#Interspecific synchrony
#first, I take each species majority state per year
#Here, per year, the majortiy state across species (species all together)
species_overall_draws<- species_year_draws%>%
  filter(!grepl("/",major_state))%>%
  mutate(state=as.numeric(major_state))%>%
  group_by(draw,year)%>%
  filter(n()>=2)%>%
  group_modify(~majority_state(.x))%>%
  ungroup()

interspecifc_summary<-summarise_draws2(species_overall_draws, "proportion", "year")

interspecifc_plot<-ggplot(interspecifc_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
  x = "year", y = "Proportion of species in the majority state",
  title = "Interspecific synchrony",
) +
  theme_minimal(base_size = 13)
interspecifc_plot
ggsave(filename= "interspecifc_plot.png", plot= interspecifc_plot,path = "Plots/Manuscript")

#Now which species were in the majority each year
#I can now compare the major state per year and check which species matched that. 
species_in_majority_interspecific<-species_year_summary%>%
  select(species, year, sp_state= majority_value)%>%
  left_join(interspecifc_summary%>%select(year, majority_value, conf, median), by= "year")%>%
  filter(sp_state== majority_value)%>%
  group_by(year,majority_value,conf,median)%>%
  summarise(n_in_majority= n(), species= paste(sort(species), collapse= ","), .groups = "drop")

interspecifc_table<-interspecifc_summary%>% 
  left_join(species_in_majority_interspecific%>%select(year,majority_value,n_in_majority, species), by=c("year","majority_value"))%>%
  select(year, median,mean, lo90, hi90, majority_value, p_tie, n_in_majority, species)

interspecifc_table<-interspecifc_table%>%
  mutate(across(where (is.numeric), ~round(.x, 4)))

write.table(interspecifc_table, file = 'data/Tables/interspecifc_table.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")


#Interspecific synchrony within a stand


#Stand loop

stand_loop<-function(st){
  
  target_meta <- meta_all %>%
    filter(stand == st)
  
  #majority state across species per year
  stand_draws <- map_dfr(draw_idx_use, function(d) {
    get_unit_year_state(d) %>%
      filter(stand == st) %>%
      group_by(year) %>%
      filter(n() >= 2) %>%
      group_modify(~majority_state(.x)) %>%
      ungroup() %>%
      mutate(draw = d)
    })
  
  if (nrow(stand_draws)==0){
    message("skipped",st, ":because only one specie so no comparison possible")
    return(NULL)
  }
  
  stand_summary <- summarise_draws2(stand_draws, "proportion", "year")
  write.table(stand_summary, file = paste0("data/Tables/",st,"_summary.txt"), col.names = TRUE,
              row.names = FALSE, sep = "\t")
  
  stand_majority_plot<-ggplot(stand_summary, aes(x = year, y = median)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 1.8) +
    geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+
    labs(
      x = "year", y = "Proportion of species in the majority state",
      title = "Interspecific synchrony across stands"
    ) +
    theme_minimal(base_size = 13)

  ggsave(filename= paste0(st,"_majority.png"), plot= stand_majority_plot,path = "Plots/Manuscript")
  
  
  target_species_state<-target_meta%>%
    mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
    select(year, species, mean_state, state)
  
  stand_majority_speciesID<-target_species_state%>%
    left_join(stand_summary%>%select(year,majority_value, conf), by= "year")%>%
    filter(as.character(state)==majority_value)%>%
    group_by(year,majority_value, conf)%>%
    summarise(n_in_majority=n(),
              species= paste(sort(species),collapse = ","),.groups = "drop")
  
  write.table(stand_majority_speciesID, file = paste0( "data/Tables/",st,"_majority_speciesID.txt"), col.names = TRUE,
              row.names = FALSE, sep = "\t")
  
  list(summary=stand_summary, plot= stand_majority_plot, species= stand_majority_speciesID)
  

}


all_stands<-unique(stand_year_all$stand)

stand_results<-map(all_stands, stand_loop)%>% #this takes some time. 
  set_names(all_stands)


stand_table_all<-imap_dfr(stand_results, function(results, st){
  results$summary%>%
    left_join(results$species %>% select (year, n_in_majority, species), by= "year")%>%
    mutate(stand= st) %>%
    select(stand, year, median, mean, lo90, hi90, majority_value, n_in_majority, species, conf)
})

stand_table_all<-stand_table_all%>%
  mutate(across(where (is.numeric), ~round(.x, 4)))

write.table(stand_table_all, "data/Tables/stand_table_all.txt", col.names = TRUE,
             row.names = FALSE, sep = "\t")

stand_order<-c("TO11" , "TA01", "TO04", "TB13","AV02" ,"AO03" ,"AG05" ,"AV06", "AB08" ,"AX15", "AV14" ,"PP17" ,"AM16", "AE10" ,"AR07", "PARA", "SPRY", "SUNR" )
stand_table_all<-stand_table_all%>%
  mutate(stand= factor(stand, levels= stand_order))

#plots
stand_interspecific_synchrony<-ggplot(stand_table_all, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  facet_wrap(~stand)+
  labs(
    x = "year", y = "Proportion of species in the majority state",
    title = "Interspecific synchrony",
    subtitle = "Within Stands"
  ) +
  theme_minimal(base_size = 13)
stand_interspecific_synchrony
ggsave(filename= "stand_interspecific_synchrony.png", plot= stand_interspecific_synchrony,path = "Plots/Manuscript")

stand_interspecific_synchrony_mean<-ggplot(stand_table_all, aes(x = stand, y = mean)) +
  geom_boxplot(fill= "lightblue",alpha = 0.7) +
  geom_jitter()+
  stat_summary(fun= mean, geom= "point", shape= 21, fill= "blue", size=3.5) +
  labs(
    x = "Stands", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony"
  ) +
  theme_minimal(base_size = 13)
stand_interspecific_synchrony_mean
ggsave(filename= "stand_interspecific_synchrony_mean.png", plot= stand_interspecific_synchrony_mean,path = "Plots/Manuscript")



stand_lines <- ggplot(stand_table_all,
                      aes(x = year, y = median, colour = stand, group = stand)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = seq(min(stand_table_all$year),
                                  max(stand_table_all$year))) +
  scale_colour_viridis_d(option = "turbo", name = "Stand") +   # 18 distinguishable colours
  coord_cartesian(ylim = c(0.5, 1)) +                          # a majority proportion is always >= 0.5
  labs(
    x = "Year", y = "Proportion of species in the majority state",
    title = "Interspecific synchrony within stands"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))     # years do not overlap
stand_lines
ggsave(filename = "stand_interspecific_synchrony_lines.png", plot = stand_lines,
       path = "Plots/Manuscript", width = 10, height = 6)



stand_overview<-stand_year_all%>%
  group_by(stand)%>%
  summarise(
    species_present= paste(sort(unique(spp)), collapse = ","),
    n_species= n_distinct(spp),
    years_surveyed= paste(sort(unique(year)), collapse= ","),
    total_years= n_distinct(year),
    .groups= "drop"
  )



write.table(stand_overview, "data/Tables/stand_overview.txt", col.names = TRUE,
            row.names = FALSE, sep = "\t")




