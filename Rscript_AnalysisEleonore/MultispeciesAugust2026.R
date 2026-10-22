##Code for multispecies HMM model 
#HMM with two states each defined by a NB distribution
#Start date: 12.03.2026

##Next things to do : 

#1 : Prepare plots for manuscript
  #Figure of the time series of seed production overall
  #Figure of the trace plots of chain convergence
  #Fifure of state probability for every stands and species
  #Figure with posteriror distribution of parameters (1:1?)
#2: Think about how to calculate the interspecies synchrony
#3: 


#Libraries
library(dplyr)
library(ggplot2)
library(rstan)
library(tidyr)
library(bayesplot)
library(tidyverse)
library(posterior)


library(purrr)


options(mc.cores = parallel::detectCores())




#Setting working directory
getwd()
setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

list.files()
list.files("Data/")
list.files("Rscript_AnalysisEleonore/")

#Adding Mike's code
setwd("Rscript_AnalysisEleonore/")
util <- new.env()
source("mcmc_visualization_tools.R", local=util)
source("mcmc_analysis_tools_rstan.R", local=util)
ls(util)

setwd("C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting")

# Data import + modification -------------------------------------------------------------
seed_data<-read.csv("Data/SeedData_all.csv")

# 1) Keeping species I'm interested in :
seed_data<-seed_data %>%
  filter(spp %in% c("ABAM","ABLA","CANO","PSME","TSHE","TSME","THPL"))
unique(seed_data$stand)


# 2) I want only possible combination of stands (ie. for CANO I don't want it at low elevation stands)
# a) Extra stand ( to add manually as there is no map/records of the species present at SPRY, PARA and CANO as it was not installed by the PSP crew but Janneke)
# Add SPRY, SUNR, PARA for CANO, ABAM, ABLA
extra_stands <- seed_data %>%
  filter(spp %in% c("CANO", "ABAM", "ABLA"),
         stand %in% c("SPRY", "SUNR", "PARA")) %>%
  filter(year != 2009) %>%
  group_by(spp, stand, year) %>%
  summarise(y    = sum(total_viable_sds, na.rm = TRUE),area = sum(size, na.rm = TRUE),.groups = "drop") %>%
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
  rename(spp = species,
         stand = stands)

seed_filtered <- seed_data %>%
  semi_join(stands_long, by = c("spp", "stand"))

str(seed_filtered)
#Seeds_filtered contains now all the species x stand combination 

#d)Stand elevation
stand_elev <- data.frame(
  stand = c("TO11","TO04","TA01","AV02","AE10","TB13",
            "AO03","AG05","AV06","AX15","AB08","PP17",
            "AV14","AM16","AR07","PARA","SPRY","SUNR"),
  elevation = c(600,668,700,850,1450,850,
                900,950,1060,1090,1100,1150,
                1150,1200,1450,1600,1700,1800)
)
saveRDS(stand_elev, file= "data/stand_elevation_table.rds")


# Data prep for Stan ---------------------------------------------------------------

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


# Model diagnostic --------------------------------------------------------

summary(fit_75)
print(fit_75)

diagnostics <- util$extract_hmc_diagnostics(fit_75)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit_75)

names(samples)
base_samples1 <- util$filter_expectands(samples,c("rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2"))
base_samples2 <- util$filter_expectands(samples,c(paste0("alpha_theta1_species_nc[",1:S, "]"),paste0("alpha_theta2_species_nc[",   1:S, "]"), "sigma_theta1_species","sigma_theta2_species"))
base_samples3 <- util$filter_expectands(samples,c("mu_log_low",paste0('alpha_low_species[',1:S, ']'),paste0('alpha_stand[',   1:S*N_stands, ']'), "sigma_low_species","sigma_stand"))
base_samples4<-util$filter_expectands(samples,c("mu_log_delta",paste0('log_delta_species_nc[',1:S, ']')),"sigma_log_delta_species","sigma_delta_stand")
base_samples5 <- util$filter_expectands(samples,c(paste0('log_phi_species[',1:S, ']'),"log_phi_state[1]","log_phi_state[2]"))
util$check_all_expectand_diagnostics(base_samples1);util$check_all_expectand_diagnostics(base_samples2)

base_samples <- util$filter_expectands(c(paste0('alpha_theta1_species_nc[',1:S, ']'),paste0('alpha_theta2_species_nc[',   1:S, ']'),"rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2",'sigma_theta1_species','sigma_theta2_species','grand_logit_theta1','grand_logit_theta2',"mu_log_low"))

base_samples <- util$filter_expectands(c(paste0('alpha_theta1_species_nc[',1:S, ']'),paste0('alpha_theta2_species_nc[',   1:S, ']'),"rho[1]", "rho[2]","grand_logit_theta1","grand_logit_theta2",'sigma_theta1_species','sigma_theta2_species','grand_logit_theta1','grand_logit_theta2','mu_log_low'))


model_parameters<- c("rho","grand_logit_theta1", "grand_logit_theta2","sigma_theta1_species", "sigma_theta2_species","mu_log_low","alpha_low_species", "sigma_low_species", "sigma_stand", "mu_log_delta","log_delta_species","sigma_log_delta_species","sigma_delta_stand","log_phi_species","log_phi_state")

summary_model<-rstan::summary(fit_75, pars= model_parameters, probs = c(0.025,0.25,0.5,0.75,0.975))$summary

model_table_summary<-as.data.frame(summary_model)%>%
  tibble::rownames_to_column("parameter")%>%
  mutate(across(c(mean, se_mean, sd, '2.5%','25%','50%','75%', '97.5%'), ~round (.x, 3)), 
         n_eff= round(n_eff),
        Rhat= round (Rhat,3))

model_table_summary

write.table(model_table_summary, file = 'data/Tables/model_table_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")



# PPC ---------------------------------------------------------------------

#Overall
# ppc for all 
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")
samples<-util$extract_expectand_vals(fit_75)

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

#Single observation example
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
par(mfrow = c(4,6))
names_yrep <- paste0("y_rep[", 1:stan_data_all$N, "]")
samples<-util$extract_expectand_vals(fit_75)

unique(stand_year_all$spp); unique(stand_year_all$stand)
species_name<-unique(stand_year_all$spp)
stand_name<-unique(stand_year_all$stand)

for (f in 1:length(stan_data_all$start_idxs)) {
  
  start_id <- stan_data_all$start_idxs[f]
  end_id   <- stan_data_all$end_idxs[f]
  
  species_id<-as.character(stand_year_all$spp[start_idxs[f]])
  stand_id<-as.character(stand_year_all$stand[start_idxs[f]])
  
  util$plot_disc_pushforward_quantiles(samples, paste0("state[", start_id:end_id, "]"),
                                       display_ylim = c(1, 2))
  
  title (main = paste(species_id, stand_id))
  
  util$plot_disc_pushforward_quantiles(samples, paste0("y_rep[", start_id:end_id, "]"))
  points(x = 1:length(start_id:end_id), y = stan_data_all$y[start_id:end_id], pch = 20)
  
  readline(prompt = paste0("Stand ", f, "/", length(stan_data_all$start_idxs), " -- Press [Enter] for next..."))
}



# Prior vs posterior ------------------------------------------------

samples <- rstan::extract(fit_75)
n <- nrow(samples$mu_log_low)

rtnorm_function<-function(n,mean, sd){
  qnorm(runif(n,pnorm(0,mean,sd),1), mean, sd)
}

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

#Reshape to long 
prior_long <- priors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(type = "Prior")

post_long <- posteriors %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(type = "Posterior")

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
ggsave(filename= "prior_posterior_plot.png", plot= prior_posterior_plot,path = "Plots/Manuscript")


species<- c("ABAM" ,"ABLA" ,"CANO" ,"PSME", "THPL" ,"TSHE" ,"TSME")
state_name<-c("non-masting", "masting")

prior_phi_sp<-data.frame(
  species= rep(rep(species, each=n),times=2),
  state= rep(state_name, each= n*length(species)),
  value= exp(rnorm(n*length(species) * 2, log(4), 0.6) + rnorm(n*length(species)*2, 0, 0.3)),
  type= "Prior"
)

posterior_phi_sp<-data.frame(
  species= rep(rep(species, each=n),times=2),
  state= rep(state_name, each= n*length(species)),
  value= c(exp(as.vector(samples$log_phi_species + samples$log_phi_state[,1])),exp(as.vector(samples$log_phi_species + samples$log_phi_state[,2]))),
  type= "Posterior"
)

df_phi<-bind_rows(prior_phi_sp, posterior_phi_sp)%>%
  mutate(type= factor(type, levels = c("Prior", "Posterior")), state= factor(state, levels= state_names))

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

species_stand<-stand_year_all%>%
  distinct(spp, stand)%>%
  group_by(spp)%>%
  summarise(stands= paste(stand, collapse = ","), .groups = "drop")
write.table(species_stand, file = 'data/Tables/species_stand.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

samples<-rstan::extract(fit_75)
species<- c("ABAM" ,"ABLA" ,"CANO" ,"PSME", "THPL" ,"TSHE" ,"TSME")
species_id<-stan_data_all$species_id
baseline_area<-min(stan_data_all$area)

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
  
  results[[s]]<-bind_rows(
    summary_masting(low, "Non-masting"),
    summary_masting(high, "Masting"),
    summary_masting(high/low, "Fold change"),
    summary_masting(theta1, "P stay low"), 
    summary_masting(theta2, "P stay high"),
    summary_masting((1-theta1)/((1-theta1)+(1-theta2)), "P time high"),
    summary_masting(phi_low, "Dispersion low"),
    summary_masting(phi_high, "Dispersion high")
  ) %>%
    mutate(species= species[s])
}

species_table<-bind_rows(results)%>% mutate (across(c(median, lo90, hi90), ~round (.x,2)))
species_table

write.table(species_table, file = 'data/Tables/species_table_model_results.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")



#checking the model seed production results vs real observed data

sp_check<-"ABAM"
str(stan_data_all)
str(samples)
selection<-which(stand_year_all$spp==sp_check)
length(selection)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
dim(y_model)
unique(samples$y_rep)

set.seed(1)
d<-(sample(nrow(y_model),200))
abam_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,]))+ labs(title= sp_check, x= "Seed count (log scale)", y= "Density")  
abam_observed_model
ggsave(filename= "abam_observed_model.png", plot= abam_observed_model,path = "Plots/Manuscript")

sp_check<-"ABLA"
selection<-which(stand_year_all$spp==sp_check)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
set.seed(1)
d<-(sample(nrow(y_model),200))
abla_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
abla_observed_model
ggsave(filename= "abla_observed_model.png", plot= abla_observed_model,path = "Plots/Manuscript") 


sp_check<-"CANO"
selection<-which(stand_year_all$spp==sp_check)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
set.seed(1)
d<-(sample(nrow(y_model),200))
cano_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
cano_observed_model
ggsave(filename= "cano_observed_model.png", plot= cano_observed_model,path = "Plots/Manuscript") 


sp_check<-"PSME"
selection<-which(stand_year_all$spp==sp_check)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
set.seed(1)
d<-(sample(nrow(y_model),200))
psme_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
psme_observed_model
ggsave(filename= "psme_observed_model.png", plot= psme_observed_model,path = "Plots/Manuscript")


sp_check<-"THPL"
selection<-which(stand_year_all$spp==sp_check)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
set.seed(1)
d<-(sample(nrow(y_model),200))
thpl_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
thpl_observed_model
ggsave(filename= "thpl_observed_model.png", plot= thpl_observed_model,path = "Plots/Manuscript")


sp_check<-"TSHE"
selection<-which(stand_year_all$spp==sp_check)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
set.seed(1)
d<-(sample(nrow(y_model),200))
tshe_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,]))+ labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
tshe_observed_model
ggsave(filename= "tshe_observed_model.png", plot= tshe_observed_model,path = "Plots/Manuscript")

sp_check<-"TSME"
selection<-which(stand_year_all$spp==sp_check)
y_obs<-stan_data_all$y[selection]
y_model<-samples$y_rep[,selection]
set.seed(1)
d<-(sample(nrow(y_model),200))
tsme_observed_model<-ppc_dens_overlay(log1p(y_obs),log1p(y_model[d,])) + labs(title= sp_check, x= "Seed count (log scale)", y= "Density") 
tsme_observed_model
ggsave(filename= "tsme_observed_model.png", plot= tsme_observed_model,path = "Plots/Manuscript")


mcmc_trace(as.array(fit_75))

# Synchrony functions-------------------------------------------------


#pulling ths state for each draws (one row per posterior draw and one colum per stand/year/specie row)
state_draws <- rstan::extract(fit_75, pars = "state")$state

meta_all <- stand_year_all %>% #I named it meta as it is overall data but can be named as we want
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
  state_counts <- table(df$state) #how many rows fall in each state
  max_count<-max(state_counts) #the amount of counts in the majority state
  major_state_ID<-names(state_counts)[state_counts==max_count] # the most frequent state ID for that year and that draw
  #majority_proportion = max(state_counts) / sum(state_counts) # the proportion of rows in that majority state
  tibble (major_state = paste(major_state_ID, collapse= "/"), proportion= max_count / sum(state_counts))
}


summarise_draws <- function(df, value_col, group_col) { #for example later with ABAM I will have df= abam draws, value_col= proportion and group_col= year
  df %>%
    group_by(.data[[group_col]]) %>% #will split the group by the variable I like with ABAM year for example.
    summarise(
      median = median(.data[[value_col]], na.rm = TRUE), #the median of proportions across all draws that year
      mean = mean(.data[[value_col]], na.rm= TRUE),
      lo90   = quantile(.data[[value_col]], 0.05, na.rm = TRUE), #the 5th percentile of proportion across draws
      hi90   = quantile(.data[[value_col]], 0.95, na.rm = TRUE), #the 95th percentile
      majority_value = names(table(.data[["major_state"]]))[which.max(table(.data[["major_state"]]))], #summarizes the major state across draws
      conf = max(table(.data[["major_state"]]))/n(), #the share of draws in which the winning state was the majority
      .groups = "drop"
    )
}


# Intraspecfic ------------------------------------------------------------


#ABAM
target_species <- "ABAM"
abam_meta <- meta_all %>%
  filter(species == target_species)

abam_state_draws <- state_draws[, abam_meta$row_id]

n_draws_total <-nrow(state_draws) #or abam state draws but I keep it like this so I can reuse it for the next species
draw_idx_use<- seq(1,n_draws_total, by=4)#so it is faster

#mapdfr from the purr package: It takes list and vector and applies a function to each element combining afterwards everything intot a single dataframe. 
abam_draws <- map_dfr(draw_idx_use, function(d) {
  get_unit_year_state(d) %>%
    filter(species == target_species) %>%
    group_by(year) %>%
    filter(n() >= 2) %>%
    group_modify(~majority_state(.x)) %>%
    ungroup()%>%
    mutate(draw = d)
})

names(abam_draws)

abam_summary <- summarise_draws(abam_draws, "proportion", "year")

write.table(abam_summary, file = 'data/Tables/abam_summary.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

abam_majority<-ggplot(abam_summary, aes(x = year, y = median))+
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+
  labs(
    x = NULL, y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "ABAM"
  ) +
  theme_minimal(base_size = 13)
abam_majority
ggsave(filename= "abam_majority.png", plot= abam_majority,path = "Plots/Manuscript")


year_stands_abam<-abam_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_abam, file = 'data/Tables/year_stands_abam.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

abam_stand_state<-abam_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop= FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)

abam_majority_standsID<-abam_stand_state%>%
  left_join(abam_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")

write.table(abam_majority_standsID, file = 'data/Tables/abam_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")


#ABLA
target_species <- "ABLA"
target_meta <- meta_all %>%
  filter(species == target_species)

year_stands_abla<-target_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_abla, file = 'data/Tables/year_stands_abla.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

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

abla_majority<-ggplot(abla_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+
  labs(
    x = "year", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "ABLA"
  ) +
  theme_minimal(base_size = 13)
abla_majority
ggsave(filename= "abla_majority.png", plot= abla_majority,path = "Plots/Manuscript")

abla_summary%>%select(year,median,majority_value)

target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)

abla_majority_standsID<-target_stand_state%>%
  left_join(abla_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")

write.table(abla_majority_standsID, file = 'data/Tables/abla_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#CANO
target_species <- "CANO"
target_meta <- meta_all %>%
  filter(species == target_species)

year_stands_cano<-target_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_cano, file = 'data/Tables/year_stands_cano.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

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

cano_majority<-ggplot(cano_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+
  labs(
    x = "year", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "CANO"
  ) +
  theme_minimal(base_size = 13)
cano_majority
ggsave(filename= "cano_majority.png", plot= cano_majority,path = "Plots/Manuscript")

cano_summary%>%select(year,median,majority_value)

target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)

cano_majority_standsID<-target_stand_state%>%
  left_join(cano_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")

write.table(cano_majority_standsID, file = 'data/Tables/cano_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#PSME
target_species <- "PSME"
target_meta <- meta_all %>%
  filter(species == target_species)

year_stands_psme<-target_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_psme, file = 'data/Tables/year_stands_psme.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

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

psme_majority<-ggplot(psme_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = "year", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "PSME"
  ) +
  theme_minimal(base_size = 13)
psme_majority
ggsave(filename= "psme_majority.png", plot= psme_majority,path = "Plots/Manuscript")

psme_summary%>%select(year,median,majority_value)

target_stand_state<-target_meta%>%
  mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
  select(year, stand, mean_state, state)

psme_majority_standsID<-target_stand_state%>%
  left_join(psme_summary%>%select(year,majority_value), by= "year")%>%
  filter(as.character(state)==majority_value)%>%
  group_by(year,majority_value)%>%
  summarise(n_in_majority=n(),
            stands= paste(sort(stand),collapse = ","),.groups = "drop")

write.table(psme_majority_standsID, file = 'data/Tables/psme_majority_standsID.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

#TSHE
target_species <- "TSHE"
target_meta <- meta_all %>%
  filter(species == target_species)

year_stands_tshe<-target_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_tshe, file = 'data/Tables/year_stands_tshe.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

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

tshe_majority<-ggplot(tshe_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+
  labs(
    x = "year", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "TSHE"
  ) +
  theme_minimal(base_size = 13)
tshe_majority
ggsave(filename= "tshe_majority.png", plot= tshe_majority,path = "Plots/Manuscript")

tshe_summary%>%select(year,median,majority_value)

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

#TSME
target_species <- "TSME"
target_meta <- meta_all %>%
  filter(species == target_species)

year_stands_tsme<-target_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_tsme, file = 'data/Tables/year_stands_tsme.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

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

tsme_majority<-ggplot(tsme_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = "year", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "TSME"
  ) +
  theme_minimal(base_size = 13)
tsme_majority
ggsave(filename= "tsme_majority.png", plot= tsme_majority,path = "Plots/Manuscript")

tsme_summary%>%select(year,median,majority_value)

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

#THPL
target_species <- "THPL"
target_meta <- meta_all %>%
  filter(species == target_species)

year_stands_thpl<-target_meta%>%
  group_by(year)%>%
  summarise(n_stands=n_distinct(stand),
            stands = paste(sort(unique(stand)), collapse = ","))

write.table(year_stands_thpl, file = 'data/Tables/year_stands_thpl.txt', col.names = TRUE,
            row.names = FALSE, sep = "\t")

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

thpl_majority<-ggplot(thpl_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  labs(
    x = "year", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony",
    subtitle = "THPL"
  ) +
  theme_minimal(base_size = 13)
thpl_majority

thpl_summary%>%select(year,median,majority_value)


ggsave(filename= "thpl_majority.png", plot= thpl_majority,path = "Plots/Manuscript")

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

species_order<-c("THPL" , "TSHE", "PSME", "ABAM","CANO" ,"TSME" ,"ABLA")
all_species_intraspecific_synchrony_table<-all_species_intraspecific_synchrony_table%>%
  mutate(species= factor(species, levels= species_order))

all_species_overall_intra_plot<-ggplot(all_species_intraspecific_synchrony_table, aes(x = species, y = mean)) +
  geom_boxplot(fill= "lightblue",alpha = 0.7) +
  geom_jitter()+
  stat_summary(fun= mean, geom= "point", shape= 21, fill= "blue", size=3.5) +
  labs(
    x = "Species", y = "Proportion of stands in the majority state",
    title = "Intraspecific synchrony"
  ) +
  theme_minimal(base_size = 13)
all_species_overall_intra_plot

ggsave(filename= "all_species_overall_intra_plot.png", plot= all_species_overall_intra_plot,path = "Plots/Manuscript")


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


#Interspecific synchrony

draw_idx_use<- seq(1,n_draws_total, by=4)#so it is faster
n_draws_total <-nrow(state_draws)

summarise_draws2 <- function(df, value_col, group_col) { #for example later with ABAM I will have df= abam draws, value_col= proportion and group_col= year
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





#overall per years
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
ggsave(filename= "overall_synchrony.png", plot= overall_synchrony,path = "Plots/Manuscript")

#interspecific synchrony. 
#1. defining for each species and year the majority state across all stands

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

species_year_summary_plot<-ggplot(species_year_summary, aes(x = year, y = median)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
  facet_wrap(~species)+
  labs(
    x = "year", y = "Proportion in the majority state",
    title= "Intraspecific synchrony"
  ) +
  theme_minimal(base_size = 13)
species_year_summary_plot
ggsave(filename= "species_year_summary_plot.png", plot= species_year_summary_plot,path = "Plots/Manuscript")
  
#2.Per year, the majortiy state across species (species all together)

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




#Per stands 
#for one stand
unique(stand_year_all$stand)
#AB08 (one stand example, same as for the species)
# target_stand <- "AB08"
# target_meta <- meta_all %>%
#   filter(stand == target_stand)
# 
# target_state_draws <- state_draws[, target_meta$row_id]
# 
# n_draws_total <-nrow(target_state_draws) 
# 
# ab08_draws <- map_dfr(draw_idx_use, function(d) {
#   get_unit_year_state(d) %>%
#     filter(stand == target_stand) %>%
#     group_by(year) %>%
#     filter(n() >= 2) %>%
#     group_modify(~majority_state(.x)) %>%
#     mutate(draw = d)
# })
# 
# ab08_summary <- summarise_draws(ab08_draws, "proportion", "year")
# write.table(ab08_summary, file = 'data/Tables/ab08_summary.txt', col.names = TRUE,
#             row.names = FALSE, sep = "\t")
# 
# ab08_majority<-ggplot(ab08_summary, aes(x = year, y = median)) +
#   geom_line(alpha = 0.7) +
#   geom_point(size = 1.8) +
#   geom_ribbon(aes(ymin=lo90,ymax=hi90),alpha=0.4, fill="lightblue")+ 
#   labs(
#     x = "year", y = "Proportion of species in the majority state",
#     title = "Interspecific synchrony",
#     subtitle = "AB08"
#   ) +
#   theme_minimal(base_size = 13)
# ab08_majority
# 
# ab08_summary%>%select(year,median,majority_value)
# 
# 
# ggsave(filename= "ab08_majority.png", plot= ab08_majority,path = "Plots/Manuscript")
# 
# target_species_state<-target_meta%>%
#   mutate(mean_state= colMeans(state_draws[draw_idx_use,row_id, drop = FALSE]), state= round(mean_state))%>%
#   select(year, species, mean_state, state)
# 
# ab08_majority_speciesID<-target_species_state%>%
#   left_join(ab08_summary%>%select(year,majority_value), by= "year")%>%
#   filter(as.character(state)==majority_value)%>%
#   group_by(year,majority_value)%>%
#   summarise(n_in_majority=n(),
#             species= paste(sort(species),collapse = ","),.groups = "drop")
# 
# write.table(ab08_majority_speciesID, file = 'data/Tables/ab08_majority_speciesID.txt', col.names = TRUE,
#             row.names = FALSE, sep = "\t")


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
stand_results<-map(all_stands, stand_loop)%>%
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




