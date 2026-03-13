#Function for plotting the posterior predictive checks

run_hmm_diagnostics <- function(fit, stan_data, stand_year_df, stand_elev, F,
                                species = species,
                                base_dir = "C:/Users/eperret/polybox - Eleonore Perret (eleonore.perret@usys.ethz.ch)@polybox.ethz.ch/phD/PhD/R/Masting_UBC/Masting/Plots") {
  

save_dir <- file.path(base_dir, paste0("plots_", species))
if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
message("All plots will be saved in: ", normalizePath(save_dir))
  
util <- new.env()
source("mcmc_analysis_tools_rstan.R", local = util)
source("mcmc_visualization_tools.R", local = util)
  

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
  

samples <- util$extract_expectand_vals(fit)
  
base_samples <- util$filter_expectands(samples,
    c("rho","theta1","theta2","log_lambda","log_mu","stand_effect_raw","phi1","phi2","sigma"),check_arrays = TRUE)

util$check_all_expectand_diagnostics(base_samples)

all_plots <- list()
  
names_yrep <- paste0("y_rep[", 1:stan_data$N, "]")
  
#PPC plots 
#1) 
png(file.path(save_dir, paste0("PPC_all_50.png")), width = 800, height = 600)
util$plot_hist_quantiles(
    samples, "y_rep",
    bin_min = 0,
    bin_max = max(stan_data$y) * 1.1, # dynamic max range
    bin_delta = 50,
    baseline_values = stan_data$y,
    main = paste0(species, " - PPC Model (bin 50)")
  )
dev.off()

#2)
png(file.path(save_dir, paste0("PPC_all_100.png")), width = 800, height = 600)
util$plot_hist_quantiles(
    samples, "y_rep",
    bin_min = 0,
    bin_max = max(stan_data$y) * 1.1,
    bin_delta = 100,
    baseline_values = stan_data$y,
    main = paste0(species, " - PPC Model (bin 100)")
  )
dev.off()
  

#Ordering the stands in the right order
original_order <- unique(stand_year_df$stand)
stand_order <- stand_elev %>%
    filter(stand %in% original_order) %>%
    arrange(elevation) %>%
    pull(stand)
perm <- match(stand_order, original_order)


#for PPC for each stands
for(s in 1:stan_data$F){
  png(file.path(save_dir, paste0("PPC_", stand_order[s], ".png")), width = 800, height = 600)
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  bin_max <- max(stan_data$y[idxs]) * 1.1
  subsamples <-  util$filter_expectands(samples,
                                        paste0("y_rep[", idxs, "]"))
  util$plot_hist_quantiles(subsamples, "y_rep",
                           bin_min = 0, bin_max = bin_max, bin_delta = 50,
                           baseline_values = stan_data$y[idxs])
  all_plots[[paste0("PPC_", stand_order[s])]] <- recordPlot()  # <--- store plot
  dev.off()
}

#masting vs non-masting
for(i in 1:stan_data$F){
  s <- perm[i]
  idxs <- stan_data$start_idxs[s]:stan_data$end_idxs[s]
  names <- paste0("y_rep[", idxs, "]")
  
  png(file.path(save_dir, paste0("MastVsNonMast_", stand_order[s], ".png")), width = 800, height = 600)
  
  util$plot_disc_pushforward_quantiles(
    samples, names,
    baseline_values = stan_data$y[idxs],
    main = paste0(stand_order[s], " (", stand_elev$elevation[stand_elev$stand == stand_order[s]], " m)")
  )
  
  all_plots[[paste0("MastVsNonMast_", stand_order[s])]] <- recordPlot()
  dev.off()
}
  
# rho
png(file.path(save_dir, "rho.png"), width = 800, height = 400)
par(mfrow = c(1,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["rho[1]"]],50,
                                  flim=c(0,1),
                                  display_name=bquote(.(species) ~ rho * "(Low)"))
util$plot_expectand_pushforward(samples[["rho[2]"]],50,
                                  flim=c(0,1),
                                  display_name=bquote(.(species) ~ rho * "(High)"))
dev.off()
  
# theta
png(file.path(save_dir, "theta.png"), width = 800, height = 400)
par(mfrow = c(1,2), mar = c(5,5,1,1))
util$plot_expectand_pushforward(samples[["theta1[1]"]],50,
                                  flim=c(0,1),
                                  display_name=bquote(.(species) ~ theta * "(Staying in Low)"))
util$plot_expectand_pushforward(samples[["theta2[1]"]],50,
                                  flim=c(0,1),
                                  display_name=bquote(.(species) ~ theta * "(Staying in Mast)"))
dev.off()
  
# log_lambda
png(file.path(save_dir, "log_lambda.png"), width = 600, height = 400)
par(mfrow = c(1,1), mar = c(5,5,1,1))
util$plot_expectand_pushforward(exp(samples[["log_lambda"]]),20,
                                  display_name=bquote(.(species) ~ lambda * "(low)"))
dev.off()
  
# log_mu
png(file.path(save_dir, "log_mu.png"), width = 600, height = 400)
util$plot_expectand_pushforward(exp(samples[["log_mu"]]),50,
                                  display_name=bquote(.(species) ~ mu * "(high)"))
dev.off()


#log_alpha  
for(f in 1:F){
  # Get the stand name for labeling and saving
  stand_name <- stand_order[f]
  
  # Prepare PNG file for this stand
  png(file.path(save_dir, paste0("log_alpha_", stand_name, ".png")), width = 600, height = 400)
  
  # Plot the posterior predictive distribution for log_alpha
  util$plot_expectand_pushforward(
    exp(samples[[paste0("log_alpha[", f, "]")]]),  # remove exp() if you want log scale
    50,                                             # number of bins
    display_name = bquote(alpha * "(" * .(stand_name) * ")")  # nicely labeled
  )
  
  # Save plot to list
  all_plots[[paste0("log_alpha_", stand_name)]] <- recordPlot()
  
  dev.off()
}



  
# phi
png(file.path(save_dir, "phi.png"), width = 600, height = 400)
par(mfrow=c(2,1))
util$plot_expectand_pushforward(samples[["phi1"]],50,
                                  display_name=bquote(.(species) ~ phi * "(NB nonmast)"))
util$plot_expectand_pushforward(samples[["phi2"]],50,
                                  display_name=bquote(.(species) ~ phi * "(NB mast)"))
dev.off()

  
message("All posterior predictive and parameter plots saved for species: ", species)
return(all_plots)
}



