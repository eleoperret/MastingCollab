#Code for prior vs posterior predctive check

# Prior vs Posterior check ------------------------------------------------
plot_prior_posterior <- function(prior_df, fit, stan_data,
                                 save_plot = FALSE,
                                 filename = "prior_posterior_plot.png"){
  
  F <- stan_data$F
  
  post_samples <- rstan::extract(fit)
  
  posterior_df <- data.frame(
    log_lambda = post_samples$log_lambda,
    log_mu = post_samples$log_mu,
    sigma = post_samples$sigma,
    phi1 = post_samples$phi1,
    phi2 = post_samples$phi2
  )
  
  # add stand effects
  
  for(f in 1:Stand){
    stand_raw <- rnorm(n_prior,0,1)
    prior_df[[paste0("stand_effect_", f)]] <- stand_raw * prior_df$sigma
  }
  
  
  
  prior_long <- reshape(
    prior_df,
    direction = "long",
    varying = list(names(prior_df)),
    v.names = "value",
    timevar = "parameter",
    times = names(prior_df),
    idvar = "draw"
  )
  
  posterior_long <- reshape(
    posterior_df,
    direction = "long",
    varying = list(names(posterior_df)),
    v.names = "value",
    timevar = "parameter",
    times = names(posterior_df),
    idvar = "draw"
  )
  
  
  # Assume species is a variable
  species <- species
  
  p <- ggplot() +
    geom_density(data = prior_long,
                 aes(x = value, colour = "Prior"),
                 linewidth = 0.8) +
    geom_density(data = posterior_long,
                 aes(x = value, colour = "Posterior"),
                 linewidth = 0.5) +
    facet_wrap(~parameter, scales = "free_x") +
    labs(
      title = paste0("Prior vs Posterior - ", species),
      x = "Parameter value",
      y = "Density",
      color = "Curve"
    ) +
    theme_minimal()
  
  print(p)
  
}
