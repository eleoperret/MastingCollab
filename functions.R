
# 1) DataCleaning ------------------------------------------------------------
filter_species_data <- function(seed_data, species){
  
  data <- seed_data %>%
    filter(spp == species)
  
  return(data)
  
}

remove_invalid_stands <- function(data, stands){
  
  data <- data %>%
    filter(!stand %in% stands)
  
  return(data)
  
}


remove_invalid_stands <- function(data, species_name, stands_per_species){
  # Get the valid stands for this species
  valid_stands <- stands_per_species$stands[stands_per_species$species == species_name][[1]]
  
  # Filter data to keep only valid stands
  data <- data %>%
    filter(stand %in% valid_stands)
  
  return(data)
}

# 2) DatasetCreation ---------------------------------------------------------

##Removing each time the year 2009

create_stand_year_dataset <- function(data){
  
  stand_year_df <- data %>%
    group_by(stand, year) %>%
    filter(year!="2009") %>%
    summarise(
      y = sum(total_viable_sds, na.rm = TRUE),
      area = sum(size, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(stand, year)
  
  return(stand_year_df)
  
}





# 3) Plots --------------------------------------------------------



plot_total_seeds_per_stand <- function(species_data){
  
  # Summarize total seeds per stand and join elevation
  total_stand <- species_data %>%
    group_by(stand) %>%
    summarise(
      total_seeds = sum(total_viable_sds, na.rm = TRUE),
      elevation = first(elevation),
      .groups = "drop"
    )
  
  # Plot
  p <- ggplot(total_stand, aes(x = reorder(stand, elevation), y = total_seeds)) +
    geom_col(fill = "steelblue") +
    theme_bw() +
    labs(
      title = paste("Total viable seeds per stand:"),
      x = "Stand",
      y = "Total viable seeds"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  print(p)
} 

#Seed density
plot_seed_density_timeseries <- function(df){
  
  df <- df %>%
    mutate(seed_density = y / area)
  
  p <- ggplot(df,
              aes(x = year,
                  y = seed_density,
                  color = stand,
                  group = stand)) +
    geom_line(size = 1) +
    geom_point() +
    theme_minimal()
  
  print(p)
  
}


# 4) StanDataPrep ------------------------------------------------------------


#stan data preparation
prepare_stan_data <- function(stand_year_df){
  
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
  
  return(stan_data)
  
}





# 5) ModelFit ----------------------------------------------------------------

#model fit
fit_hmm_model <- function(stan_file, stan_data){
  
  mod <- stan_model(stan_file)
  
  fit <- sampling(
    mod,
    data = stan_data,
    chains = 4,
    cores = 4,
    iter = 2000,
    warmup = 1000
  )
  
  return(fit)
  
}