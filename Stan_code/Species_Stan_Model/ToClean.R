

















# Checking the 0's --------------------------------------------------------

# --- Extract from Stan fit ---
# y_rep: matrix [n_draws x N]
# state: matrix [n_draws x N]
# y, sp, stand_id, series_id: your original data vectors


draws <- rstan::extract(fit_allPhi3)

y_rep <- draws$y_rep   # matrix [draws x N]
state <- draws$state   # matrix [draws x N]

n_draws <- nrow(y_rep)
N      <- ncol(y_rep)

# Expand stand_id from series-level to observation-level
stand_obs <- integer(stan_data_all$N)
for (f in seq_len(stan_data_all$F)) {
  stand_obs[start_idxs[f]:end_idxs[f]] <- years_per_series$stand_id[f]
}

obs <- tibble(
  t     = 1:stan_data_all$N,
  y     = stan_data_all$y,
  sp    = stan_data_all$sp,
  stand = stand_obs
)

# ── 1. Overall zero frequency ──────────────────────────────────────────────
prop_zero_rep <- rowMeans(y_rep == 0)   # one value per draw

tibble(prop_zero_rep) |>
  ggplot(aes(prop_zero_rep)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = mean(stan_data_all$y == 0), color = "red", linewidth = 1) +
  labs(title = "PPC: overall zero frequency",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()

# ── 2. Zero frequency by species ───────────────────────────────────────────
species_list <- sort(unique(stan_data_all$sp))

ppc_zeros_species <- map_dfr(species_list, function(s) {
  idx <- which(stan_data_all$sp == s)
  rep_props <- rowMeans(y_rep[, idx, drop = FALSE] == 0)
  tibble(
    species      = s,
    prop_zero    = rep_props,
    obs_prop_zero = mean(stan_data_all$y[idx] == 0)
  )
})

ppc_zeros_species |>
  ggplot(aes(prop_zero)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  geom_vline(aes(xintercept = obs_prop_zero), color = "red", linewidth = 0.8) +
  facet_wrap(~species, scales = "free_x") +
  labs(title = "PPC: zero frequency by species",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()

# ── 3. Zero frequency by inferred state ────────────────────────────────────
# Use posterior mode state per observation (or loop over draws for full uncertainty)

# Quick version: use most frequent state per observation
state_mode <- apply(state, 2, \(x) as.integer(names(which.max(table(x)))))

ppc_zeros_state <- map_dfr(1:2, function(k) {
  idx <- which(state_mode == k)
  if (length(idx) == 0) return(NULL)
  rep_props <- rowMeans(y_rep[, idx, drop = FALSE] == 0)
  tibble(
    state         = paste("State", k),
    prop_zero     = rep_props,
    obs_prop_zero = mean(stan_data_all$y[idx] == 0)
  )
})

ppc_zeros_state |>
  ggplot(aes(prop_zero)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  geom_vline(aes(xintercept = obs_prop_zero), color = "red", linewidth = 0.8) +
  facet_wrap(~state) +
  labs(title = "PPC: zero frequency by inferred state",
       subtitle = "Red = observed. Excess zeros in State 1 → consider hurdle",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()

# ── 4. Zero frequency by species x state ───────────────────────────────────
ppc_zeros_sp_state <- map_dfr(species_list, function(s) {
  map_dfr(1:2, function(k) {
    idx <- which(stan_data_all$sp == s & state_mode == k)
    if (length(idx) < 5) return(NULL)   # skip sparse combos
    rep_props <- rowMeans(y_rep[, idx, drop = FALSE] == 0)
    tibble(
      species       = s,
      state         = paste("State", k),
      prop_zero     = rep_props,
      obs_prop_zero = mean(stan_data_all$y[idx] == 0)
    )
  })
})

ppc_zeros_sp_state |>
  ggplot(aes(prop_zero)) +
  geom_histogram(bins = 25, fill = "steelblue", alpha = 0.8) +
  geom_vline(aes(xintercept = obs_prop_zero), color = "red", linewidth = 0.8) +
  facet_grid(species ~ state, scales = "free_x") +
  labs(title = "PPC: zero frequency by species and state",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()



# Checking improvement ----------------------------------------------------


draws <- rstan::extract(fit_allPhi3)

# Compare phi by species and state
phi_df <- expand.grid(species = 1:6, state = 1:2) |>
  mutate(
    mean_phi = mapply(function(s, k) mean(exp(draws$log_phi_species_state[, s, k])),
                      species, state),
    lower = mapply(function(s, k) quantile(exp(draws$log_phi_species_state[, s, k]), 0.1),
                   species, state),
    upper = mapply(function(s, k) quantile(exp(draws$log_phi_species_state[, s, k]), 0.9),
                   species, state),
    sp_name = c("ABAM","ABLA","CANO","PSME","THPL","TSHE")[species],
    state = paste("State", state)
  )

ggplot(phi_df, aes(sp_name, mean_phi, color = state)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(0.4), width = 0.2) +
  labs(title = "Phi by species and state",
       y = "Dispersion (phi)", x = "Species") +
  theme_minimal()


state_old <- apply(rstan::extract(fit_allPhi2)$state, 2, 
                   \(x) mean(x == 2))  # prob of state 2 per observation
state_new <- apply(rstan::extract(fit_allPhi3)$state, 2, 
                   \(x) mean(x == 2))

# Quick visual comparison for a problematic series
par(mfrow = c(2,1))
plot(state_old[start_idxs[9]:end_idxs[9]], type = "l", 
     main = "Old model - P(State 2) Series 9", ylim = c(0,1))
plot(state_new[start_idxs[9]:end_idxs[9]], type = "l", 
     main = "New model - P(State 2) Series 9", ylim = c(0,1))


draws <- rstan::extract(fit_allPhi3)

delta_df <- expand.grid(species = 1:6) |>
  mutate(
    sp_name = c("ABAM","ABLA","CANO","PSME","THPL","TSHE"),
    mean_delta = sapply(1:6, \(s) mean(exp(draws$log_delta_high_species[,s]))),
    lower = sapply(1:6, \(s) quantile(exp(draws$log_delta_high_species[,s]), 0.1)),
    upper = sapply(1:6, \(s) quantile(exp(draws$log_delta_high_species[,s]), 0.9))
  )

ggplot(delta_df, aes(sp_name, mean_delta)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(title = "Estimated delta (high - low state) by species",
       y = "exp(log_delta_high_species)", x = "Species") +
  theme_minimal()




library(tidyverse)

# Add species names to your observation data
sp_names <- c("ABAM", "ABLA", "CANO", "PSME", "THPL", "TSHE")

obs <- tibble(
  y      = stan_data_all$y,
  sp     = stan_data_all$sp,
  sp_name = sp_names[stan_data_all$sp],
  stand  = stand_obs,  # your expanded stand vector from earlier
  area   = stan_data_all$area
) |>
  mutate(y_per_area = y / area)  # standardize by trap area

# ── 1. Raw time series per stand and species ───────────────────────────────
# First need year index - assuming your data has a year column
obs <- obs |>
  mutate(year = stand_year_all$year)  # adjust to your actual year column name

ggplot(obs, aes(x = year, y = y_per_area, group = stand)) +
  geom_line(alpha = 0.4, color = "steelblue") +
  geom_point(alpha = 0.6, size = 0.8) +
  facet_wrap(~sp_name, scales = "free_y") +
  labs(title = "Raw seed counts per unit area by species",
       subtitle = "Each line = one stand",
       x = "Year", y = "Seeds / area") +
  theme_minimal()

# ── 2. Median across stands per year ──────────────────────────────────────
obs_summary <- obs |>
  group_by(sp_name, year) |>
  summarise(
    median_y  = median(y_per_area),
    q10       = quantile(y_per_area, 0.1),
    q90       = quantile(y_per_area, 0.9),
    .groups = "drop"
  )

ggplot(obs_summary, aes(x = year)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = median_y), color = "steelblue") +
  geom_point(aes(y = median_y), size = 1.5) +
  facet_wrap(~sp_name, scales = "free_y") +
  labs(title = "Seed production over time by species",
       subtitle = "Line = median across stands, ribbon = 10th-90th percentile",
       x = "Year", y = "Seeds / area") +
  theme_minimal()

# ── 3. Coefficient of variation per species ────────────────────────────────
# High CV = high year-to-year variability = likely mast cycling
obs |>
  group_by(sp_name, year) |>
  summarise(mean_y = mean(y_per_area), .groups = "drop") |>
  group_by(sp_name) |>
  summarise(
    cv = sd(mean_y) / mean(mean_y),
    mean = mean(mean_y),
    .groups = "drop"
  ) |>
  ggplot(aes(sp_name, cv)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(title = "Coefficient of variation by species",
       subtitle = "Higher CV = more year-to-year variability = stronger mast signal",
       x = "Species", y = "CV") +
  theme_minimal()



draws <- rstan::extract(fit_new)

# Total delta = grand mean + species effect
delta_total <- sweep(draws$log_delta_high_species, 1, 
                     draws$log_delta_high_grand_mean, "+")

delta_df <- tibble(
  sp_name = rep(sp_names, each = nrow(draws$log_delta_high_species)),
  log_delta = as.vector(delta_total)
) |>
  group_by(sp_name) |>
  summarise(
    mean_delta = mean(exp(log_delta)),
    lower = quantile(exp(log_delta), 0.1),
    upper = quantile(exp(log_delta), 0.9)
  )

ggplot(delta_df, aes(sp_name, mean_delta)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(title = "Total estimated delta by species (grand mean + species effect)",
       y = "exp(log_delta_high_grand_mean + log_delta_high_species)",
       x = "Species") +
  theme_minimal()



# General diagnostic plots: -----------------------------------------------

# ============================================================
# HMM Seed Production Model — Diagnostics & Synchrony Analysis
# ============================================================
# Requires: rstan, tidyverse, bayesplot, patchwork

library(tidyverse)
library(rstan)
library(bayesplot)
library(patchwork)

# Colour palette (consistent throughout)
clr_state1  <- "#6BAED6"   # blue  = low state
clr_state2  <- "#E6550D"   # orange = mast state
clr_obs     <- "black"
sp_names    <- c("ABAM", "ABLA", "CANO", "PSME", "THPL", "TSHE")

# Pull draws once
draws      <- rstan::extract(fit_allPhi3)
y_rep      <- draws$y_rep          # [draws x N]
state_draws <- draws$state         # [draws x N]
y_obs      <- stan_data_all$y
sp_vec     <- stan_data_all$sp
N          <- stan_data_all$N
S          <- stan_data_all$S
F_         <- stan_data_all$F

# Posterior mode state per observation
state_mode <- apply(state_draws, 2, \(x) as.integer(names(which.max(table(x)))))

# ── Helper: expand stand_id to observation level ──────────────────────────
stand_obs <- integer(N)
for (f in seq_len(F_)) {
  stand_obs[start_idxs[f]:end_idxs[f]] <- stan_data_all$stand_id[f]
}

obs <- tibble(
  t        = seq_len(N),
  y        = y_obs,
  sp       = sp_vec,
  sp_name  = sp_names[sp_vec],
  stand    = stand_obs,
  year     = stand_year_all$year,   # adjust if column name differs
  area     = stan_data_all$area,
  state    = state_mode
)

# ============================================================
# 1.  PRIOR vs POSTERIOR
# ============================================================
n_draws <- nrow(draws$log_phi_species_state)
n_prior <- 4000

prior_vs_posterior <- function(posterior_vec, prior_vec, param_name) {
  tibble(
    value  = c(posterior_vec, prior_vec),
    source = rep(c("Posterior", "Prior"), c(length(posterior_vec), length(prior_vec)))
  ) |>
    ggplot(aes(value, fill = source, colour = source)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    scale_fill_manual(values  = c(Posterior = "#2171B5", Prior = "#BDBDBD")) +
    scale_colour_manual(values = c(Posterior = "#2171B5", Prior = "#636363")) +
    labs(title = param_name, x = NULL, y = "Density") +
    theme_minimal(base_size = 11) +
    theme(legend.title = element_blank(),
          plot.title   = element_text(size = 10, face = "bold"))
}

# Grand means
p_pvp1 <- prior_vs_posterior(
  draws$grand_logit_theta1,
  rnorm(n_prior, 1, 0.7),
  "grand_logit_theta1"
)
p_pvp2 <- prior_vs_posterior(
  draws$grand_logit_theta2,
  rnorm(n_prior, 0, 0.7),
  "grand_logit_theta2"
)
p_pvp3 <- prior_vs_posterior(
  draws$grand_mean_low,
  rnorm(n_prior, 2.6, 1.0),
  "grand_mean_low"
)
p_pvp4 <- prior_vs_posterior(
  draws$log_delta_high_grand_mean,
  rnorm(n_prior, 3, 1),
  "log_delta_high_grand_mean"
)

# Sigma hyperparameters
p_pvp5 <- prior_vs_posterior(
  draws$sigma_theta1_species,
  abs(rnorm(n_prior, 0, 0.7)),
  "sigma_theta1_species"
)
p_pvp6 <- prior_vs_posterior(
  draws$sigma_low_species,
  abs(rnorm(n_prior, 0, 0.5)),
  "sigma_low_species"
)
p_pvp7 <- prior_vs_posterior(
  draws$mu_log_phi,
  rnorm(n_prior, log(4), 0.6),
  "mu_log_phi"
)
p_pvp8 <- prior_vs_posterior(
  draws$sigma_log_phi,
  abs(rnorm(n_prior, 0, 0.5)),
  "sigma_log_phi"
)

pvp_plot <- (p_pvp1 | p_pvp2 | p_pvp3 | p_pvp4) /
  (p_pvp5 | p_pvp6 | p_pvp7 | p_pvp8) +
  plot_annotation(
    title    = "Prior vs Posterior — key parameters",
    subtitle = "Separation indicates data are informative",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave("diag_01_prior_vs_posterior.pdf", pvp_plot,
       width = 14, height = 7, device = "pdf")

# ============================================================
# 2.  PPC — OVERALL COUNT DISTRIBUTION
# ============================================================
ppc_dist <- ppc_dens_overlay(
  y     = y_obs,
  yrep  = y_rep[sample(nrow(y_rep), 100), ]
) +
  coord_cartesian(xlim = c(0, quantile(y_obs, 0.995))) +
  labs(title = "PPC: count distribution (observed vs 100 posterior draws)") +
  theme_minimal()

ggsave("diag_02_ppc_count_dist.pdf", ppc_dist,
       width = 8, height = 5, device = "pdf")

# ============================================================
# 3.  PPC — ZERO FREQUENCY (overall, by species, by state)
# ============================================================
prop_zero_rep <- rowMeans(y_rep == 0)

p_z1 <- tibble(prop_zero_rep) |>
  ggplot(aes(prop_zero_rep)) +
  geom_histogram(bins = 40, fill = "#6BAED6", alpha = 0.85) +
  geom_vline(xintercept = mean(y_obs == 0), colour = "red", linewidth = 0.9) +
  labs(title = "PPC: overall zero frequency",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()

# By species
ppc_zero_sp <- map_dfr(seq_len(S), function(s) {
  idx <- which(sp_vec == s)
  tibble(
    species       = sp_names[s],
    prop_zero     = rowMeans(y_rep[, idx, drop = FALSE] == 0),
    obs_prop_zero = mean(y_obs[idx] == 0)
  )
})

p_z2 <- ppc_zero_sp |>
  ggplot(aes(prop_zero)) +
  geom_histogram(bins = 30, fill = "#6BAED6", alpha = 0.85) +
  geom_vline(aes(xintercept = obs_prop_zero), colour = "red", linewidth = 0.8) +
  facet_wrap(~species, scales = "free_x") +
  labs(title = "PPC: zero frequency by species",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()

# By state
ppc_zero_state <- map_dfr(1:2, function(k) {
  idx <- which(state_mode == k)
  tibble(
    state         = paste("State", k),
    prop_zero     = rowMeans(y_rep[, idx, drop = FALSE] == 0),
    obs_prop_zero = mean(y_obs[idx] == 0)
  )
})

p_z3 <- ppc_zero_state |>
  ggplot(aes(prop_zero)) +
  geom_histogram(bins = 30, fill = "#6BAED6", alpha = 0.85) +
  geom_vline(aes(xintercept = obs_prop_zero), colour = "red", linewidth = 0.8) +
  facet_wrap(~state) +
  labs(title = "PPC: zero frequency by state",
       subtitle = "Red = observed",
       x = "Proportion zeros", y = "Draws") +
  theme_minimal()

ggsave("diag_03_ppc_zeros.pdf",
       p_z1 / p_z2 / p_z3,
       width = 10, height = 14, device = "pdf")

# ============================================================
# 4.  PPC — MEAN AND SD BY SPECIES
# ============================================================
ppc_mean_sp <- map_dfr(seq_len(S), function(s) {
  idx <- which(sp_vec == s)
  tibble(
    species   = sp_names[s],
    stat      = "Mean",
    rep_val   = rowMeans(y_rep[, idx, drop = FALSE]),
    obs_val   = mean(y_obs[idx])
  )
}) |>
  bind_rows(
    map_dfr(seq_len(S), function(s) {
      idx <- which(sp_vec == s)
      tibble(
        species  = sp_names[s],
        stat     = "SD",
        rep_val  = apply(y_rep[, idx, drop = FALSE], 1, sd),
        obs_val  = sd(y_obs[idx])
      )
    })
  )

p_meansd <- ppc_mean_sp |>
  ggplot(aes(rep_val)) +
  geom_histogram(bins = 30, fill = "#6BAED6", alpha = 0.85) +
  geom_vline(aes(xintercept = obs_val), colour = "red", linewidth = 0.8) +
  facet_grid(stat ~ species, scales = "free") +
  labs(title = "PPC: mean and SD by species",
       subtitle = "Red = observed",
       x = "Value", y = "Draws") +
  theme_minimal()

ggsave("diag_04_ppc_mean_sd.pdf", p_meansd,
       width = 14, height = 6, device = "pdf")

# ============================================================
# 5.  STATE IDENTIFIABILITY — P(State 2) per series
# ============================================================
state_prob_series <- map_dfr(seq_len(F_), function(f) {
  idx   <- start_idxs[f]:end_idxs[f]
  s_id  <- sp_vec[start_idxs[f]]
  tibble(
    series   = f,
    sp_name  = sp_names[s_id],
    year     = obs$year[idx],
    p_state2 = colMeans(state_draws[, idx, drop = FALSE] == 2)
  )
})

p_state_id <- state_prob_series |>
  ggplot(aes(year, p_state2)) +
  geom_line(colour = clr_state2, linewidth = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  facet_wrap(~sp_name + series, scales = "free_x",
             labeller = label_wrap_gen(width = 15)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title    = "State identifiability: P(mast state) over time",
       subtitle = "Values near 0 or 1 indicate confident state assignments",
       x = "Year", y = "P(State 2 = mast)") +
  theme_minimal(base_size = 9)

ggsave("diag_05_state_identifiability.pdf", p_state_id,
       width = 16, height = 20, device = "pdf")

# ============================================================
# 6.  STATE SEPARATION — total delta by species
# ============================================================
delta_total <- sweep(draws$log_delta_high_species, 1,
                     draws$log_delta_high_grand_mean, "+")

delta_df <- tibble(
  sp_name    = rep(sp_names, each = nrow(delta_total)),
  log_delta  = as.vector(delta_total)
) |>
  group_by(sp_name) |>
  summarise(
    mean_delta = mean(exp(log_delta)),
    lower      = quantile(exp(log_delta), 0.1),
    upper      = quantile(exp(log_delta), 0.9),
    .groups    = "drop"
  )

p_delta <- delta_df |>
  ggplot(aes(fct_reorder(sp_name, mean_delta), mean_delta)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8) +
  coord_flip() +
  labs(title    = "State separation: total delta (high − low) by species",
       subtitle = "80% credible interval",
       x = NULL, y = "Mast state excess seeds (delta)") +
  theme_minimal()

ggsave("diag_06_state_separation.pdf", p_delta,
       width = 7, height = 4, device = "pdf")

# ============================================================
# 7.  TRANSITION PROBABILITIES by species
# ============================================================
theta_df <- map_dfr(seq_len(S), function(s) {
  # Average theta across series for this species
  f_idx <- which(sp_vec[start_idxs] == s)
  tibble(
    sp_name = sp_names[s],
    theta1  = mean(draws$theta1[, f_idx]),   # P(stay low | low)
    theta2  = mean(draws$theta2[, f_idx]),   # P(stay mast | mast)
    theta1_lo = quantile(rowMeans(draws$theta1[, f_idx, drop=FALSE]), 0.1),
    theta1_hi = quantile(rowMeans(draws$theta1[, f_idx, drop=FALSE]), 0.9),
    theta2_lo = quantile(rowMeans(draws$theta2[, f_idx, drop=FALSE]), 0.1),
    theta2_hi = quantile(rowMeans(draws$theta2[, f_idx, drop=FALSE]), 0.9)
  )
})

p_theta <- theta_df |>
  pivot_longer(c(theta1, theta2),
               names_to = "param",
               values_to = "value") |>
  mutate(
    label = if_else(param == "theta1",
                    "P(stay low | low)",
                    "P(stay mast | mast)"),
    lo = if_else(param == "theta1", theta1_lo, theta2_lo),
    hi = if_else(param == "theta1", theta1_hi, theta2_hi)
  ) |>
  ggplot(aes(sp_name, value, colour = label)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(0.4), size = 0.8) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_manual(values = c(clr_state1, clr_state2)) +
  labs(title    = "Transition probabilities by species",
       subtitle = "80% credible interval",
       x = NULL, y = "Probability", colour = NULL) +
  theme_minimal()

ggsave("diag_07_transition_probs.pdf", p_theta,
       width = 8, height = 5, device = "pdf")

# ============================================================
# 8.  SYNCHRONY ANALYSIS
# ============================================================

# -- 8a. Pairwise stand synchrony: proportion of years both stands
#        are in the same state
stand_year_state <- obs |>
  select(year, stand, sp_name, state) |>
  distinct()

# For each species, compute fraction of years stands share same state
synchrony_between_stands <- stand_year_state |>
  group_by(sp_name, year) |>
  summarise(
    n_stands    = n(),
    n_mast      = sum(state == 2),
    prop_mast   = mean(state == 2),
    # synchrony = proportion of stand pairs in same state
    sync        = (choose(n_mast, 2) + choose(n_stands - n_mast, 2)) /
      choose(n_stands, 2),
    .groups = "drop"
  ) |>
  filter(n_stands > 1)  # need at least 2 stands

# Annual synchrony across species
p_sync_annual <- synchrony_between_stands |>
  ggplot(aes(year, sync, colour = sp_name)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_brewer(palette = "Dark2") +
  facet_wrap(~sp_name) +
  labs(title    = "Between-stand synchrony over time by species",
       subtitle = "Proportion of stand pairs in same state each year",
       x = "Year", y = "Synchrony", colour = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("diag_08a_synchrony_annual.pdf", p_sync_annual,
       width = 10, height = 7, device = "pdf")

# -- 8b. Overall synchrony summary by species
sync_summary <- synchrony_between_stands |>
  group_by(sp_name) |>
  summarise(
    mean_sync = mean(sync),
    lo        = quantile(sync, 0.25),
    hi        = quantile(sync, 0.75),
    .groups   = "drop"
  )

p_sync_summary <- sync_summary |>
  ggplot(aes(fct_reorder(sp_name, mean_sync), mean_sync)) +
  geom_col(fill = "#6BAED6", alpha = 0.85, width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  coord_flip() +
  labs(title    = "Mean between-stand synchrony by species",
       x = NULL, y = "Synchrony") +
  theme_minimal()

ggsave("diag_08b_synchrony_summary.pdf", p_sync_summary,
       width = 7, height = 4, device = "pdf")

# -- 8c. Cross-species synchrony — are mast years shared across species?
# For each year and stand, how many species are in mast state?
cross_sp_sync <- obs |>
  group_by(year, stand) |>
  summarise(
    n_sp_mast = sum(state == 2),
    n_sp_obs  = n(),
    prop_mast = n_sp_mast / n_sp_obs,
    .groups = "drop"
  ) |>
  group_by(year) |>
  summarise(
    mean_prop_mast = mean(prop_mast),
    lo = quantile(prop_mast, 0.25),
    hi = quantile(prop_mast, 0.75),
    .groups = "drop"
  )

p_cross_sp <- cross_sp_sync |>
  ggplot(aes(year, mean_prop_mast)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = clr_state2, alpha = 0.25) +
  geom_line(colour = clr_state2, linewidth = 1) +
  geom_point(colour = clr_state2, size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title    = "Cross-species synchrony over time",
       subtitle = "Mean proportion of species in mast state per year",
       x = "Year", y = "Proportion of species in mast state") +
  theme_minimal()

ggsave("diag_08c_cross_species_synchrony.pdf", p_cross_sp,
       width = 8, height = 4, device = "pdf")

# ============================================================
# 9.  POSTERIOR UNCERTAINTY ON STATE PROBABILITY — with raw data
# ============================================================
# Pick a few representative series to show to colleagues

plot_series_with_states <- function(f_id) {
  idx    <- start_idxs[f_id]:end_idxs[f_id]
  s_name <- sp_names[sp_vec[start_idxs[f_id]]]
  st_id  <- stan_data_all$stand_id[f_id]
  
  p_state <- tibble(
    year     = obs$year[idx],
    p_mast   = colMeans(state_draws[, idx, drop = FALSE] == 2),
    lo       = apply(state_draws[, idx, drop=FALSE] == 2, 2,
                     quantile, probs = 0.1),
    hi       = apply(state_draws[, idx, drop=FALSE] == 2, 2,
                     quantile, probs = 0.9)
  ) |>
    ggplot(aes(year, p_mast)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = clr_state2, alpha = 0.2) +
    geom_line(colour = clr_state2) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    scale_y_continuous(limits = c(0,1), labels = scales::percent) +
    labs(y = "P(mast)", x = NULL,
         title = glue::glue("{s_name} — Stand {st_id}")) +
    theme_minimal(base_size = 10)
  
  p_counts <- tibble(
    year  = obs$year[idx],
    y     = obs$y[idx],
    state = state_mode[idx]
  ) |>
    ggplot(aes(year, y, colour = factor(state))) +
    geom_point(size = 2.5) +
    geom_line(aes(group = 1), colour = "grey70", linewidth = 0.4) +
    scale_colour_manual(values = c("1" = clr_state1, "2" = clr_state2),
                        labels = c("1" = "Low", "2" = "Mast")) +
    labs(y = "Seed count", x = "Year", colour = "State") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  
  p_state / p_counts
}

# Save one per species (first series of each)
first_series_per_sp <- map_int(seq_len(S), function(s) {
  which(sp_vec[start_idxs] == s)[1]
})

series_plots <- map(first_series_per_sp, plot_series_with_states)
combined_series <- wrap_plots(series_plots, ncol = 2)

ggsave("diag_09_series_state_examples.pdf", combined_series,
       width = 12, height = 18, device = "pdf")

# ============================================================
# 10.  DISPERSION (phi) BY SPECIES AND STATE
# ============================================================
phi_df <- expand.grid(species = seq_len(S), state = 1:2) |>
  as_tibble() |>
  mutate(
    sp_name    = sp_names[species],
    mean_phi   = mapply(\(s,k) mean(exp(draws$log_phi_species_state[,s,k])),
                        species, state),
    lower      = mapply(\(s,k) quantile(exp(draws$log_phi_species_state[,s,k]),0.1),
                        species, state),
    upper      = mapply(\(s,k) quantile(exp(draws$log_phi_species_state[,s,k]),0.9),
                        species, state),
    state_label = paste("State", state)
  )

p_phi <- phi_df |>
  ggplot(aes(sp_name, mean_phi, colour = state_label)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(0.4), size = 0.8) +
  scale_colour_manual(values = c("State 1" = clr_state1,
                                 "State 2" = clr_state2)) +
  labs(title    = "Overdispersion (phi) by species and state",
       subtitle = "80% credible interval",
       x = NULL, y = "phi", colour = NULL) +
  theme_minimal()

ggsave("diag_10_phi_by_species_state.pdf", p_phi,
       width = 8, height = 4, device = "pdf")

# ============================================================
# DONE — list output files
# ============================================================
message("\n✓ All diagnostics saved:")
message("  diag_01_prior_vs_posterior.pdf")
message("  diag_02_ppc_count_dist.pdf")
message("  diag_03_ppc_zeros.pdf")
message("  diag_04_ppc_mean_sd.pdf")
message("  diag_05_state_identifiability.pdf")
message("  diag_06_state_separation.pdf")
message("  diag_07_transition_probs.pdf")
message("  diag_08a_synchrony_annual.pdf")
message("  diag_08b_synchrony_summary.pdf")
message("  diag_08c_cross_species_synchrony.pdf")
message("  diag_09_series_state_examples.pdf")
message("  diag_10_phi_by_species_state.pdf")
