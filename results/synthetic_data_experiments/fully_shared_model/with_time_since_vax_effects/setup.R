source("R/default_synthetic_setup.R")

RVE_model <- fully_shared_model

params$RVE <- RVE_model$form %>%
  select(RVE_level, default_synthetic_value) %>%
  rename(r = default_synthetic_value) %>%
  unique()

params$universal$gamma <- -0.5 
params$universal$beta_time_under28 <- -0.5
params$universal$beta_time_over42 <- 0.3
include_time_since_vax_synth_data <- T
