source("R/default_synthetic_setup.R")

RVE_model <- no_decay_model

params$RVE <- RVE_model$form %>%
  select(RVE_level, default_synthetic_value) %>%
  rename(r = default_synthetic_value) %>%
  unique()