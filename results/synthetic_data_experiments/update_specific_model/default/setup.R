source("R/default_synthetic_setup.R")

RVE_model <- update_specific_model

params$RVE <- RVE_model$form %>%
  select(RVE_level, default_synthetic_value) %>%
  rename(r = default_synthetic_value) %>%
  unique()

params$universal$gamma <- -0.5
