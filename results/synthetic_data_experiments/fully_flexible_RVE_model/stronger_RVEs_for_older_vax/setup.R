source("R/default_synthetic_setup.R")


# For this scenario, we make it so that for any given subtype/observation year
# more recent vaccines have a *weaker* effect on the current response than older vaccines
# We achieve that by reordering the default values.

RVE_model <- fully_flexible_RVE_model

RVE_model$form <- RVE_model$form %>%
    group_by(subtype, year) %>%
    arrange(subtype, year, previous_year, .by_group = TRUE) %>%
    mutate(r = sort(default_synthetic_value, decreasing = FALSE)) %>%
    group_by(subtype) %>%
    # (For year 2, since there are no prior effects to reorder, we just set 
    #  the effect of year 1 to the minimum across all years for the subtype)
    mutate(r = if_else(year == 2, max(r), r)) %>%
    ungroup()

RVE_params <- RVE_model$form %>%
  select(RVE_level, r) %>%
  unique()

params$RVE <- RVE_params
