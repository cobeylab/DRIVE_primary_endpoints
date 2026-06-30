source("R/default_synthetic_setup.R")

# This uses the default model in default_synthetic_setup 
# (changed so that r = 0) to simulate without RVEs
params$RVE <- params$RVE %>% mutate(r = 0)


# This makes it so that repeat vaccination effects are constrained
# at 0 and not estimated in stan.
include_RVE <- F