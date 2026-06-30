library(tidyverse)
source("R/base_functions.R")
source("R/bayesian_model_functions.R")

# Load data scaffold
data_scaffold <- read_csv("results/synthetic_data_experiments/data_scaffold.csv") 

RVE_model <- fully_flexible_RVE_model

synth_experiment_iterations <- 4000
synth_experiment_chains <- 4

# Vaccine strain names to use in simulations. In some years, equivalent strains
# with different names are listed. We take just one vaccine strain per subtype
# per year
core_vaccine_strains <- vaccine_strains %>%
  group_by(year, subtype) %>%
  slice(1) %>%
  ungroup() %>%
  filter(year != 1)

# Default value of parameters shared across strains
universal_params <- 
  tibble(
    omega = 0.005,
    tau = 35,  # Baseline tau (irrelevant with include_kinetics = F, as are other kinetic parameters)
    f0 = 0.3,  # Baseline long-term fraction
    r_LT = -0.5,  # Repeat vaccination effect on long-term fraction
    b_LT = 0, # Effect of Hpre on long-term fraction
    k_peak_flu = 1,
    k_peak_cov2_M = 0, # effect for sarscov2 infection, male
    k_peak_cov2_F = 0, # effect for sarscov2 infection, female
    sigma_HAI = 1.2, # Measurement error for HAI
    sigma_FRNT = 1, # Measurement error for FRNT
    sigma_upeak_shared = 0.5, # Standard deviation of shared individual effects on peak titer
    sigma_upeak_subtype = 0, # Standard deviation of subtype-specific individual effects on peak titer
    sigma_upeak_year = 0,
    beta_M = 0, # Effect of male sex on Hpeak
    beta_age3140 = -0.5, # Effect of age group 31-40 on Hpeak
    beta_age4150 = -0.7,  # Effect of age group 41-50 on Hpeak
    beta_smoking = -1,     # Effect of smoking
    beta_asthma = 0,       # Effect of asthma
    beta_BMI_under = 0,    # Effect for underweight BMI relative to healthy
    beta_BMI_over = -0.6,      # Effect for overweight BMI relative healthy,
    gamma = 0, # Decay rate for RVE
    beta_time_under28 = 0, # Effect of time-since-vax group under 28 days on Hpeak
    beta_time_over42 = 0   # Effect of time-since-vax group over 42 days on Hpeak
  )

# Default true values of strain-specific parameters
strain_specific_params <- core_vaccine_strains %>%
  select(strain) %>%
  unique() %>%
  mutate(a = case_match(strain,
    "B/Washington/02/2019" ~ 4.04,
    "B/Phuket/3073/2013" ~ 6.5,
    "A/Wisconsin/588/2019" ~ 5.15,
    "A/Cambodia/e0826360/2020" ~ 5,
    "B/Austria/1359417/2021" ~ 4.72,
    "A/Darwin/9/2021" ~ 5,
    "A/Wisconsin/67/2022" ~ 5.8
  )) %>%
  mutate(b_peak = case_match(strain,
    "B/Washington/02/2019" ~ 0.35,
    "B/Phuket/3073/2013" ~ 0.65,
    "A/Wisconsin/588/2019" ~ 0.4,
    "A/Cambodia/e0826360/2020" ~ 0.4,
    "B/Austria/1359417/2021" ~ 0.35,
    "A/Darwin/9/2021" ~ 0.45,
    "A/Wisconsin/67/2022" ~ 0.5
  )) %>% 
  filter(strain %in% unique(data_scaffold$strain))

# Default values of repeat vaccination effect parameters
RVE_params <- RVE_model$form %>%
  select(RVE_level, default_synthetic_value) %>%
  rename(r = default_synthetic_value) %>%
  # This unique needs to be here for models other than the fully flexible.
  unique()
  
params <- list(universal = universal_params, strain_specific = strain_specific_params,
              RVE = RVE_params)

censoring = T

include_prevax_error_synth_data <- T
include_kinetics_synth_data <- F 
include_RVE <- T
scaffold_multiples <- 1 # How many times to replicate the scaffold for generating synthetic data
include_upeak_shared <- T
include_upeak_subtype <- F
include_upeak_year <- F
include_k_peak_flu_synth_data <- T
include_k_peak_cov2 <- F
include_NAI_infections_synth_data <- F # By default, we simulate/fit effects of PCR-confirmed flu infections only.
include_time_since_vax_synth_data <- F
include_sex_synth_data <- T
include_age_synth_data <- T
include_BMI_synth_data <- T
include_smoking_synth_data <- T
include_asthma_synth_data <- T

if(any(universal_params$k_peak_cov2_F != 0 | universal_params$k_peak_cov2_M != 0)){
  if(!include_k_peak_cov2){
    stop("Simulation assumes a non-zero effect of prior covid infection, but inference won't include that effect")
  }
}

 # Generate a new realization of the synthetic data each time, instead of running replicate searches
 # on a single fixed synthetic data set.
fixed_synth_data <- F
