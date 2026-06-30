library(tidyverse)

source("R/base_functions.R")
source("R/bayesian_model_functions.R")
source("R/load_data.R")

dir.create("results/bayesian_fits_real_data", showWarnings = F)
dir.create("results/synthetic_data_experiments", showWarnings = F)

processed_data <- titers %>%
  flag_vaccine_strains() %>%
  filter(is_vaccine_strain) %>%
  # Prioritizing HAI batches that measured the same person in multiple years
  select_titer_batches() %>%
  prepare_input_data_bayesian()

# Currently there's no BMI > 30 in the titer data. 
# So we can't fit a parameter for the effect of BMI >= 30+
# Revise this categorization (and the stan code) if that ever changes
stopifnot(all(processed_data$BMI < 30))

processed_data <- processed_data %>%
  mutate(BMI_group = case_when(
    BMI < 18.5 ~ "under",
    BMI >= 18.5 & BMI < 25 ~ "healthy",
    BMI >= 25 & BMI < 30 ~ "over"
  ))

time_since_vax_group <- processed_data %>%
  select(individual, year, timepoint, t) %>%
  unique() %>%
  mutate(time_since_vax_group = case_when(
    t < 28  ~ "under28",
    t >= 28 & t <= 42 ~ "28to42",
    t > 42  ~ "over42"
  )) 

processed_data <- left_join(processed_data, time_since_vax_group)

# Use titer data to create a scaffold for simulating synthetic data
data_scaffold <- processed_data %>%
  select(sample_id, measurement_replicate, year, timepoint, t, time_since_vax_group, individual, treatment, sex,
         age_group, asthma, smoking, BMI, BMI_group, n_previous_vax, matches("vax_Y"),
         strain, subtype, titer_type, undetectable_value, infection_before_sample_subtype_matched_PCR,
         infection_before_sample_subtype_matched_NAI, recent_infection_before_sample_cov2_PCR, titer_type, Hpre)

# Exporting the scaffold (for synthetic data experiments)
write_csv(data_scaffold %>%
            select(-matches("RVE_level")),
          "results/synthetic_data_experiments/data_scaffold.csv")

# Export the processed real data 
# (Make sure there are no dates)
stopifnot(any(str_detect(names(processed_data), "date"))) 
write_csv(processed_data, "results/bayesian_fits_real_data/processed_data.csv")
