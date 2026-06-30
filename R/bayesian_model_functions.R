library(GGally)
library(ggh4x)
library(boot) # has logit and inv.logit functions.
source("R/base_functions.R")

# Parameter names for a model with longitudinal titer dynamics,
# (not used in this manuscript)
kinetic_par_names <- c("tau", "f0", "omega", "r_LT", "b_LT")

assign_wide_format_vax_indicators <- function(data){
  # For each year, annotate with vaccination status in years 1-3 in wide format
  # (This will include non-sensical year pairs, e.g. a row with year == 2 gets annotations for years 2, 3 and 4)
  # Downstream code prevents current and future years from counting toward repeated vaccination effects
  data <- data %>%
    left_join(
      treatment_assignment %>%
        rename(individual = pID) %>%
        select(individual, matches("^Y[1-3]")) %>%
        mutate(across(matches("^Y[1-3]"), ~.x == "vaccine")) %>%
        rename_with(.fn = ~str_replace(.x, "^Y([1-3])$", "vax_Y\\1"), matches("^Y[1-3]"))
    )
}

assign_RVE_levels <- function(data, RVE_model){
  # Given an input data set, annotates RVE levels together with prior vaccination status in wide format
  # This means adding columns indicating which effect prior vaccinations have in the current year vaccination
  # E.g., for a row with year = 4, annotates which effect applies to the y4 response from vaccines from years 1:3 
  # The wide format is convenient for the Stan fit but includes pairs the current and future years
  # E.g., a row with year = 2 gets annotated with RVE in years 2:3 in wide format, not just year 1.
  # These current- and future-year columns need to be defined such that current and future years don't count toward RVEs
  # Normally we'd assign NA to the RVE levels of future-year vaccinations.
  # Because Stan doesn't allow missing values, instead we assign them a special level (RVE_special_level)
  # In the Stan code, an effect of zero is hard-coded to apply to that level (see prepare_stan_input)

  if("drive" %in% names(data)){
    if(any(data$drive == 2)){
      stop("This function needs to be changed to use calendar years once we're analyzing DRIVE II")
  }
}

# Annotate RVE levels
annotated_data <- data %>%
  mutate(year = factor(year)) %>%
  left_join(
    RVE_model$form %>%
      mutate(year = factor(year)) %>%
      select(year, previous_year, subtype, RVE_level) %>%
      arrange(year, previous_year, subtype) %>%
      pivot_wider(
      id_cols = c(year, subtype),
      names_from = previous_year,
      values_from = RVE_level,
      names_glue = "{.value}_Y{previous_year}"
      ) %>%
      mutate(across(matches("(_Y[1-3])$"), ~replace_na(.x, RVE_special_level))),
    by = c("year", "subtype")
  )

# If fitting this later to data from year 5, revise the function
stopifnot(all(as.integer(data$year) <= 4))

# We must order the levels of the RVE_level columns such that RVE_special_level is the last level.
# (the stan model assigns a zero RVE effect to the last level, preventing current and future year vaccination
# from counting toward repeat vaccination effects)
RVE_levels <- get_RVE_model_levels(RVE_model)
  
annotated_data <- annotated_data %>%
  mutate(
    RVE_level_Y1 = factor(RVE_level_Y1, levels = RVE_levels),
    RVE_level_Y2 = factor(RVE_level_Y2, levels = RVE_levels),
    RVE_level_Y3 = factor(RVE_level_Y3, levels = RVE_levels)
  )
  stopifnot(nrow(annotated_data) == nrow(data))
  return(annotated_data)
}


# Prepare input data for fitting Bayesian models (a separate function prepares the final stan input list)
prepare_input_data_bayesian <- function(data){

  input_data <- data %>%
    # Keep only vaccine strains
    flag_vaccine_strains() %>%
    filter(is_vaccine_strain) %>%
    # Keep only people in the vaccinated groups
    filter(str_detect(treatment, "V")) %>%
    annotate_years_with_vaccine_updates() %>%
    annotate_with_n_prior_vax() %>%
    select(pID, age_group, sex, smoking, asthma, BMI, year,
           timepoint, subtype, strain,
           treatment, matches("titer_type"), matches('n_previous_vax'),
           matches('proportion_IgG'), matches('infection_before'),
           matches('infection_since'), matches('ndays'),
           titer,
           matches("FRNT_assay"))  %>%
    # For individual effects...
    mutate(pID_subtype_combination = paste(pID, subtype, sep = '-'),
           pID_year_combination = paste(pID, year, sep = '-')) %>%
    # Not fitted to year 1 data (when everyone is a first time vaccinee)
    filter(year != 1) %>%
    mutate(year = factor(year))

  # Puts pre-vaccination titer in wide format for each post-vaccination timepoint
  # For HAI titers this is straightforward...
  pre_HAI <- input_data %>% filter(str_detect(titer_type, "HAI"), timepoint == 0) %>%
    rename(prevax_titer = titer) %>%
    select(pID, strain, year, prevax_titer)

  post_HAI <- input_data %>% filter(str_detect(titer_type, "HAI"), timepoint != 0)  %>%
    rename(postvax_titer = titer)

  paired_HAI <- post_HAI %>%
    left_join(pre_HAI) %>%
    select(-matches('assay')) %>%
    mutate(measurement_replicate = 1)

  # For FRNT there are repeated measurements for some samples
  # Put repeated post-vaccination FRNT titers in long format, include them
  # if, and only if, there's a repeated measurement of the pre-vaccination titer
  processed_FRNT <- input_data %>%
    filter(str_detect(titer_type, "FRNT")) %>%
    # So far, titer column will have the GMT for repeated measurements. Remove that column.
    select(-matches("titer")) %>%
    pivot_longer(cols = matches("assay"), names_to = "measurement_replicate", names_prefix = "FRNT_assay_",
                 values_to = "titer") %>%
    mutate(measurement_replicate = as.integer(measurement_replicate))
                
  pre_FRNT <- processed_FRNT %>%
    filter(timepoint == 0) %>%
    rename(prevax_titer = titer) %>%
    select(pID, strain, year, prevax_titer, measurement_replicate) %>%
    filter(!is.na(prevax_titer))

  post_FRNT <- processed_FRNT %>%
    filter(timepoint != 0) %>%
    rename(postvax_titer = titer) %>%
    filter(!is.na(postvax_titer))

  paired_FRNT <- post_FRNT %>%
    left_join(pre_FRNT) %>%
    filter(!is.na(prevax_titer)) %>%
    mutate(titer_type = "FRNT")

  # Combining FRNT and HAI titers back into a single tibble
  paired_titers <- bind_rows(paired_HAI, paired_FRNT) %>%
    mutate(Hpre = log2(prevax_titer),
           Ht = log2(postvax_titer)) %>%
    rename(t = ndays_since_year_vax,
           individual = pID)  %>%
    select(individual, timepoint, subtype,
           strain, matches("titer_type"),
           Ht, Hpre, everything()) %>%
     mutate(undetectable_value =
              case_when(
                str_detect(titer_type, "FRNT") ~ 10,
                str_detect(titer_type, "HAI") ~ 5
  )) 

  # Because we have repeated measurements for the same sample (same person/year/timepoint/strain)
  # We must uniquely link each measurement (row) with a sample
  sample_id <- paired_titers %>%
    select(individual, year, timepoint, strain) %>%
    unique() %>%
    mutate(sample_id = 1:n())

  paired_titers <- paired_titers %>%
    left_join(sample_id) %>%
    select(sample_id, measurement_replicate, everything())

  # Vaccine indicators in wide format
  paired_titers <- paired_titers %>%
    assign_wide_format_vax_indicators()

  return(paired_titers)
}

# Artificially de-censors prevaccination titers (for use in simulation scaffold)
decensor_Hpre <- function(data, sigma_HAI, sigma_FRNT) {
  decensored_prevax_titers <- data %>%
    # For this, take only the first replicate measurement if there's more than one
    filter(measurement_replicate == 1) %>%
    select(individual, strain, year, undetectable_value, Hpre, titer_type) %>%
    unique() %>%
    rowwise() %>%
    mutate(
      model_mean = ifelse(
        near(Hpre, log2(undetectable_value)),
        # If observed prevax titer undetectable,
        # set mean to half the log of the smallest actual dilution.
        0.5 * log2(2 * undetectable_value), 
        # Otherwise, set mean to halfway between observed titer and next possible dilution
        Hpre + 0.5
      ),
      sigma = ifelse(str_detect(titer_type, "HAI"), sigma_HAI, sigma_FRNT),
      Hpre_decensored = rnorm(1, mean = model_mean, sd = sigma)
    ) %>%
    ungroup() %>%
    mutate(Hpre = Hpre_decensored) %>%
    select(-Hpre_decensored, -model_mean, -sigma)

  output_data <- data %>%
    select(-Hpre) %>%
    left_join(decensored_prevax_titers %>%
                select(individual, strain, year, Hpre),
              by = c("individual", "strain", "year"))

  # Check that each individual/strain/year combination (represented as multiple timepoints) was
  # assigned a single decensored Hpre value.
  stopifnot(
    output_data %>%
      group_by(individual, strain, year) %>%
      summarise(n_unique_Hpre = length(unique(Hpre))) %>%
      ungroup() %>%
      pull(n_unique_Hpre) %>%
      unique() == 1
  )
  return(output_data)
}

plot_simulated_responses <- function(simulated_data, params, include_NAI_infections) {

  # If strains have different taus, needs to be revised
  plot_time <- params$universal$tau
  stopifnot(length(unique(plot_time)) == 1) 

  plot_data_in_peak_plot <- nrow(simulated_data %>% filter(t == plot_time)) > 0

  peak_plot <- plot_Ht_vs_Hpre_predictions(params_summary = params, data = simulated_data,
                                           t = plot_time, plot_data = plot_data_in_peak_plot,
                                           include_NAI_infections = include_NAI_infections) +
                  ylab("Peak post-vaccination titer")

  # Add baseline (time 0) data for trajectory plotting
  trajectory_pl_data <- simulated_data %>%
    bind_rows(
      simulated_data %>% 
        select(year, individual, treatment, n_previous_vax, strain, subtype, measurement_replicate, Hpre) %>%
        unique() %>%
        mutate(t = 0, 
               timepoint = 0,
               Ht = Hpre)
    )
  
  # Create trajectory plot showing observed values over time
  trajectory_plot <- ggplot(trajectory_pl_data, aes(x = t, y = Ht, 
                                                    group = interaction(individual, measurement_replicate), color = factor(n_previous_vax))) +
    geom_line(alpha = 0.3, linewidth = 0.5) +  # Thin and semi-transparent
    facet_grid(rows = vars(subtype), cols = vars(year)) +  # Facet by year and subtype
    scale_color_manual(values = scales::colour_ramp(c("lightblue", "darkblue"))(
      seq(0, 1, length.out = length(unique(trajectory_pl_data$n_previous_vax)))
    )) +  # Consistent color mapping
    theme_minimal() +
    labs(x = "Days since vaccination", y = "Ht", 
         title = "Observed Trajectories Over Time", 
         color = "Prior Vaccinations")
  
  return(list(peak_plot = peak_plot, trajectory_plot = trajectory_plot))
}

calculate_RVE <- function(data) {
  data %>%
    # RVE = sum of vax_Yk * r_Yk * (1 + gamma)^(year - k - 1) for k < year of current observation
    # cases where k >= year are handled by r_Y being set to zero.
    # so current and future years don't count toward RVEs.
    mutate(RVE =
      # as.integer(as.character()) safely converts year to integer when it's input as a factor
      vax_Y1 * r_Y1 * (1 + gamma)^(as.integer(as.character(year)) - 1 - 1) +
      vax_Y2 * r_Y2 * (1 + gamma)^(as.integer(as.character(year)) - 2 - 1) +
      vax_Y3 * r_Y3 * (1 + gamma)^(as.integer(as.character(year)) - 3 - 1)
    )
}

calculate_Hpeak <- function(data) {
  data %>%
    mutate(
      Hpeak = a * (1 + RVE) + 
        (1 - b_peak * (1 + RVE)) * Hpre + 
        infection_before_sample_subtype_matched * k_peak_flu +
        recent_infection_before_sample_cov2 * k_peak_cov2 + 
        u_peak_shared +
        u_peak_subtype +
        u_peak_year +
        sex_effect +
        age_effect +
        smoking_effect +
        asthma_effect +
        bmi_effect +
        time_effect
    )
}
 
calculate_HLT <- function(data) {
  data %>%
    mutate(
      logit_f = logit(f0) + r_LT * n_previous_vax + b_LT * Hpre,
      f = inv.logit(logit_f),
      HLT = Hpre + f * (Hpeak - Hpre)
    )
}

calculate_Ht <- function(data) {
  data %>%
    mutate(
      Ht = ifelse(
        t <= tau,
        Hpre + (Hpeak - Hpre) * (t / tau),
        HLT + (Hpeak - HLT) * exp(-omega * (t - tau))
      )
    )
}


# R version of the stan function for calculating censored likelihood.
# (for cross-checking)
log_Pr_H <- function(Ht, H_hat, sigma, L, U) {
  # Special case: Ht = L-1 → P(H_hat < L)
  if (abs(Ht - (L - 1)) < 1e-9) {
    return(pnorm(L, mean = H_hat, sd = sigma, log.p = TRUE))
  }
  # Special case: Ht = U → P(H_hat ≥ U)
  else if (abs(Ht - U) < 1e-9) {
    # note lower.tail = F is what gives us the probability to the right.
    return(pnorm(U, mean = H_hat, sd = sigma, lower.tail = FALSE, log.p = TRUE))
  }
  # General case: P(Ht ≤ H_hat < Ht + 1)
  else {
    # To compute the probability in the interval Ht - Ht + 1, compare 2 approaches

    # Regular (lower tail) CDF
    log_cdf_lower <- pnorm(Ht, mean = H_hat, sd = sigma, log.p = TRUE)
    log_cdf_upper <- pnorm(Ht + 1, mean = H_hat, sd = sigma, log.p = TRUE)

    # Complementary (upper tail) CDF
    log_complementary_cdf_lower <- pnorm(Ht, mean = H_hat, sd = sigma, lower.tail = FALSE, log.p = TRUE)
    log_complementary_cdf_upper <- pnorm(Ht + 1, mean = H_hat, sd = sigma, lower.tail = FALSE, log.p = TRUE)

    # For numerical stability, use the approach whose factors have the largest absolute sum
    # (i.e., where we're less likely to be subtracting very small values)
    sum_regular <- log_cdf_upper + log_cdf_lower
    sum_complementary <- log_complementary_cdf_lower + log_complementary_cdf_upper

    if (abs(sum_regular) > abs(sum_complementary)) {
      result <- log(exp(log_cdf_upper) - exp(log_cdf_lower))
    } else {
      # Notice that here we submit lower - upper
      result <- log(exp(log_complementary_cdf_lower) - exp(log_complementary_cdf_upper))
    }
    return(result)
  }
}

distribute_parameter_values <- function(data, universal_params, strain_specific_params, RVE_params){
  
  data %>%
    # Combine universal parameters (identical for all rows)
    bind_cols(universal_params) %>%
    # Combine strain-specific parameters
    left_join(strain_specific_params, by = c("strain")) %>%
    # Distribute RVE parameters into columns r_Y1:3, linking values of r to the RVE level in each column using RVE_params
    # (assumes data already annotated with assign_RVE_levels)
    rowwise() %>%
    mutate(
      across(
      .cols = starts_with("RVE_level_Y"),
      # Assigns 0 if RVE level is equal to the special level (no RVE from vaccination in current & future years)
      .fns = ~ifelse(.x == RVE_special_level, 0, RVE_params$r[match(.x, RVE_params$RVE_level)]),
      .names = "r_{.col}"
      )
    ) %>%
    rename_with(~str_replace(.x, "r_RVE_level_", "r_"), starts_with("r_RVE_level_")) %>%
    ungroup()

}

simulate_vaccine_responses <- function(data_scaffold, RVE_model, params, censoring, include_kinetics,
                                       include_prevax_measurement_error, include_NAI_infections,
                                       scaffold_multiples) {

  # If !include_kinetics, set tau to 0 (instant peak), f0 to 1, omega, r_LT and bL_t to 0 (no waning)
  # (whatever the values of those parameters were in params)
  if(!include_kinetics){
    params$universal <- params$universal %>%
      mutate(tau = 0, f0 = 1, omega = 0, r_LT = 0, b_LT = 0)
  }
  
  universal_params <- params$universal
  strain_specific_params <- params$strain_specific
  RVE_params <- params$RVE

  sigma_HAI <- universal_params$sigma_HAI
  sigma_FRNT <- universal_params$sigma_FRNT
  sigma_upeak_shared <- universal_params$sigma_upeak_shared
  sigma_upeak_subtype <- universal_params$sigma_upeak_subtype
  sigma_upeak_year <- universal_params$sigma_upeak_year

  # Add vaccination indicator variables and assign RVE levels according to RVE model
  # (must be done before potentially expanding the scaffold)
  data_scaffold <- data_scaffold %>%
    mutate(year = factor(year)) %>%
    assign_RVE_levels(RVE_model)

  # Determine which prior infection columns to use
  # If not including NAI infections, define prior subtype-matched infections based on PCR only
  if(!include_NAI_infections){
    data_scaffold <- data_scaffold %>%
      mutate(infection_before_sample_subtype_matched = infection_before_sample_subtype_matched_PCR)
  }else{
    data_scaffold <- data_scaffold %>%
      mutate(infection_before_sample_subtype_matched = infection_before_sample_subtype_matched_PCR | infection_before_sample_subtype_matched_NAI)
  }

  # COVID infections based on PCR only
  data_scaffold <- data_scaffold %>%
    mutate(recent_infection_before_sample_cov2 = recent_infection_before_sample_cov2_PCR)

  expanded_scaffold <- c()
  for(i in 1:scaffold_multiples){
    expanded_scaffold <- bind_rows(expanded_scaffold,
                                    data_scaffold %>%
                                      mutate(year = factor(year), 
                                             individual = paste0(individual, "_", i)))
  }

  # Assign sex as factor (assume levels are "M" and "F")
  expanded_scaffold <- expanded_scaffold %>%
    mutate(sex = factor(sex, levels = c("M", "F")))

  # Make sure age_group as factor
  expanded_scaffold <- expanded_scaffold %>%
    mutate(age_group = factor(age_group, levels = c("18-30", "31-40", "41-50")))

  # Ensure BMI_group is a factor:
  # (Order of levels doesn't need to match what Stan assumes, but making it so just in case)
  expanded_scaffold <- expanded_scaffold %>%
    mutate(BMI_group = factor(BMI_group, levels = c("healthy", "under", "over")))

  # Decensor Hpre in the scaffold to generate realistic pre-vaccination titer data
  decensored_scaffold <- decensor_Hpre(expanded_scaffold, sigma_HAI, sigma_FRNT)

  # Using the decensored scaffold, compute the mean and sd for Hpre for each combination of strain and n_previous_vax
  # (matches assumption of the Bayesian model)
  Hpre_means_and_sds <- decensored_scaffold %>%
    mutate(strain_npriorvax_group_index = as.integer(factor(paste(strain, n_previous_vax)))) %>%
    select(strain, n_previous_vax, strain_npriorvax_group_index, individual, Hpre) %>%
    unique() %>%
    group_by(strain, n_previous_vax, strain_npriorvax_group_index) %>%
    summarise(Hpre_mean = mean(Hpre),
              Hpre_sd = sd(Hpre)) %>%
    ungroup()
            
  # Replace Hpre's with normally distributed values 
  true_Hpres <- decensored_scaffold %>%
    select(individual, strain, year, n_previous_vax) %>%
    unique() %>%
    left_join(Hpre_means_and_sds %>%
                select(-strain_npriorvax_group_index)) %>%
    rowwise() %>%
    mutate(Hpre = rnorm(1, mean = Hpre_mean, sd = Hpre_sd)) %>%
    ungroup()

  decensored_scaffold <- decensored_scaffold %>%
    select(-Hpre) %>%
    left_join(true_Hpres)

  # For samples with multiple measurements, this checks that replicated measurements
  # will have the same underlying true latent titer
  stopifnot(
    decensored_scaffold %>%
      group_by(sample_id) %>%
      summarise(S = length(unique(Hpre))) %>%
      ungroup() %>%
      pull(S) %>%
      unique() == 1
  )

  # (Process Hpre_means_and_sds for output)
  Hpre_means_and_sds <- Hpre_means_and_sds %>%
    pivot_longer(cols = matches("Hpre"), names_to = "parameter", values_to = "true_value") %>%
    mutate(parameter = paste0(parameter, "[", strain_npriorvax_group_index, "]")) %>%
    select(parameter, strain, n_previous_vax, strain_npriorvax_group_index, true_value)

  # Generate shared individual effects (shared across years and subtypes)
  individual_effects_shared <- decensored_scaffold %>%
    select(individual) %>%
    unique() %>%
    mutate(u_peak_shared = rnorm(n(), mean = 0, sd = sigma_upeak_shared))
  
  # Generate subtype-specific random effects (shared across years)
  individual_effects_subtype <- decensored_scaffold %>%
    select(individual, subtype) %>%
    distinct() %>%
    mutate(u_peak_subtype = rnorm(n(), mean = 0, sd = sigma_upeak_subtype))

  # Generate year-specific random effects (shared across subtypes)
  individual_effects_year <- decensored_scaffold %>%
    select(individual, year) %>%
    distinct() %>%
    mutate(u_peak_year = rnorm(n(), mean = 0, sd = sigma_upeak_year)) %>%
    mutate(year = factor(year))
  
  # Assign parameters and simulate Hpeak
  simulated_data <- decensored_scaffold %>%
    distribute_parameter_values(
      universal_params = universal_params,
      strain_specific_params = strain_specific_params,
      RVE_params = RVE_params) %>%
    left_join(individual_effects_shared, by = "individual") %>%
    left_join(individual_effects_subtype, by = c("individual", "subtype")) %>%
    left_join(individual_effects_year, by = c("individual", "year")) %>%
    # Use sex-specific k_peak_cov2
    mutate(k_peak_cov2 = case_when(
      sex == "M" ~ k_peak_cov2_M,
      sex == "F" ~ k_peak_cov2_F
    )) %>%
    # Assign sex, age and bmi effects
    mutate(sex_effect = beta_M * (sex == "M"),
           age_effect = case_when(age_group == "31-40" ~ beta_age3140,
                                  age_group == "41-50" ~ beta_age4150,
                                  TRUE ~ 0),
           bmi_effect = case_when(BMI_group == "under" ~ beta_BMI_under,
                                  BMI_group == "over" ~ beta_BMI_over,
                                  TRUE ~ 0),
           smoking_effect = beta_smoking * smoking,
           asthma_effect = beta_asthma * asthma,
           time_effect = case_when(time_since_vax_group == "under28" ~ beta_time_under28,
                                   time_since_vax_group == "over42"  ~ beta_time_over42,
                                   TRUE ~ 0)
    ) %>%
    # Assign sigma specific to HAI/FRNT
    mutate(sigma = ifelse(str_detect(titer_type, "HAI"), sigma_HAI, sigma_FRNT)) %>%
    calculate_RVE() %>%
    calculate_Hpeak() %>%
    calculate_HLT() %>%
    calculate_Ht() %>%
    # Simulate an imperfect observation of the true post-vaccination titer
    mutate(Ht_latent = Ht,
           Ht = Ht_latent + rnorm(n(), mean = 0, sd = sigma)) %>%
    mutate(
      strain = factor(strain, levels = unique(strain_specific_params$strain)),
      subtype = factor(subtype, levels = subtype_levels)
    )

  if(include_prevax_measurement_error){
    # Simulate imperfect observation of the continuous pre-vaccination titer
    # that gave rise to the (imperfectly observed) post-vaccination titers
    imperfectly_observed_Hpre <- simulated_data %>%
      # Note that replicate measurements of the same pre-vaccination titer get independent measurement errors here
      select(individual, strain, year, Hpre, titer_type, measurement_replicate, sigma) %>%
      unique() %>%
      mutate(Hpre_latent = Hpre,
             Hpre = Hpre_latent + rnorm(n(), mean = 0, sd = sigma)) %>%
      select(-sigma)

    simulated_data <- simulated_data %>%
      select(-Hpre) %>%
      left_join(imperfectly_observed_Hpre, by = c("individual", "strain", "year", "measurement_replicate", "titer_type"))
  }

  plots <- plot_simulated_responses(simulated_data, params, include_NAI_infections)
  
  # Apply censoring if enabled (after plotting)
  if (censoring) {
    simulated_data <- simulated_data %>%
      mutate(
        Hpre_uncensored = Hpre,
        Ht_uncensored = Ht,
        Hpre = log(censor_titers(titer = 2^Hpre_uncensored, assay = titer_type), base = 2),
        Ht = log(censor_titers(titer = 2^Ht_uncensored, assay = titer_type), base = 2)
      )
  }

  # Additional plot of day 30 vs Hpre
  day30_vs_Hpre_plot <- plot_d30_vs_Hpre_synth_data(simulated_data)

  # If !include_titer_kinetics, output only day 30 data
  if(!include_kinetics){
    simulated_data <- simulated_data %>%
      filter(timepoint == 30)
  }

  return(list(data = simulated_data, peak_plot = plots$peak_plot,
              trajectory_plot = plots$trajectory_plot,
              day30_vs_Hpre_plot = day30_vs_Hpre_plot,
              Hpre_means_and_sds = Hpre_means_and_sds
            ))
}

# Get a unique integer index for each Hpre group
# (combination of individual, strain, year)
get_Hpre_group_index <- function(data){
  as.integer(factor(paste(data$individual, data$strain, data$year)))
}

get_strain_npriorvax_group_index <- function(data){
  as.integer(factor(paste(data$strain, data$n_previous_vax)))
}

# Given input data, finds empirical means/sds for pre-vaccination titers
# for each combination of strain / n. prior vaccinations
# Uses those values to empirically set prior means for the distribution of latent pre-vaccination titers
# (also defines SDs for the priors of both the mean and the sd of those distributions)
get_empirical_Hpre_means <- function(data){
  data %>%
    mutate(strain_npriorvax_group_index = get_strain_npriorvax_group_index(data)) %>%
    group_by(strain, n_previous_vax, strain_npriorvax_group_index) %>%
    # Set prior means for Hpre_mean and Hpre_sd using empirical values
    summarise(prior_mean_Hpre_mean = mean(Hpre),
              prior_mean_Hpre_sd = sd(Hpre),
              fraction_undetectable = mean(Hpre < log2(2 * undetectable_value))) %>%
    ungroup() %>%
    arrange(strain_npriorvax_group_index) %>%
    # For combinations with > 0.75 undetectable values, replace
    # the empirical SD with an arbitrary value of 0.5
    mutate(prior_mean_Hpre_sd = ifelse(fraction_undetectable > 0.75, 0.5, prior_mean_Hpre_sd)) %>%
    # We've set the prior means, now set the prior SDs for Hpre_mean and Hpre_sd
    mutate(
      prior_sd_Hpre_mean = 0.1,
      prior_sd_Hpre_sd = 0.1
    ) %>%
    select(strain_npriorvax_group_index, matches("prior_mean"), matches("prior_sd")) 
}


prepare_stan_input <- function(data, censoring, include_upeak_subtype, include_upeak_shared,
                               include_kinetics, include_RVE, include_upeak_year, include_k_peak_flu,
                               include_k_peak_cov2, RVE_model,
                               loglik_unit = 1, holdout_index = 0, include_Hpre_priors, include_NAI_infections,
                               include_time_since_vax, include_sex, include_age, include_BMI,
                               include_smoking, include_asthma){

  # Priors on the population distributions of latent pre-vaccination titers
  # (will be ignored if include_Hpre_priors = F)                              
  Hpre_priors <- get_empirical_Hpre_means(data)

  # Format prior infections column depending on whether or not we're including infections detected by NAI
  if(include_NAI_infections){
    infection_before_sample_subtype_matched = data$infection_before_sample_subtype_matched_PCR | data$infection_before_sample_subtype_matched_NAI
  }else{
    infection_before_sample_subtype_matched = data$infection_before_sample_subtype_matched_PCR
  }
  
  list(
    N_total = nrow(data),  # Total number of observations
    measurement_replicate = data$measurement_replicate, # Vector of measurement replicate number (1 or 2)
    Hpre = data$Hpre, # Vector of observed pre-vaccination titers (log2)
    Ht = data$Ht, # Vector of observed post-vaccination titers (log2)
    prior_mean_Hpre_mean =  Hpre_priors$prior_mean_Hpre_mean, # Priors for the mean and sd of pop. distribution of latent pre-vax titers
    prior_mean_Hpre_sd = Hpre_priors$prior_mean_Hpre_sd,
    prior_sd_Hpre_mean = Hpre_priors$prior_sd_Hpre_mean,
    prior_sd_Hpre_sd = Hpre_priors$prior_sd_Hpre_sd,
    t = data$t, # Vector with the number of days since vaccination in the current year.
    strain = as.integer(data$strain), # Vector with the vaccine strain for each observation, converted to an integer.
    K = length(unique(data$strain)), # Number of strains
    n_previous_vax = data$n_previous_vax, # Number vaccinations prior to year's vaccine
    I = length(unique(data$individual)),  # Total number of unique individuals
    individual = as.integer(factor(data$individual)),  # Convert individual IDs to integer indices
    N_ind_sub_combinations = n_distinct(data$individual, data$subtype),  # Number of unique individual-subtype combinations
    ind_subtype_combination = as.integer(factor(paste(data$individual, data$subtype))),  # Convert individual-subtype IDs to integer indices
    undetectable_value = data$undetectable_value, # Undetectable value (5 for rows that are HAI, 10 for FRNT) 
    max_dilution = as.integer(max(original_dilutions)), # Maximum dilution (non-log)
    infection_before_sample_subtype_matched = infection_before_sample_subtype_matched,
    recent_infection_before_sample_cov2 = data$recent_infection_before_sample_cov2_PCR, # SC2 infections ascertained by PCR. 
    # Number and indexing of individual/strain/year combinations (groups of post-vaccination titers with the same pre-vaccination titer)
    N_Hpre_groups = n_distinct(data$individual, data$strain, data$year),
    Hpre_group_index = get_Hpre_group_index(data),
    # Number and indexing of strain/n_priorvax combinations
    N_strain_npriorvax_groups = n_distinct(data$strain, data$n_previous_vax),
    strain_npriorvax_group_index = get_strain_npriorvax_group_index(data),
    # Individual-year combinations
    N_ind_year_combinations = n_distinct(data$individual, data$year),
    ind_year_combination = as.integer(factor(paste(data$individual, data$year))),
    # RVE levels for years 1:3 in wide format (to determine which repeat vaccination effect applies)
    # (values for years greater than the observation's year handled by a special level)
    RVE_level_Y1 = as.integer(data$RVE_level_Y1),
    RVE_level_Y2 = as.integer(data$RVE_level_Y2),
    RVE_level_Y3 = as.integer(data$RVE_level_Y3),
    # Wide vaccination status in years 1:3 for all observations
    vax_Y1 = as.integer(data$vax_Y1),
    vax_Y2 = as.integer(data$vax_Y2),
    vax_Y3 = as.integer(data$vax_Y3),
    N_RVE_effects = length(unique(RVE_model$form$RVE_level)),
    # Covariates
    sex = as.integer(factor(data$sex, levels = c("M", "F"))), # 1=male, 2=female
    age_group = as.integer(factor(data$age_group, levels = c("18-30", "31-40", "41-50"))),
    BMI_group = as.integer(factor(data$BMI_group, levels = c("healthy", "under", "over"))), # so 1 = healthy,2 = under,3 = over in Stan
    smoking = as.integer(data$smoking),
    asthma = as.integer(data$asthma),
    time_since_vax_group = as.integer(factor(data$time_since_vax_group, levels = c("28to42", "under28", "over42"))), # 1=reference(28-42), 2=under28, 3=over42
    titer_type = as.integer(factor(data$titer_type, levels = c("HAI", setdiff(unique(data$titer_type), "HAI")))), # 1=HAI, 2=other (FRNT)
    # Fit specifications
    include_upeak_subtype = as.integer(include_upeak_subtype),
    include_upeak_shared = as.integer(include_upeak_shared),
    include_kinetics = as.integer(include_kinetics),
    include_RVE = as.integer(include_RVE),
    include_upeak_year = as.integer(include_upeak_year),
    include_k_peak_flu = as.integer(include_k_peak_flu),
    include_k_peak_cov2 = as.integer(include_k_peak_cov2),
    censoring = as.integer(censoring),
    year = as.integer(as.character(data$year)), # Pass year as integer vector (as.character safely converts from factor)
    RVE_decay = as.integer(RVE_model$RVE_decay), # Pass RVE_decay switch
    loglik_unit = as.integer(loglik_unit), # 1 = factorize log-lik by individual/strain/year, 2 = factorize by individual (only matters for cross-validation but a value is required)
    holdout_index = as.integer(holdout_index), # Index of Hpre group (if loglik_unit is 1) or individual (if loglik_unit is 2) to be held-out. If 0, no unit is held out.
    include_Hpre_priors = as.integer(include_Hpre_priors),
    include_time_since_vax = as.integer(include_time_since_vax),
    include_sex = as.integer(include_sex),
    include_age = as.integer(include_age),
    include_BMI = as.integer(include_BMI),
    include_smoking = as.integer(include_smoking),
    include_asthma = as.integer(include_asthma)
  )
}

get_stan_diagnostics <- function(fit){
  
  diagnostics <- monitor(fit, print = F, warmup = 0) # warmup = 0 because warmup already excluded

  diagnostics <- as_tibble(diagnostics) %>%
    mutate(par = rownames(diagnostics)) %>%
    select(par, Rhat, valid, Bulk_ESS, Tail_ESS) %>%
    filter(!str_detect(par, "_transformed")) %>%
    filter(!str_detect(par, "LOO_LPD")) %>%  # exclude LOO LPD from parameter diagnostics
    mutate(par = case_when(
       par %in% c("sigma_upeak_shared[1]", "sigma_upeak_subtype[1]", "sigma_upeak_year[1]",
                  "f0[1]", "omega[1]", "tau[1]", "b_LT[1]", "r_LT[1]", "gamma[1]",
                  "beta_time_under28[1]", "beta_time_over42[1]",
                  "beta_M[1]", "beta_age3140[1]", "beta_age4150[1]",
                  "beta_BMI_under[1]", "beta_BMI_over[1]",
                  "beta_smoking[1]", "beta_asthma[1]",
                  "k_peak_flu[1]") ~ str_remove(par, "\\[1\\]"),
       T ~ par
    )) %>%
    mutate(par = case_when(
      par == "k_peak_cov2[1]" ~ "k_peak_cov2_M",
      par == "k_peak_cov2[2]" ~ "k_peak_cov2_F",
      T ~ par
    ))

  return(diagnostics)

}

# Function to extract and process posterior summaries for all parameters
process_stan_fit <- function(fit, data, RVE_model, chains = NULL, predicted_values_sample_size = 1000) {

  RVE_levels = get_RVE_model_levels(RVE_model)

  # = INDEX RETRIEVAL (linking categorical numerical indexing to variable levels) 

  # Data is used to determine the mapping of numbers to strains
  # Retrieve strain mapping
  strain_mapping <- data %>%
    select(subtype, strain) %>%
    mutate(strain_index = as.integer(strain)) %>%
    select(strain, strain_index, subtype) %>%
    unique()

  # Check that the number of strain mappings matches the number of strain-specific parameters
  stopifnot(length(unique(strain_mapping$strain_index)) == sum(str_detect(names(fit), "^b_peak")))
  stopifnot(length(unique(strain_mapping$strain_index)) == sum(str_detect(names(fit), "^a\\[")))

  # Retrieve individual mapping
  individual_mapping <- data %>%
    select(individual) %>%
    mutate(individual_index = as.integer(factor(individual))) %>%
    select(individual, individual_index) %>%
    unique()
  # Check that the number of individual mappings matches the number of individual effects
  stopifnot(length(unique(individual_mapping$individual_index)) == sum(str_detect(names(fit), "^u_peak_shared_transformed")))

  # Retrieve mapping of individual-by-subtype combinations
  ind_by_subtype_mapping <- data %>%
    select(individual, subtype) %>%
    mutate(ind_subtype_index = as.integer(factor(paste(individual, subtype)))) %>%
    select(individual, subtype, ind_subtype_index) %>%
    unique()
  
  # Check that the number of ind-by-subtype mappings matches the number of ind-by-subtype effects
  stopifnot(nrow(ind_by_subtype_mapping) == sum(str_detect(names(fit), "^u_peak_subtype_transformed")))
 
  # Retrieve mapping of individual-year combinations
  ind_year_mapping <- data %>%
    mutate(ind_year_index = as.integer(factor(paste(individual, year)))) %>%
    select(individual, year, ind_year_index) %>%
    unique()
  
  # Check that the number of ind-year mappings matches the number of ind-year effects
  stopifnot(nrow(ind_year_mapping) == sum(str_detect(names(fit), "^u_peak_year_transformed")))

  # = DEFAULT DIAGNOSTICS 
  # Note: diagnostic is done based on all chains, even if chains is not NULL.
  diagnostics <- get_stan_diagnostics(fit)

  # ===== EXTRACTING/PROCESSING POSTERIOR SAMPLE

  # Extract samples from all chains
  posterior_samples <- rstan::extract(fit, permuted = FALSE)

  # Keep only the specified chains
  if(is.null(chains)){
    # If chains is NULL, use all chains
    chains <- seq(1, dim(posterior_samples)[2])
  }

  # Tibble to store posterior samples from selected chains
  selected_samples <- tibble()

  for(chain in chains){
    samples_subset <- posterior_samples[, chain, ]
    samples_subset <- lapply(1:tail(dim(samples_subset), 1), function(i) samples_subset[, i])
    names(samples_subset) <- dimnames(posterior_samples)[[3]]

    param_names <- names(samples_subset)
    #param_names <- setdiff(names(samples_subset), "lp__")  # Exclude lp__


    # No need to process "_transformed" parameters
    # (they exist to handle cases where certain parameters are not estimated)
    param_names <- param_names[!str_detect(param_names, "_transformed")]
    
    # Convert posterior samples for this chain into a long-format tibble
    posterior_tibble <- bind_rows(
      lapply(param_names, function(param) {
        tibble(chain = chain, parameter = param, value = samples_subset[[param]]) %>%
          mutate(iteration = 1:n()) %>%
          select(chain, iteration, parameter, value)
        }
      )
    )
    
    selected_samples <- selected_samples %>%
      bind_rows(posterior_tibble)
  }

  # Remove the [1] from parameter names introduced by the on/off option in stan
  # (case_when produced memory errors on large fits, hence this approach)
  strip_index1 <- c("sigma_upeak_shared", "sigma_upeak_subtype", "sigma_upeak_year",
                    "f0", "omega", "tau", "b_LT", "r_LT", "gamma",
                    "beta_time_under28", "beta_time_over42",
                    "beta_M", "beta_age3140", "beta_age4150",
                    "beta_BMI_under", "beta_BMI_over",
                    "beta_smoking", "beta_asthma",
                    "k_peak_flu")
  for (p in strip_index1) {
    selected_samples$parameter[selected_samples$parameter == paste0(p, "[1]")] <- p
  }
  
  # Compute posterior summary statistics for the parameters
  posterior_summary_pars <-  selected_samples %>%
    filter(!str_detect(parameter, "^H_hat"), # Exclude predicted values...
           !str_detect(parameter, "^u_")) %>% # and individual effects
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    mutate(strain_index = ifelse(str_detect(parameter, "^(a|b_peak)\\["),
                   as.integer(str_extract(parameter, "(?<=\\[)[0-9]+")),
                   NA_integer_)) %>%
    mutate(RVE_index = ifelse(str_detect(parameter, "^r\\["),
                   as.integer(str_extract(parameter, "(?<=\\[)[0-9]+")),
                   NA_integer_)) %>%
    left_join(strain_mapping, by = c("strain_index")) %>%
    mutate(RVE_level = RVE_levels[RVE_index]) %>%
    select(-strain_index) %>%
    select(-RVE_index) %>%
    mutate(parameter = case_when(
      parameter == "k_peak_cov2[1]" ~ "k_peak_cov2_M",
      parameter == "k_peak_cov2[2]" ~ "k_peak_cov2_F",
      T ~ parameter
    ))
  
  #...and a posterior summary for each predicted value
  posterior_summary_pred <- selected_samples %>%
    filter(str_detect(parameter, "^H_hat")) %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975)
    ) %>%
    # This gets the predicted values in the right order
    mutate(observation_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
    arrange(observation_index) %>%
    select(-observation_index)

  # Return the raw parameter posterior samples in long format
  parameter_samples <- selected_samples %>%
    filter(!str_detect(parameter, "^H_hat"),
           !str_detect(parameter, "^u_")) %>% # exclude predicted values and random effects
    mutate(sample_index = paste0(chain, "-", iteration)) %>%
    select(-chain, -iteration) %>%
    left_join(
      posterior_summary_pars %>% filter(!is.na(strain)) %>% select(parameter, strain) %>% unique()
    ) %>%
    left_join(
      posterior_summary_pars %>% filter(!is.na(RVE_level)) %>% select(parameter, RVE_level) %>% unique(),
    )

  # Return individual effect samples in long format
  individual_effect_samples <- selected_samples %>%
    filter(str_detect(parameter, "^u_peak_shared")) %>%
    mutate(sample_index = paste0(chain, "-", iteration)) %>%
    select(-chain, -iteration) %>% 
    mutate(individual_index = str_extract(parameter, "[0-9]+") %>% as.integer()) %>%
    rename(u_peak_shared = value) %>%
    select(sample_index, individual_index, u_peak_shared) %>%
    left_join(individual_mapping, by = "individual_index")
    
  # Return individual-by-subtype effect samples in long format
  ind_by_subtype_effect_samples <- selected_samples %>%
    filter(str_detect(parameter, "^u_peak_subtype")) %>%
    mutate(sample_index = paste0(chain, "-", iteration)) %>%
    select(-chain, -iteration) %>%
    mutate(ind_subtype_index = str_extract(parameter, "[0-9]+") %>% as.integer()) %>%
    rename(u_peak_subtype = value) %>%
    select(sample_index, ind_subtype_index, u_peak_subtype) %>%
    left_join(ind_by_subtype_mapping, by = "ind_subtype_index")

  # Return individual-year effect samples in long format
  ind_year_effect_samples <- selected_samples %>%
    filter(str_detect(parameter, "^u_peak_year")) %>%
    mutate(sample_index = paste0(chain, "-", iteration)) %>%
    select(-chain, -iteration) %>%
    mutate(ind_year_index = str_extract(parameter, "[0-9]+") %>% as.integer()) %>%
    rename(u_peak_year = value) %>%
    select(sample_index, ind_year_index, u_peak_year) %>%
    left_join(ind_year_mapping, by = "ind_year_index")

  # Return the raw posterior sample of predicted titers, including full prediction (H_hat)
  # and fixed-effects (H_hat_fixed) prediction
  # (NOTE: For these, return keep only predicted_values_sample_size samples per chain)
  iterations_subsample <- sample(1:max(selected_samples$iteration), 
                               size = predicted_values_sample_size, 
                               replace = FALSE)

  predicted_value_samples <- selected_samples %>%
    filter(iteration %in% iterations_subsample) %>%
    filter(str_detect(parameter, "^H_hat")) %>%
    mutate(sample_index = paste0(chain, "-", iteration)) %>%
    mutate(observation_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
    mutate(parameter = str_remove(parameter, "\\[[0-9]+\\]")) %>%
    select(chain, iteration, sample_index, observation_index, parameter, value)
    
  predicted_value_samples <- predicted_value_samples %>%    
    pivot_wider(names_from = "parameter", values_from = "value")

  return(list(
    diagnostics = diagnostics,
    posterior_summary_pars = posterior_summary_pars,
    posterior_summary_pred = posterior_summary_pred,
    parameter_samples = parameter_samples,
    individual_effect_samples = individual_effect_samples,
    ind_by_subtype_effect_samples = ind_by_subtype_effect_samples,
    ind_year_effect_samples = ind_year_effect_samples,
    predicted_value_samples = predicted_value_samples,
    RVE_model = RVE_model
  ))
}

# Individual effects are exported as posterior samples from process_stan_fit
# This function summarizes those samples
summarize_individual_effect_samples <- function(model_results){

  # Summarizing u_peak_shared
  u_peak_shared_summary <- model_results$individual_effect_samples %>%
    group_by(individual) %>%
    summarise(
      u_peak_shared_mean = mean(u_peak_shared),
      u_peak_shared_lower = quantile(u_peak_shared, 0.025),
      u_peak_shared_upper = quantile(u_peak_shared, 0.975)
    )

  # Summarizing u_peak_subtype
  u_peak_subtype_summary <- model_results$ind_by_subtype_effect_samples %>%
    group_by(individual, subtype) %>%
    summarise(
      u_peak_subtype_mean = mean(u_peak_subtype),
      u_peak_subtype_lower = quantile(u_peak_subtype, 0.025),
      u_peak_subtype_upper = quantile(u_peak_subtype, 0.975)
    )

  # Summarizing u_peak_year
  u_peak_year_summary <- model_results$ind_year_effect_samples %>%
    group_by(individual, year) %>%
    summarise(
      u_peak_year_mean = mean(u_peak_year),
      u_peak_year_lower = quantile(u_peak_year, 0.025),
      u_peak_year_upper = quantile(u_peak_year, 0.975)
    )
  return(
    list(
      u_peak_shared = u_peak_shared_summary,
      u_peak_subtype = u_peak_subtype_summary,
      u_peak_year = u_peak_year_summary
    )
  )

}

plot_parameter_estimates <- function(model_results, include_Hpre_means_and_sd = T){

  input <- model_results$posterior_summary_pars %>%
            filter(!str_detect(parameter, "^u_"),  # Exclude individual effects
                   !str_detect(parameter, "Hpre_continuous"), # And continuous pre-vaccination titers
                   !str_detect(parameter, "LOO_LPD"),
                   !str_detect(parameter, "log_lik")) 
  if(!include_Hpre_means_and_sd){
    input <- input %>%
      filter(
        !str_detect(parameter, "Hpre_mean"),
        !str_detect(parameter, "Hpre_sd")
      )
  }

  input <- input %>%
    mutate(strain_label = strain_labeller(strain),
           #pair_label = vax_pair_labeller(consecutive_vax_pair),
           param_group = case_when(
                        # This groups strain-specific parameters into a single panel
                        !is.na(strain) ~ str_extract(parameter, "^[a-z|0-9|_]+"),
                        !is.na(RVE_level) ~ str_extract(parameter, "^[a-z|0-9|_]+"),
                        is.na(strain) & str_detect(parameter, "sigma") ~ "sigma",
                        is.na(strain) & str_detect(parameter, "Hpre_mean") ~ "Hpre_mean",
                        is.na(strain) & str_detect(parameter, "Hpre_sd") ~ "Hpre_sd",
                        is.na(strain) & str_detect(parameter, "k_peak") ~ "k_peak",
                        is.na(strain) & str_detect(parameter, "beta") ~ "covariates",
                        T ~ parameter)) %>%
                mutate(parameter_label = case_when(
                  !is.na(strain) ~ paste0("<", param_group, ">", strain_label),
                  #!is.na(consecutive_vax_pair) ~ pair_label,
                  str_detect(parameter, "Hpre_(sd|mean)\\[") ~ paste0("<", str_extract(parameter, "Hpre_(sd|mean)"), ">", str_extract(parameter, "(?<=\\[)[0-9]+(?=\\])")),
                  T ~ parameter
                ))

  parameter_order <- input %>%
    arrange(
      desc(str_detect(parameter_label, "H3")),
      desc(str_detect(parameter_label, "H1")),
      desc(str_detect(parameter_label, "B/Vic")),
      desc(str_detect(parameter_label, "B/Yam"))
    ) %>%
    pull(parameter_label)

  input <- input %>% mutate(parameter_label = factor(parameter_label, levels = parameter_order))

  pl <- ggplot(input,
                aes(x = parameter_label, y = mean)) +
                geom_point(color = "black", size = 3) +  # Posterior mean
                geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "black") +  # 95% CrI
                facet_wrap(~param_group, scales = "free") +  # Separate panels for parameter types, free y-axis
                theme_cowplot() +  # Use cowplot theme
                labs(x = "Parameter", y = "Estimate") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for clarity
                scale_x_discrete(labels = ~str_remove(.x, "<.+>"))  

  if("true_value" %in% names(model_results$posterior_summary_pars)){
    pl <- pl + 
      geom_point(aes(y = true_value), color = "red", size = 3)
  }

  return(pl)
}

combine_data_and_predictions <- function(data, model_results, censor_predictions){

  # This subsets posterior_summary_pred to keep just the full model predictions (excluding predictions for fixed-effects only)
  full_predictions <- model_results$posterior_summary_pred %>%
    filter(str_detect(parameter, "^H_hat\\[[0-9]+")) %>%
    mutate(observation_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
    arrange(observation_index) %>%
    select(-observation_index) 

  stopifnot(nrow(data) == nrow(full_predictions))

  combined_data <- data %>%
    bind_cols(full_predictions)

  if(censor_predictions){
    combined_data <- combined_data %>%
      mutate(across(all_of(c("mean", "lower", "upper")), ~log(censor_titers(titer = 2^.x, assay = titer_type), base = 2)))
  }

  combined_data <- combined_data %>%
    rename(H_hat_mean = mean,
           H_hat_upper = upper,
           H_hat_lower = lower)
  
  return(combined_data)

}

plot_obs_vs_predicted_titers <- function(data, model_results, censor_predictions){

  plot_data <- combine_data_and_predictions(data, model_results, censor_predictions)

  if(censor_predictions){
    ylab <- "Predicted titer (censored posterior mean)"
  }else{
    ylab <- "Predicted titer (posterior mean)"
  }
  
  plot_data %>%
    select(year, subtype, treatment, n_previous_vax, Ht, H_hat_mean, H_hat_upper, H_hat_lower) %>%
    mutate(n_previous_vax = factor(n_previous_vax),
           subtype = factor(subtype, levels = subtype_levels)) %>%
    ggplot(aes(x = Ht, y = H_hat_mean, color = n_previous_vax)) +
    geom_point(position = position_jitter(width = ifelse(censor_predictions, 0.1, 0),
                                          height = ifelse(censor_predictions, 0.1, 0))) +
    facet_grid(year ~ subtype) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_smooth(se = F, method = "lm") +
    theme_cowplot() +
    baseline_figure_settings +
    scale_color_manual(name = "Prior vaccinations", values = year_4_nprior_vax_colors) +
    scale_x_continuous(breaks = log(original_dilutions, base = 2),
                       labels = ~2^.x) +
    scale_y_continuous(breaks = log(original_dilutions, base = 2),
                       labels = ~2^.x,
                       sec.axis = sec_axis(~.x, name = 'Study year\n', breaks = NULL)) + 
    xlab("Observed titer") +
    ylab(ylab) +
    theme(legend.position = "top",
          axis.text.x = element_text(size = default_figure_font_size, angle = 45, hjust = 1))
    
}

plot_variance_partition <- function(data, model_results, by_individual_effect_type = F){


  # Proportion of observations that are HAI vs FRNT
  # (For the purposes of partitioning the variance, we'll define the
  #  standard deviation of residuals as the weighted average of sigma_HAI and sigma_FRNT)
  prop_HAI <- mean(data$titer_type == "HAI")
  prop_FRNT <- 1 - prop_HAI

  observed_titers <- data %>%
    mutate(observation_index = 1:n()) %>%
    select(observation_index, Ht)

  # Start with the posterior sample of predicted titers (full prediction and fixed-effects predictions)
  variance_partition <- model_results$predicted_value_samples %>%
    # Left join the observed titers
    left_join(observed_titers, by = "observation_index") %>%
    # For each posterior sample, estimate sigma2_u (individual effects variance) and sigma2_fixed
    # (variance explained by fixed effects predictions)
    group_by(sample_index) %>%
    summarise(
      sigma2_fixed = var(H_hat_fixed),
      # Below is a post-hoc estimate of the variance of individual effects using deviations 
      # between H_hat and H_hat_fixed. We use it when the model is non-linear (includes kinetics)
      # (in which case the individual effects don't modify the observations directly,
      # so we can't estimate the total attributable to individual effects by summing their
      # variance parameters)
      # If model didn't assume kinetics, this column will be ignored below
      sigma2_u_posthoc = sum((H_hat - H_hat_fixed)^2) / (n() - 1))

  # For the same posterior samples, left join the values of standard deviation parameters
  # (residual sigma, sigma for individual effects)
  variance_partition <- variance_partition %>%
    left_join(model_results$parameter_samples %>%
                filter(str_detect(parameter, "sigma")) %>%
                pivot_wider(names_from = parameter, values_from = value) %>%
                # Define sigma (residual variance standard deviation) as weighted sum across FRNT and HAI
                mutate(sigma = sigma_HAI * prop_HAI + sigma_FRNT * prop_FRNT) %>%
                select(-sigma_HAI, -sigma_FRNT) %>%
                mutate(across(matches("sigma"), ~ .x^2, .names = "{.col}2")) %>%
                select(sample_index, sigma2, matches("^sigma_.*2$")),
              by = "sample_index")

  # If the model included kinetics...
  if(any(kinetic_par_names %in% model_results$posterior_summary_pars$parameter)){
    # Use the post-hoc estimate of individual effects variance
    variance_partition <- variance_partition %>%
      rename(sigma2_u = sigma2_u_posthoc) %>%
      select(-matches('sigma_upeak'))
  }else{
    # If model didn't include kinetics, can estimate sigma u as sum of sigma upeaks
    variance_partition <- variance_partition %>%
      mutate(sigma2_u = rowSums(select(., matches('sigma_upeak.*2')), na.rm = T)) %>%
      select(-sigma2_u_posthoc)
  }

  fraction_variance <- variance_partition %>%
    mutate(total_variance = sigma2_fixed + sigma2_u + sigma2) %>%
    pivot_longer(cols = matches("sigma"),
                  names_to = "variance_type",
                  values_to = "variance") %>%
    group_by(sample_index) %>%
    mutate(fraction = variance / total_variance) %>%
    mutate(variance_type = factor(variance_type))
                          
  if(by_individual_effect_type){
    input_tibble <- fraction_variance %>%
      filter(variance_type != "sigma2_u")

    partition_levels <- c("sigma2_fixed",
                          input_tibble %>% 
                            filter(str_detect(variance_type, "upeak")) %>%
                            pull(variance_type) %>%
                            unique() %>%
                            as.character(),
                          "sigma2"
                          )    
  }else{
    input_tibble <- fraction_variance %>%
      filter(!str_detect(variance_type, "upeak"))
    
    partition_levels <- c("sigma2_fixed", "sigma2_u", "sigma2")
  }
  input_tibble <- input_tibble %>%
    mutate(variance_type = factor(variance_type, levels = partition_levels))
      
  CrIs <- input_tibble %>%
    group_by(variance_type) %>%
    summarise(mean = mean(fraction),
              lower = quantile(fraction, 0.025),
              upper = quantile(fraction, 0.975))

  pl <- input_tibble %>%
    ggplot(aes(x = fraction)) +
    geom_histogram() +
    facet_wrap(~variance_type, scales = "free", nrow = 1) +
    theme_cowplot() +
    geom_vline(data = CrIs, aes(xintercept = mean), color = "red") +
    geom_vline(data = CrIs, aes(xintercept = lower), color = "red", linetype = 2) +
    geom_vline(data = CrIs, aes(xintercept = upper), color = "red", linetype = 2) +
    geom_label(data = CrIs, aes(x = mean, y = 0, label = paste0(round(mean, 2), " (", round(lower, 2), ", ", round(upper, 2), ")")),
              vjust = -1, hjust = 0.5) +
    xlab("Fraction") +
    ylab("Number of posterior samples")

  return(pl)

}

compute_average_covariates <- function(data, include_NAI_infections) {

  input_data <- data

  if(include_NAI_infections){
    input_data <- input_data %>%
      mutate(infection_before_sample_subtype_matched = infection_before_sample_subtype_matched_PCR | infection_before_sample_subtype_matched_NAI)
  }else{
    input_data <- input_data %>%
      mutate(infection_before_sample_subtype_matched = infection_before_sample_subtype_matched_PCR)
  }

  input_data <- input_data %>%
    mutate(recent_infection_before_sample_cov2 = recent_infection_before_sample_cov2_PCR)

  input_data %>%
    select(year, subtype, infection_before_sample_subtype_matched, recent_infection_before_sample_cov2, sex,
           age_group, smoking, asthma, BMI_group, time_since_vax_group) %>%
    unique() %>%
    group_by(year, subtype) %>%
    summarise(
      infection_before_sample_subtype_matched = mean(infection_before_sample_subtype_matched),
      recent_infection_before_sample_cov2 = mean(recent_infection_before_sample_cov2),
      fraction_male = sum(sex == "M") / n(),
      fraction_female = 1 - fraction_male,
      fraction_age3140 = sum(age_group == "31-40") / n(),
      fraction_age4150 = sum(age_group == "41-50") / n(),
      fraction_smoking = sum(smoking) / n(),
      fraction_asthma = sum(asthma) / n(),
      fraction_BMI_under = sum(BMI_group == "under") / n(),
      fraction_BMI_over = sum(BMI_group == "over") / n(),
      fraction_time_under28 = sum(time_since_vax_group == "under28") / n(),
      fraction_time_over42 = sum(time_since_vax_group == "over42") / n()
    ) %>%
    ungroup()
}

get_Ht_vs_Hpre_fixed_effects <- function(params_summary, params_sample, param_subsample_size,
                                         data, t, use_oficial_timepoint, data_time_tolerance,
                                         remove_covariate_effects, include_NAI_infections){
                                                                                
  # Get predicted relationship between Ht and Hpre for people
  # with different numbers of previous vaccinations (for plotting with plot_Ht_vs_Hpre_predictions)

  # t can be an integer representing arbitrary number of days after vaccination,
  # or a character string equal to either "Hpeak" or "HLT"

  # If use_oficial_timepoint is TRUE, we'll take all data points from the official timepoint 
  # (e.g., we take all "day 30" samples regardless of how many days after vaccination they were taken exactly)
  if(use_oficial_timepoint){
    stopifnot(t %in% c(30, 182))
    if(!is.null(data_time_tolerance)){
      stop("data_time_tolerance must be NULL if use_oficial_timepoint is TRUE")
    }
  }


  # We ignore individual effects. For all other fixed-effects acting on Hpeak, we
  # compute the average value of the predictor in the population, then apply the corresponding effect
  

  average_predictors <- compute_average_covariates(data, include_NAI_infections)

  if(remove_covariate_effects){
    average_predictors <- average_predictors %>%
      mutate(across(!any_of(c("year", "subtype")), ~ 0))
  }

  prevax_titer_range <- seq(
    min(log2(original_dilutions)),
    max(log2(original_dilutions)),
    length.out = 100
  )
  
  scaffold <- tibble(Hpre = prevax_titer_range) %>%
    expand_grid(data %>%
                  select(year, strain, subtype, n_previous_vax, matches("RVE_level"), matches("vax_Y")) %>%
                  unique()) %>%
    arrange(year, subtype, strain) %>%
    left_join(average_predictors, by = c("year", "subtype"))

  # Create tibble for plotting the posterior average of the Ht vs Hpre lines
  input_tibble_posterior_average <- scaffold %>%
    distribute_parameter_values(universal_params = params_summary$universal,
                                strain_specific_params = params_summary$strain_specific,
                                RVE_params = params_summary$RVE) %>%
    # Give it a "sample_index" of "posterior_average", in case we're also plotting 
    # uncertainty using a posterior sample
    mutate(sample_index = "posterior_average") %>%
    select(sample_index, everything())

  # If a sample of parameter values from the posterior was provided (and not just the posterior summary)
  if(!is.null(params_sample)){

    # Downsample posterior samples to param_subsample_size if necessary
    if(length(unique(params_sample[[1]]$sample_index)) > param_subsample_size){
      posterior_subsample_indices <- sample(unique(params_sample[[1]]$sample_index),
                                  size = param_subsample_size,
                                  replace = F)
    }else{
      posterior_subsample_indices <- unique(params_sample[[1]]$sample_index)
    }

    params_subsample <- lapply(params_sample, function(df) {
      if(!is.null(df)){
        # For each of the objects in params_subsample, subset to index subsample
        df %>% filter(sample_index %in% posterior_subsample_indices)
      }
    })

    # Create an input tibble for plotting a sample of predicted Ht vs Hpre lines
    input_tibble_posterior_sample <- bind_rows(
      # For each sampled index
      lapply(posterior_subsample_indices, function(idx) {
        # Take the parameter values for that sample...
        params_slice <- lapply(params_subsample, function(df){
          if(!is.null(df)){
            df %>% filter(sample_index == idx)
          }}) 
        # Distribute them in the scaffold
        scaffold %>%
          distribute_parameter_values(
            universal_params = params_slice$universal,
            strain_specific_params = params_slice$strain_specific,
            RVE_params = params_slice$RVE) %>%
          mutate(sample_index = idx)
        })
      ) %>%
      select(sample_index, everything())
  }else{
    input_tibble_posterior_sample <- c()
  }

  input_tibble <- bind_rows(
    input_tibble_posterior_average,
    input_tibble_posterior_sample
   ) %>%
   # setting individual effects to 0 
    mutate(u_peak_shared = 0, u_peak_subtype = 0, u_peak_year = 0)
    
  # If kinetics are not included, set tau to 0 (instant peak), f0 to 1, omega, r_LT and b_LT to 0 (no waning)
  if(!all(c("f0", "omega", "b_LT", "r_LT", "tau") %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(tau = 0, f0 = 1, omega = 0, r_LT = 0, b_LT = 0)
  }

  # If k_peak_flu not included, set to 0
  if(!("k_peak_flu" %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(k_peak_flu = 0)
  }

  # if beta_M not included, set to 0
  if(!("beta_M" %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(beta_M = 0)
  }    

  # If beta_age parameters not included, set to 0
  if(!all(c("beta_age3140", "beta_age4150") %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(beta_age3140 = 0, beta_age4150 = 0)
  }

  # If beta_smoking not included, set to 0
  if(!("beta_smoking" %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(beta_smoking = 0)
  }

  # If beta_asthma not included, set to 0
  if(!("beta_asthma" %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(beta_asthma = 0)
  }

  # If beta_BMI_under/over not included, set to 0
  if(!all(c("beta_BMI_under", "beta_BMI_over") %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(beta_BMI_under = 0, beta_BMI_over = 0)
  }

  # If k_peak_cov_2_M/F parameters not included, set them to 0
  if(!all(c("k_peak_cov2_M", "k_peak_cov2_F") %in% names(params_summary))){
    input_tibble <- input_tibble %>%
      mutate(k_peak_cov2_M = 0, k_peak_cov2_F = 0)
  }

  # If RVE decay not included, set to zero
  if(!("gamma" %in% names(params_summary$universal))){
    input_tibble <- input_tibble %>%
      mutate(gamma = 0)
  }

  # If time-since-vax parameters not included, set to zero
  if(!all(c("beta_time_under28", "beta_time_over42") %in% names(params_summary$universal))){
    input_tibble <- input_tibble %>%
      mutate(beta_time_under28 = 0, beta_time_over42 = 0)
  }
 
  # Compute the average effects of covariates
  input_tibble <- input_tibble %>%
    mutate(
      k_peak_cov2 = k_peak_cov2_M * fraction_male + k_peak_cov2_F * fraction_female,
      sex_effect = beta_M * fraction_male,
      age_effect = beta_age3140 * fraction_age3140 + beta_age4150 * fraction_age4150,
      smoking_effect = beta_smoking * fraction_smoking,
      asthma_effect = beta_asthma * fraction_asthma,
      bmi_effect = beta_BMI_under * fraction_BMI_under + beta_BMI_over * fraction_BMI_over,
      time_effect = beta_time_under28 * fraction_time_under28 + beta_time_over42 * fraction_time_over42
    )

  input_tibble <- input_tibble %>%
    calculate_RVE() %>%
    calculate_Hpeak() %>%
    calculate_HLT()


    if(t == "Hpeak"){
      input_tibble <- input_tibble %>%
        mutate(predicted_titer = Hpeak)
    }else{
      if(t == "HLT"){
        input_tibble <- input_tibble %>%
          mutate(predicted_titer = HLT)
    }else{
      stopifnot(is.numeric(t))
      input_tibble <- input_tibble %>%
        mutate(t = ifelse(t == 0, 1e-6, t)) %>% # Assume we're observing everyone this many days after vaccination
        calculate_Ht() %>%
        mutate(predicted_titer = Ht)
    }
  }


  # Now that predictions have been computed, if we're plotting uncertainty, get the 95% range at each pre-vaccination titer
  # (This will be empty if we're not plotting uncertainty, but that's okay)
  prediction_CrI <- input_tibble %>%
    filter(sample_index != "posterior_average") %>%
    group_by(year, subtype, strain, n_previous_vax, Hpre) %>%
    summarise(predicted_titer_upper = quantile(predicted_titer, 0.975),
              predicted_titer_lower = quantile(predicted_titer, 0.025))

  return(list(posterior_average = input_tibble %>% filter(sample_index == "posterior_average"),
              prediction_CrI = prediction_CrI))

}


plot_Ht_vs_Hpre_predictions <- function(params_summary, params_sample = NULL, param_subsample_size = 500,
                                        data, t, use_oficial_timepoint = F, plot_data = F, data_time_tolerance = 3,
                                        include_NAI_infections, remove_covariate_effects = F, use_latent_prevax_titers = F){
  # Plot predicted relationship between Ht and Hpre as separate lines for people
  # with different numbers of previous vaccinations

  # t can be an integer representing arbitrary number of days after vaccination,
  # or a character string equal to either "Hpeak" or "HLT"

  # We ignore individual effects. For all other fixed-effects acting on Hpeak, we
  # compute the average value of the predictor in the population, then apply the corresponding effect 
  # (see get_Ht_vs_Hpre_fixed_effects)

  fixed_effects_predictions <- get_Ht_vs_Hpre_fixed_effects(
    params_summary = params_summary,
    params_sample = params_sample,
    param_subsample_size = param_subsample_size,
    data = data,
    t = t,
    use_oficial_timepoint = use_oficial_timepoint,
    data_time_tolerance = data_time_tolerance,
    remove_covariate_effects = remove_covariate_effects,
    include_NAI_infections = include_NAI_infections)                                         

  x_axis_empirical_range <- data %>%
    group_by(strain, year, n_previous_vax) %>%
    summarise(Hpre_empirical_min = min(Hpre), Hpre_empirical_max = max(Hpre)) %>%
    ungroup()

  # Subset to plot only predictions up to 2 * max dilution
  fixed_effects_predictions$posterior_average <- fixed_effects_predictions$posterior_average %>%
    filter(2^predicted_titer <= 2 * max(original_dilutions)) %>%
    # Subset x-axis range to range of empirical range, + 1 on either side
    left_join(x_axis_empirical_range) %>%
    filter(Hpre >= Hpre_empirical_min - 1, Hpre <= Hpre_empirical_max + 1) %>%
    select(-Hpre_empirical_min, -Hpre_empirical_max)

  fixed_effects_predictions$prediction_CrI <- fixed_effects_predictions$prediction_CrI %>%
    semi_join(fixed_effects_predictions$posterior_average)

  input_tibble <- fixed_effects_predictions$posterior_average %>%
    mutate(solid_line = T) %>%
    mutate(subtype = factor(subtype, levels = subtype_levels))

  fixed_effects_predictions <- fixed_effects_predictions$prediction_CrI %>%
    mutate(subtype = factor(subtype, levels = subtype_levels))

  if(t == "Hpeak"){
    y_axis_label = "Peak titer"
    # If plotting Hpeak, will make lines dashed if r effect not significant
    geom_line_mapping = aes(linetype = solid_line)
    line_guides_spec = 
      guides(linetype = guide_legend(override.aes = list(linewidth = 0.7, size = 2)))
    if("upper_r" %in% names(input_tibble)){
      scale_ltype_spec = 
        scale_linetype_manual(name = "Significant RVE\non intercept", values = c(2, 1))
    }else{
      scale_ltype_spec = 
        scale_linetype_manual(name = "Significant RVE\non intercept", values = c(1))
    }
  }else{
    scale_ltype_spec = NULL
    geom_line_mapping = NULL
    line_guides_spec = NULL
    if(t == "HLT"){
      y_axis_label = "Long-term stable titer"
    }else{
      y_axis_label = paste0("Observed titer ", t, " days after vaccination")
    }
  }

  if(use_latent_prevax_titers){
    stopifnot(plot_data = T)
    x_axis_label <- "Latent pre-vaccination titer (posterior average)"
    x_axis_breaks <- c(0.01, 0.1, 1, original_dilutions[c(2, 4, 6, 8)])
  }else{
    x_axis_label <- "Pre-vaccination titer"
    x_axis_breaks <- original_dilutions[c(2, 4, 6, 8)]
  }
  y_axis_breaks <- original_dilutions[c(2, 4, 6, 8)]

  pl <- input_tibble %>%
    filter(sample_index == "posterior_average") %>%
    ggplot(aes(x = Hpre, y = predicted_titer, color = factor(n_previous_vax))) +
    geom_ribbon(
      data = fixed_effects_predictions,
      aes(x = Hpre, ymin = predicted_titer_lower,
          ymax = predicted_titer_upper, fill = factor(n_previous_vax)),
      alpha = 0.3,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_line(linewidth = 0.7, mapping = geom_line_mapping) +
    facet_grid(year ~ subtype, labeller = labeller(year = function(x) paste("Year", x))) +
    theme_cowplot() +
    baseline_figure_settings + 
    labs(x = x_axis_label, y = y_axis_label, color = "Prior vaccinations") +
    scale_x_continuous(breaks = log(x_axis_breaks, base = 2),
                       labels = ~2^.x) +
    scale_y_continuous(breaks = log(y_axis_breaks, base = 2),
                       labels = ~2^.x) +
    baseline_figure_settings + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') +
    scale_ltype_spec +
    line_guides_spec + 
    theme(legend.key.width = unit(0.5, "cm"), plot.margin = margin(t = 1, l = 8),
         legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) +
    scale_color_manual(values = year_4_nprior_vax_colors) +
    scale_fill_manual(values = year_4_nprior_vax_colors)

  if(plot_data){

    data_points <- data %>%
      # In case we'll use latent pre-vaccination titers.
      mutate(Hpre_group_index = get_Hpre_group_index(data))

    if(use_oficial_timepoint){
      data_points <- data_points %>%
        filter(timepoint == !!t) 
      title <- ""
    }else{
      # Plot data points within data_time_tolerance days of t.
      data_points <- data_points %>%
        filter(abs(t - !!t) <= data_time_tolerance) %>%
      title <- paste0("Note: observed data points are within ", data_time_tolerance, " days of day ", t)
    }

    data_points <- data_points %>%
        select(year, strain, subtype, n_previous_vax, matches("Hpre_group"), Hpre, Ht)

    if(use_latent_prevax_titers){
      # If this is true, we'll plot the estimated latent, rather than the observed
      # pre-vaccination titer for each observation.
      data_points <- data_points %>%
        left_join(params_summary$latent_prevax_titers %>%
                   rename(latent_prevax_titer = mean) %>%
                   select(Hpre_group_index, latent_prevax_titer)) %>%
        mutate(Hpre = latent_prevax_titer)  
    }

    if(nrow(data_points) > 0){
     pl <- pl +
        geom_point(data = data_points %>%
                            mutate(subtype = factor(subtype, levels = subtype_levels)),
                   aes(x = Hpre, y = Ht),
                  position = position_jitter(width = 0.1, height = 0.1),
                  alpha = 1, size = main_text_point_size,
                  shape = boxplot_point_shape, stroke = boxplot_point_stroke) +
                  ggtitle(title)

     if(use_latent_prevax_titers){
      pl <- pl + geom_vline(xintercept = 0, linetype = 2)
     }
    }else{
      print(paste("No data points within", data_time_tolerance, "days of", t))
    }
  }

  return(pl)
}

# Takes summary of parameter estimates produced by process_stan_fit and formats it for plotting
# (i.e., makes it so that all strains have a value for each parameter, including parameters shared across strains)
# (similar to strain_params for generating synthetic data)
format_params_for_plotting <- function(model_results, posterior_summary = T){

  # If posterior_summary is TRUE, will take a single set of values corresponding to the posterior means and CrI bounds,
  # then format it for plotting.
  # If posterior_summary is FALSE, will take each sample from the posterior sample and format the parameters for plotting

  if(posterior_summary == T){
    parameter_estimates <- model_results$posterior_summary_pars
    values_from_cols = c("mean", "lower", "upper")
  }else{
    parameter_estimates <- model_results$parameter_samples
    values_from_cols = c("value")
  }

  strain_specific_params <- parameter_estimates %>%
    select(-matches("true_value")) %>% # Handles fits to synthetic data
    filter(!is.na(strain)) %>%
    mutate(parameter = str_remove(parameter, "\\[[0-9]+\\]")) %>%
    pivot_wider(names_from = parameter, values_from = any_of(values_from_cols)) %>%
    select(-matches("RVE_level"), -matches("subtype")) %>%
    rename_with(~str_remove(.x, "mean_"), starts_with("mean_"))

  RVE_params <- parameter_estimates %>%
    select(-matches("true_value")) %>% # Handles fits to synthetic data
    filter(!is.na(RVE_level)) %>%
    mutate(parameter = str_remove(parameter, "\\[[0-9]+\\]")) %>%
    pivot_wider(names_from = parameter, values_from = any_of(values_from_cols)) %>%
    select(-strain, -matches('subtype')) %>%
    rename_with(~str_remove(.x, "mean_"), starts_with("mean_"))

  universal_params <- parameter_estimates %>%
    select(-matches("true_value")) %>% # Handles fits to synthetic data
    filter(!str_detect(parameter, "^u_"), parameter != "lp__") %>%  # Exclude individual effects and the log-posterior
    # Because we plot predictions for an arbitrary range of pre-vaccination titers, no need to keep
    # Hpre parameters
    filter(!str_detect(parameter, "Hpre_"),
           !str_detect(parameter, "log_lik")) %>%
    filter(is.na(strain), is.na(RVE_level)) %>%
    select(parameter, matches("sample_index"), any_of(values_from_cols)) %>%
    pivot_wider(names_from = parameter, values_from = any_of(values_from_cols)) %>%
    rename_with(~str_remove(.x, "mean_"), starts_with("mean_"))

  if(posterior_summary == T){
    latent_prevax_titers <- parameter_estimates %>%
      select(-matches("true_value"), -matches("strain"), -matches("RVE_level")) %>%
      filter(str_detect(parameter, "Hpre_continuous")) %>%
      mutate(Hpre_group_index = as.integer(str_extract(parameter, "[0-9]+")))

  }else{
    # No need to export this when exporting a raw sample from the posterior.
    latent_prevax_titers <- NULL
  }

  # If RVE parameters were not estimated (model with no RVEs), set them to 0
  if(nrow(RVE_params) == 0){
    RVE_params <- tibble(
      # No-RVE models is built on the scaffold of one of the RVE models
      RVE_level = get_RVE_model_levels(model_results$RVE_model)
      ) %>%
      filter(RVE_level != "NA") %>% # Ignore special level passed to stan coded as "NA"
      mutate(r = 0, lower_r = 0, upper_r = 0)

    if(!posterior_summary){
      RVE_params <- RVE_params %>%
        cross_join(universal_params %>%
                    select(sample_index)) %>%
        select(sample_index, RVE_level, r)
    }
  }

  formatted_params <- list(
    universal = universal_params,
    strain_specific = strain_specific_params,
    RVE = RVE_params,
    latent_prevax_titers = latent_prevax_titers
  )

  return(formatted_params)
}

# Scatterplot of posterior parameter samples (for non-strain-specific parameters)
plot_parameter_correlations <- function(model_results){
  wide_format_data <- model_results$parameter_samples %>%
    # Exclude predicted values, individual effects, 
    filter(!str_detect(parameter, "^H_hat"),
           !str_detect(parameter, "^u_")) %>%
    select(-strain, -subtype) %>%
    pivot_wider(names_from = parameter, values_from = value) %>% 
    select(-sample_index)

  pl <- wide_format_data %>%
    ggpairs(columns = names(wide_format_data)[!str_detect(names(wide_format_data), "\\[")])

  return(pl)
}

plot_model_predictions_by_timepoint <- function(data, model_results, year, censor_predictions, transpose = F, subtype = NULL){

  predicted_value_posterior_averages <- model_results$predicted_value_samples %>%
    group_by(observation_index) %>%
    summarise(predicted_mean = mean(H_hat),
              lower = quantile(H_hat, 0.025),
              upper = quantile(H_hat, 0.975)) %>%
    ungroup()

  stopifnot(nrow(predicted_value_posterior_averages) == nrow(data))

  input_data <- data %>%
    rename(log2_titer = Ht) %>%
    mutate(observation_index = 1:n()) %>%
    left_join(predicted_value_posterior_averages,
              by = "observation_index") %>%
    mutate(subtype = factor(subtype, levels = subtype_levels))

  if(censor_predictions){
    input_data <- input_data %>%
      mutate(across(all_of(c("predicted_mean", "lower", "upper")), ~log(censor_titers(titer = 2^.x, assay = titer_type), base = 2)))
  }

  if(!is.null(subtype)){
    input_data <- input_data %>% filter(subtype %in% !!subtype)
  }

  pl <- plot_postvax_response(
          response_data = input_data, response_var = "log2_titer", 
          year = year,
          prior_infection_timing = "before_year_vdate",
          prior_infection_virus = "subtype_matched",
          vaccine_strains_only = T,
          show_points = F,
          show_pairwise_comparisons = F,
          transpose = transpose) +
    geom_boxplot(aes(y = predicted_mean, color = "Model prediction"),
                 alpha = 0.3,
                 # Tweak x to shift the boxplot of model predictions left or right
                 position = position_nudge(x = 0.15, y = 0), outlier.alpha = 0,
                 box.linewidth = boxplot_line_width * 0.4, median.linewidth =  boxplot_line_width * 1.5) +
    scale_color_manual(name = NULL, values = c("Model prediction" = "#fe7b3f")) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme(legend.position = "top")

  return(pl)
              
}

plot_d30_vs_Hpre_synth_data <- function(simulated_data) {
  # Helper function to create each plot
  make_plot <- function(data, x_var, y_var) {
    ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = factor(n_previous_vax))) +
      geom_point() +
      facet_grid(year ~ subtype) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      theme(legend.position = "top", axis.text.x = element_text(angle = -90, vjust = 0.5)) +
      scale_x_continuous(breaks = log2(original_dilutions), labels = function(x){2^x}) +
      scale_y_continuous(breaks = log2(original_dilutions), labels = function(x){2^x}) +
      geom_hline(aes(yintercept = log2(undetectable_value)), color = 'gray80') +
      geom_hline(aes(yintercept = log2(1280)), color = "gray80") +
      geom_vline(aes(xintercept = log2(undetectable_value)), color = "gray80")
  }

  d <- simulated_data %>% filter(timepoint == 30)
  p1 <- make_plot(d, "Hpre_latent", "Ht_latent")

  if (all(c("Hpre_uncensored", "Ht_uncensored") %in% names(d))) {
    p2 <- make_plot(d, "Hpre_uncensored", "Ht_uncensored")
    p3 <- make_plot(d, "Hpre", "Ht")
    pl_grid <- plot_grid(p1, p2, p3, nrow = 1)
  } else {
    p2 <- make_plot(d, "Hpre", "Ht")
    pl_grid <- plot_grid(p1, p2, nrow = 1)
  }
  return(pl_grid)
}

plot_RVE_estimates <- function(model_results){

  RVE_levels <- get_RVE_model_levels(model_results$RVE_model)

  params <- format_params_for_plotting(model_results, posterior_summary = T)
  RVE_params <- params$RVE

  # If RVE_model is the fully flexible model...
  if(identical(RVE_levels, get_RVE_model_levels(fully_flexible_RVE_model))){

    # Plot separate lines for each observation year/subtype
    input <- model_results$RVE_model$form %>%
      select(-matches("default_synthetic_value")) %>%
      left_join(RVE_params) %>%
      mutate(time_difference = year - previous_year) %>%
      mutate(plot_group = paste(subtype, year, sep = "-")) %>%
      mutate(subtype = factor(subtype, levels = rev(subtype_levels)))

    input <- input %>%
      left_join(
        input %>%
          select(plot_group) %>%
          unique() %>%
          mutate(manual_jitter = rnorm(n(), sd = 0.1))
      )


    input %>%
      ggplot(aes(x = time_difference + manual_jitter,
                 y = r, color = factor(year),
                 group = factor(year))) +
      geom_point() +
      geom_linerange(aes(ymin = lower_r, ymax = upper_r)) +
      geom_line() +
      facet_wrap('subtype', nrow = 1) +
      geom_text(aes(label = pair_label, color = NULL),
                 position = position_nudge(x = 0.6, y = 0),
                 show.legend = F,
                 size = 5) +
      geom_hline(yintercept = 0, linetype = 2) +
      xlab("Prior vaccination\n(most recent to oldest)") +
      scale_color_discrete(name = "Observation year") +
      theme(legend.position = 'top') +
      scale_x_continuous(breaks = unique(input$time_difference), limits = c(0.5,4))
  
  
  
  }else{
    if(identical(get_RVE_model_levels(model_results$RVE_model), get_RVE_model_levels(fully_shared_model)) ||
       identical(get_RVE_model_levels(model_results$RVE_model), get_RVE_model_levels(update_specific_model)) ||
       identical(get_RVE_model_levels(model_results$RVE_model), get_RVE_model_levels(subtype_specific_model))){

      # If model is the fully shared RVE model, plot the decay of the shared baseline r over time
      # with uncertainty
      scaffold <- tibble(years_since_prior_vaccine = 1:3)

      r_sample <- model_results$parameter_samples %>%
        filter(str_detect(parameter, "r\\[")) %>%
        mutate(parameter = "r") %>%
        select(parameter, RVE_level, value, sample_index) %>%
        pivot_wider(names_from = parameter, values_from = value)

      gamma_sample <- model_results$parameter_samples %>%
        filter(parameter == "gamma") %>%
        select(value, sample_index) %>%
        rename(gamma = value)


      RVE_trajectory <- r_sample %>%
        left_join(gamma_sample, by = join_by(sample_index)) %>%
        expand_grid(scaffold) %>%
        arrange(years_since_prior_vaccine) %>%
        mutate(R = r * (1 + gamma) ^ (years_since_prior_vaccine - 1)) %>%
        group_by(years_since_prior_vaccine, RVE_level) %>%
        summarise(mean_R = mean(R),
                  lower_R = quantile(R, 0.025),
                  upper_R = quantile(R, 0.975))
      
      # If model is subtype-specific RVE model, order subtypes using default ordering.
      if(identical(get_RVE_model_levels(model_results$RVE_model), get_RVE_model_levels(subtype_specific_model))){
        RVE_trajectory <- RVE_trajectory %>%
          mutate(RVE_level = factor(RVE_level, levels = subtype_levels))
      }else{
        RVE_trajectory <- RVE_trajectory %>%
          mutate(RVE_level = str_replace(RVE_level, "_", " ")) 
      }

      # If model is the update-specific model, use hand-picked colors
      if(identical(get_RVE_model_levels(model_results$RVE_model), get_RVE_model_levels(update_specific_model))){
        color_scale <- scale_color_manual(values = c("#33691E", "#BB9216"))
      }else{
        color_scale <-  scale_color_brewer(type = "qual", palette = "Dark2")
      }
    
      # Plot RVE trajectory, handling single or multiple RVE_levels
      RVE_trajectory %>%
        ggplot(aes(x = years_since_prior_vaccine, y = mean_R,
           color = if ("RVE_level" %in% names(RVE_trajectory) && length(unique(RVE_trajectory$RVE_level)) > 1) RVE_level else NULL,
           group = if ("RVE_level" %in% names(RVE_trajectory) && length(unique(RVE_trajectory$RVE_level)) > 1) RVE_level else NULL)) +
        geom_line(linewidth = 0.7) + 
        geom_errorbar(aes(ymin = lower_R, ymax = upper_R), width = 0.1) +
        geom_point(size = 2) +
        geom_hline(yintercept = 0, linetype = 2) +
        xlab("Prior vaccination\n(most recent to oldest)") +
        ylab("Effect of prior vaccination") +
        scale_x_continuous(breaks = unique(RVE_trajectory$years_since_prior_vaccine)) +
        scale_y_continuous(breaks = c(0, -0.2, -0.4, -0.6), limits = c(-0.65, 0)) +
        guides(color = guide_legend(title = ""), group = guide_legend(title = "")) +
        theme(legend.position = c(0.4, 0.3)) +
        color_scale +
        baseline_figure_settings

    }else{
      if(identical(get_RVE_model_levels(model_results$RVE_model), get_RVE_model_levels(no_decay_model))){
        model_results$posterior_summary_pars %>%
          filter(parameter == "r[1]") %>%
          ggplot(aes(x = RVE_level)) +
          geom_point(aes(y = mean)) +
          geom_linerange(aes(ymin = lower, ymax = upper)) +
          geom_hline(yintercept = 0, linetype = 2) +
          xlab("") +
          ylab("Posterior mean and 95% CrI")
      }else{
        stop("Function not implemented for this RVE model")
      }
    }
  }
}

plot_model_residuals <- function(data, model_results, censor_predictions){
  if(!censor_predictions){
    ylabel <- "Posterior average of residual"
    point_position <- position_identity()
  }else{
    ylabel <- "Posterior average of censored residuals\n(+ slight vertical jitter)"
    point_position <- position_jitter(width = 0, height = 0.1)
  }

  plot_data <- combine_data_and_predictions(data, model_results, censor_predictions) %>%
    mutate(mean_residual = Ht - H_hat_mean)

  # Calculate Pearson correlation coefficient and p-value between t and mean_residual
  cor_test <- cor.test(plot_data$t, plot_data$mean_residual)
  pearson_r <- cor_test$estimate
  p_value <- cor_test$p.value

  plot <- plot_data %>%
    ggplot(aes(x = t, y = mean_residual, color = factor(n_previous_vax))) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(position = point_position) +
    #facet_grid(subtype ~ year) +
    geom_smooth(aes(fill = factor(n_previous_vax)), show.legend = F) +
    geom_smooth(aes(color = NULL), color = "black") +
    xlab("Days since vaccination") +
    ylab(ylabel) +
    scale_color_discrete(name = "Number of prior vaccinations") +
    theme(legend.position = 'top') +
    annotate("text", x = Inf, y = Inf, label = paste0("italic(r) == ", round(pearson_r, 2)),
         hjust = 1.1, vjust = 1.5, size = 5, color = "black", parse = TRUE) +
    annotate("text", x = Inf, y = Inf, label = paste0("italic(p) == ", format.pval(p_value, digits = 2, eps = 1e-3)),
         hjust = 1.1, vjust = 3, size = 5, color = "black", parse = TRUE)

  return(plot)
}

plot_obs_vs_predicted_prevax_titers <- function(model_results, data, censor_predictions){

  plot_data <- data %>%
    mutate(Hpre_group_index = get_Hpre_group_index(data)) %>%
    select(individual, year, strain, subtype, n_previous_vax, Hpre_group_index, titer_type, Hpre) %>%
    unique() %>%
    left_join(
      model_results$posterior_summary_pars %>%
        filter(str_detect(parameter, "Hpre_continuous")) %>%
        mutate(Hpre_group_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
        select(Hpre_group_index, mean, lower, upper)
    )

  if(censor_predictions){
    ylab <- "Predicted pre-vaccination titer\n(censored posterior mean)"

    plot_data <- plot_data %>%
      mutate(across(all_of(c("mean", "lower", "upper")), ~log(censor_titers(titer = 2^.x, assay = titer_type), base = 2)))

  }else{
    ylab <- "Predicted pre-vaccination titer\n(posterior mean)"
  }
  
  plot_data %>%
    mutate(n_previous_vax = factor(n_previous_vax),
           subtype = factor(subtype, levels = subtype_levels)) %>%
    ggplot(aes(x = Hpre, y = mean, color = n_previous_vax)) +
    geom_point(position = position_jitter(width = ifelse(censor_predictions, 0.1, 0),
                                          height = ifelse(censor_predictions, 0.1, 0))) +
    facet_grid(year ~ subtype) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_smooth(se = F, method = "lm") +
    theme_cowplot() +
    scale_color_manual(name = "Prior vaccinations", values = year_4_nprior_vax_colors) +
    scale_x_continuous(breaks = log(original_dilutions, base = 2),
                       labels = ~2^.x) +
    scale_y_continuous(breaks = log(original_dilutions, base = 2),
                       labels = ~2^.x,
                       sec.axis = sec_axis(~.x, name = 'Study year\n', breaks = NULL)) + 
    xlab("Observed pre-vaccination titer") +
    ylab(ylab) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

}

plot_Hpre_means_and_sds <- function(model_results, input_data){

  input_tibble <- model_results$posterior_summary_pars %>%
    filter(str_detect(parameter, "Hpre_mean") |
           str_detect(parameter, "Hpre_sd")) %>%
    select(parameter, mean, lower, upper) %>%
    mutate(strain_npriorvax_group_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
    mutate(parameter = str_remove(parameter, "\\[[0-9]+\\]")) %>%
    left_join(
      input_data %>%
        mutate(strain_npriorvax_group_index = get_strain_npriorvax_group_index(input_data)) %>%
        select(strain, n_previous_vax, strain_npriorvax_group_index) %>%
        unique() %>%
        mutate(subtype = subtype_labeller(strain))
    ) %>%
    mutate(label = paste0(strain_labeller(strain), " - ", n_previous_vax)) %>%
    mutate(subtype = factor(subtype, levels = subtype_levels))

  data_points <- input_data %>%
        select(n_previous_vax, strain, subtype, Hpre) %>%
        # Just to match the right panel
        mutate(parameter = "Hpre_mean") %>%
        mutate(subtype = factor(subtype, levels = subtype_levels))

  # Empirical SDs computed from observed censored pre-vaccination titers
  empirical_sds <- data_points %>%
    group_by(strain, subtype, n_previous_vax) %>%
    summarise(Hpre_empirical_sd = sd(Hpre)) %>%
    ungroup() %>%
    # Just to match the right panel
    mutate(parameter = "Hpre_sd") %>%
    mutate(subtype = factor(subtype, levels = subtype_levels))

  pl <- input_tibble %>%
    ggplot(aes(x = n_previous_vax, y = mean)) +
    geom_boxplot( 
      data = data_points,
      aes(y = Hpre, group = n_previous_vax),
      outliers = F) +    
    geom_point(color = "red", size = 3) +
    geom_linerange(aes(ymin = lower, ymax = upper),
                   color = "red", linewidth = 2) +
    facet_nested(parameter ~ subtype + strain) +
    theme(panel.border = element_rect(color = "black")) +
    geom_point(
      data = data_points,
        aes(y = Hpre),
        position = position_jitter(height = 0, width = 0.1),
        alpha = 0.2
    ) +
    geom_point(data = empirical_sds,
               aes(y = Hpre_empirical_sd)) +
    ylab("Pre-vaccination titer") +
    xlab("Number of prior vaccinations")

  return(pl)
}

# The shell script for LOO counts the number of units based on processed_data.csv before to generate a job array.
# Some units may be filtered out from input_data, rendering some array indices irrelevant.
# This skips array indices in excess of the number of actual LOO units.
skip_missing_LOO_indices <- function(){
  
  if (holdout_index > 0) {
    max_index <- ifelse(loglik_unit == 1,
                        max(input_list_real$Hpre_group_index),
                        max(input_list_real$individual))

    if (holdout_index > max_index) {
      message(sprintf("Array index exceeds effective number of LOO units. Ending job."))
      # Remove temporary lock file if it exists
      if (file.exists(LOO_tmp_path)) {
        file.remove(LOO_tmp_path)
      }
      quit(save = "no")
    }
  }
}

# Retrieves model information from directory containing fit_specs.R
retrieve_model_info <- function(model_dir){

  fit_specs_file <- paste0(model_dir, "/fit_specs.R")

  if(file.exists(fit_specs_file)){
    source(fit_specs_file, local = environment())

    RVE_model <- fit_specs$RVE_model

    if(is.null(fit_specs$RVE_model)){
      model_name <- "no_RVE_model"
    }else{
      for(name in names(models_list)){
        if(identical(get_RVE_model_levels(RVE_model), get_RVE_model_levels(models_list[[name]]))){
          model_name <- name
          break
        } 
      }
    }

    # Revise this function if we fit models with more than 1 kind of individual effect
    stopifnot((fit_specs$include_upeak_shared + fit_specs$include_upeak_subtype + fit_specs$include_upeak_year) <= 1)

    individual_effects <- case_when(
      fit_specs$include_upeak_shared == 1 ~ "global",
      fit_specs$include_upeak_subtype == 1 ~ "subtype-specific",
      fit_specs$include_upeak_year == 1 ~ "year-specific",
      !fit_specs$include_upeak_shared & !fit_specs$include_upeak_subtype & !fit_specs$include_upeak_year ~ "none"
    )

    model_info <- tibble(
      model = model_name,
      individual_effects = individual_effects,
      Hpre_priors = as.logical(fit_specs$include_Hpre_priors),
      NAI_infections = as.logical(fit_specs$include_NAI_infections),
      include_k_peak_flu = as.logical(fit_specs$include_k_peak_flu),
      include_age = as.logical(fit_specs$include_age),
      include_sex = as.logical(fit_specs$include_sex),
      include_BMI = as.logical(fit_specs$include_BMI),
      include_smoking = as.logical(fit_specs$include_smoking),
      include_asthma = as.logical(fit_specs$include_asthma)
    )

  }else{
    # If fit_specs.R doesn't exist, check this is the null model
    stopifnot(str_detect(model_dir, "null_specific_model/no_individual_effects"))

    model_info <- tibble(
      model = "null_specific_model",
      individual_effects = "none",
      Hpre_priors = NA,
      NAI_infections = F,
      include_k_peak_flu = F,
      include_age = F,
      include_sex = F,
      include_BMI = F,
      include_smoking = F,
      include_asthma = F
    )
    
  }
  return(model_info)

}


assign_model_labels <- function(data){

  # This sets the order and labels of models in model comparison plots

  model_labels <- tibble(
    model = c("null_specific_model",
              "no_RVE_model",
              "no_decay_model",
              "fully_shared_model",
              "update_specific_model",
              "subtype_specific_model",
              "fully_flexible_RVE_model")) %>%
    mutate(model_type = case_when(
        str_detect(model, "null") ~ "Null",
        str_detect(model, "no_RVE") ~ "No RVE",
        T ~ "RVE"
    )) %>%
    mutate(model_label = case_when(
      str_detect(model, "null") ~ "Null",
      str_detect(model, "no_RVE") ~ "No RVE",
      str_detect(model, "no_decay_model") ~ "Non-saturating\nshared RVE",
      model == "fully_shared_model" ~ "Saturating\nshared RVE",
      str_detect(model, "update_specific_model") ~ "Saturating shared\nRVE with update\nmodifier",
      str_detect(model, "subtype_specific_model") ~ "Saturating\nsubtype-specific\nRVE",
      str_detect(model, "fully_flexible_RVE_model") ~ "Fully flexible RVE"
  ))

  individual_effect_levels <- c("none", "subtype-specific", "year-specific", "global")

  data %>%
    left_join(model_labels) %>%
    mutate(model = factor(model, levels = model_labels$model),
           model_type = factor(model_type, levels = unique(model_labels$model_type)),
           model_label = factor(model_label, levels = model_labels$model_label)) %>%
    mutate(individual_effects = factor(individual_effects, levels = individual_effect_levels)) %>%
    select(matches("model"), everything())
 
}

summarize_model_residuals <- function(model_residuals, by_subtype){

    grouping_vars <- c("model", "model_type", "individual_effects", "Hpre_priors", "NAI_infections",  "model_label")

    grouping_vars <- c(grouping_vars, model_residuals %>% select(matches("include")) %>% names())

    if(by_subtype){
      grouping_vars <- c(grouping_vars, "subtype")
    }

    model_residuals %>%
      group_by(across(any_of(grouping_vars))) %>%
      summarise(mean_absolute_error = mean(abs(censored_mean_residual_post)),
               lower = quantile(abs(censored_mean_residual_post), 0.25),
                upper = quantile(abs(censored_mean_residual_post), 0.75),
              fraction_within_2fold = sum(abs(censored_mean_residual_post) <= 1) / n(),
              fraction_identical = sum(abs(censored_mean_residual_post) == 0) / n(),
              n = n()) %>%
      # Computing lower and upper limits of an approximate 95% binomial CI
      mutate(fraction_within_2fold_SD = sqrt(fraction_within_2fold * (1 - fraction_within_2fold)),
             fraction_identical_SD = sqrt(fraction_identical * (1 - fraction_identical))) %>%
      mutate(fraction_within_2fold_lower = fraction_within_2fold - 1.96 * fraction_within_2fold_SD / sqrt(n),
             fraction_within_2fold_upper = fraction_within_2fold + 1.96 * fraction_within_2fold_SD / sqrt(n),
             fraction_identical_lower = fraction_identical - 1.96 * fraction_identical_SD / sqrt(n),
             fraction_identical_upper = fraction_identical + 1.96 * fraction_identical_SD / sqrt(n)) %>%
      ungroup() %>%
      select(-matches("_SD")) %>%
      mutate(percent_within_2fold_label = round(fraction_within_2fold * 100),
             percent_identical_label = round(fraction_identical * 100))
}

make_abs_error_pl <- function(model_residuals_summary, xvar){
  pl <- model_residuals_summary %>%
    ggplot(aes(x = .data[[xvar]], y = mean_absolute_error, color = Hpre_priors, group = Hpre_priors)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
                  position = position_dodge(width = 0.5),
                  linewidth = boxplot_line_width) +
    baseline_figure_settings +
    ylab("Absolute censored error\n(mean, IQR)")

  if(length(unique(model_residuals_summary$Hpre_priors)) == 1){
    pl <- pl + scale_color_manual(values = "gray40") +
      theme(legend.position = "none")
  }else{
    pl <- pl + 
      scale_color_brewer(
      name = "Priors on latent pre-vaccination titers",
      type = "qual",
      palette = "Set1",
      na.value = "gray60") +
      theme(legend.position = "top")
  }
  return(pl)
}

make_percent_within_2fold_pl <- function(model_residuals_summary, xvar){
  pl <- model_residuals_summary %>%
    ggplot(aes(x = .data[[xvar]], y = fraction_within_2fold * 100, fill = Hpre_priors, group = Hpre_priors)) +
    geom_col(position = position_dodge(width = 0.5), width = default_col_width) +
    geom_errorbar(aes(ymin = fraction_within_2fold_lower * 100, ymax = fraction_within_2fold_upper * 100), position = position_dodge(width = 0.5),
                  width = 0.3, linewidth = boxplot_line_width) + 
    geom_text(aes(y = fraction_within_2fold_upper * 100 + 3,
              label = percent_within_2fold_label),
              position = position_dodge(width = 0.5),
              size = default_figure_font_size, size.unit = "pt") +
    ylab("Percent of observations within\n 2-fold of prediction") +
    baseline_figure_settings +
    theme(axis.text.x = element_text(size = small_font_size))

  if(length(unique(model_residuals_summary$Hpre_priors)) == 1){
    pl <- pl + scale_fill_manual(values = "gray40") +
      theme(legend.position = "none")
  }else{
    pl <- pl + 
      scale_fill_brewer(
      name = "Priors on latent pre-vaccination titers",
      type = "qual",
      palette = "Set1",
      na.value = "gray60") +
      theme(legend.position = "top")
  }
  
  return(pl)

}

make_percent_identical_pl <- function(model_residuals_summary, xvar){
  pl <- model_residuals_summary %>%
    ggplot(aes(x = .data[[xvar]], y = fraction_identical * 100, fill = Hpre_priors, group = Hpre_priors)) +
    geom_col(position = position_dodge(width = 0.5), width = default_col_width) +
    geom_errorbar(aes(ymin = fraction_identical_lower * 100, ymax = fraction_identical_upper * 100), position = position_dodge(width = 0.5),
                  width = 0.3, linewidth = boxplot_line_width) +
    geom_text(aes(y = fraction_identical_upper * 100 + 3, label = percent_identical_label), ,
              position = position_dodge(width = 0.5),
              size = default_figure_font_size, size.unit = "pt") +
    ylab("Percent of observations with\n identical prediction") +
    baseline_figure_settings +
    theme(axis.text.x = element_text(size = small_font_size))

  if(length(unique(model_residuals_summary$Hpre_priors)) == 1){
    pl <- pl + scale_fill_manual(values = "gray40") +
      theme(legend.position = "none")
  }else{
    pl <- pl + 
      scale_fill_brewer(
      name = "Priors on latent pre-vaccination titers",
      type = "qual",
      palette = "Set1",
      na.value = "gray60") +
      theme(legend.position = "top")
  }
}

plot_fixed_effects_comparison <- function(model_residuals, by_subtype){

  model_residuals_summary <-  summarize_model_residuals(model_residuals = model_residuals, by_subtype = by_subtype)

  abs_error_pl <- make_abs_error_pl(model_residuals_summary, xvar = "model_label") +
    xlab("Model")

  percent_within_2fold_plot <- make_percent_within_2fold_pl(model_residuals_summary, xvar = "model_label") +
    xlab("Model")

  percent_identical_plot <- make_percent_identical_pl(model_residuals_summary, xvar = "model_label")  +
    xlab("Model")
 
  if(by_subtype){
    abs_error_pl <- abs_error_pl + facet_wrap("subtype")
    percent_within_2fold_plot <- percent_within_2fold_plot + facet_wrap("subtype", nrow = 1)
    percent_identical_plot <- percent_identical_plot + facet_wrap("subtype", nrow = 1)
  }
  percents_plot <- plot_grid(
    percent_within_2fold_plot + xlab("") + theme(axis.text.x = element_blank()),
    percent_identical_plot + theme(legend.position = "none"), nrow = 2)

  # Return a list of plots
  return(list(
    abs_error_plot = abs_error_pl,
    percents_plot = percents_plot
  ))

}

plot_individual_effects_comparison <- function(model_residuals, by_subtype){
  # Similar code to function above, but not worth generalizing.

  model_residuals_summary <- summarize_model_residuals(model_residuals = model_residuals, by_subtype = by_subtype)

  # If not plotting by subtype, change newline structure of some of the model labels
  model_residuals_summary <- model_residuals_summary %>%
    mutate(model_label = str_replace_all(as.character(model_label), "\n", " "))

  abs_error_pl <- make_abs_error_pl(model_residuals_summary, xvar = "individual_effects") +
    xlab("Individual effects")

  percent_within_2fold_plot <- make_percent_within_2fold_pl(model_residuals_summary, xvar = "individual_effects") +
    xlab("Individual effects")

  percent_identical_plot <- make_percent_identical_pl(model_residuals_summary, xvar = "individual_effects")

  if(by_subtype){
    abs_error_pl <- abs_error_pl + facet_grid(model_label ~ subtype)
    percent_within_2fold_plot <- percent_within_2fold_plot + facet_grid(model_label ~ subtype)
    percent_identical_plot <- percent_identical_plot + facet_grid(model_label ~ subtype)
  }else{
    abs_error_pl <- abs_error_pl + facet_wrap("model_label")
    percent_within_2fold_plot <- percent_within_2fold_plot + facet_wrap("model_label", nrow = 1)
    percent_identical_plot <- percent_identical_plot + facet_wrap("model_label", nrow = 1)
  }

  percents_plot <- plot_grid(
    percent_within_2fold_plot + xlab("") + theme(axis.text.x = element_blank()),
    percent_identical_plot + xlab("Individual effects") +
      theme(legend.position = "none",
            axis.text.x = element_text(size = default_figure_font_size, angle = 45, hjust = 1)),
      nrow = 2
  )

  return(list(
    abs_error_plot = abs_error_pl,
    percents_plot = percents_plot
  ))
}

plot_covariate_set_comparison <- function(model_residuals, by_subtype){
  model_residuals_summary <- summarize_model_residuals(model_residuals = model_residuals, by_subtype = by_subtype)

  model_residuals_summary <- model_residuals_summary %>%
    mutate(predictors = case_when(
      if_all(starts_with("include_"), ~ .x == TRUE) ~ "Default predictors",
      if_all(starts_with("include_"), ~ .x == FALSE) ~ "Pre-vaccination titer\nand vaccination history only",
      if_all(all_of(c("include_k_peak_flu", "include_BMI", "include_smoking")), ~ .x == FALSE) ~ "Significant predictors\nonly",
    ))

  model_residuals_summary$predictors <- factor(model_residuals_summary$predictors,
    levels = c(
      "Default predictors",
      "Pre-vaccination titer\nand vaccination history only",
      "Significant predictors\nonly"
    ))

  abs_error_pl <- make_abs_error_pl(model_residuals_summary, xvar = "predictors") +
    xlab("Predictor set")

  percent_within_2fold_plot <- make_percent_within_2fold_pl(model_residuals_summary, xvar = "predictors") +
    xlab("Predictor set")

  percent_identical_plot <- make_percent_identical_pl(model_residuals_summary, xvar = "predictors") +
    xlab("Predictor set")

  if(by_subtype){
    abs_error_pl <- abs_error_pl + facet_wrap("subtype")
    percent_within_2fold_plot <- percent_within_2fold_plot + facet_wrap("subtype", nrow = 1)
    percent_identical_plot <- percent_identical_plot + facet_wrap("subtype", nrow = 1)
  }

  percents_plot <- plot_grid(
    percent_within_2fold_plot + xlab("") + theme(axis.text.x = element_blank()),
    percent_identical_plot + xlab("Predictor set") + theme(legend.position = "none"), nrow = 2
  )

  return(list(
    abs_error_plot = abs_error_pl,
    percents_plot = percents_plot
  ))
}

# ========= RVE models =========

# RVE model objects consist of lists of two elements
# $form indicates, for each combination of prior year / current year / subtype
# Which "level" repeat vaccination effect applies (before accounting for any time decay)
# For combinations with the same level, the prior vaccine has the same effect on the current
# vaccine's response (again, before accounting for any decay)

# Because the Stan model uses vaccination status in each season in wide format,
# we use a special level to make it so that current/future years don't count
# toward repeat vaccination effects
RVE_special_level <- "NA"

# This function retrieves the unique levels of any given RVE model (for ordering purposes)
get_RVE_model_levels <- function(RVE_model){c(unique(RVE_model$form$RVE_level), RVE_special_level)}

# Fully flexible RVE model: 
fully_flexible_RVE_model <-
  list(
    form = vax_strain_pairs %>%
      # In the fully flexible model, each combination of subtype, current year and prior year
      # gets its own RVE effect (with no separate time-decay, which would not be identifiable)
      mutate(RVE_level = paste(previous_year, year, subtype, sep = '-')) %>%
      # Setting some default values for synthetic value experiments
      mutate(
        default_synthetic_value = case_when(
          RVE_level == "1-2-H3N2" ~ -0.2,
          RVE_level == "1-3-H3N2" ~ -0.1,
          RVE_level == "2-3-H3N2" ~ -0.3,
          RVE_level == "1-4-H3N2" ~ -0.05,
          RVE_level == "2-4-H3N2" ~ -0.1,
          RVE_level == "3-4-H3N2" ~ -0.3,
          RVE_level == "1-2-H1N1" ~ -0.1,
          RVE_level == "1-3-H1N1" ~ -0.05,
          RVE_level == "2-3-H1N1" ~ -0.25,
          RVE_level == "1-4-H1N1" ~ 0,
          RVE_level == "2-4-H1N1" ~ -0.15,
          RVE_level == "3-4-H1N1" ~ -0.3,
          RVE_level == "1-2-B/Yamagata" ~ -0.2,
          RVE_level == "1-3-B/Yamagata" ~ -0.1,
          RVE_level == "2-3-B/Yamagata" ~ -0.2,
          RVE_level == "1-4-B/Yamagata" ~ 0,
          RVE_level == "2-4-B/Yamagata" ~ -0.1,
          RVE_level == "3-4-B/Yamagata" ~ -0.2,
          RVE_level == "1-2-B/Victoria" ~ -0.3,
          RVE_level == "1-3-B/Victoria" ~ -0.1,
          RVE_level == "2-3-B/Victoria" ~ -0.2,
          RVE_level == "1-4-B/Victoria" ~ 0,
          RVE_level == "2-4-B/Victoria" ~ -0.1,
          RVE_level == "3-4-B/Victoria" ~ -0.3
        )
      ),
      RVE_decay = F
  )

# Strain-pair decay model
strain_pair_model <- list(
  form = vax_strain_pairs %>%
      # In the strain-pair decay model, each combination of current and past strains gets its own
      # baseline repeat vaccination effect.
      # The *baseline* effect is the same regardless how far apart the 2 strains were in time.
      # but the effect is then modified to decay exponentially with time
      # (controlled by parameter gamma and not reflected here)
      mutate(RVE_level = pair_label) %>%
      # Setting some default values for synthetic value experiments
      mutate(
        default_synthetic_value = case_when(
          RVE_level == "HK-Camb" ~ -0.1,
          RVE_level == "HK-Darwin" ~ -0.1,
          RVE_level == "Camb-Darwin" ~ -0.2,
          RVE_level == "Darwin-Darwin" ~ -0.3,
          RVE_level == "Hawaii-Wis588" ~ -0.1,
          RVE_level == "Wis588-Wis588" ~ -0.3,
          RVE_level == "Hawaii-Wis67" ~ -0.1,
          RVE_level == "Wis588-Wis67" ~ -0.2,
          RVE_level == "Wash-Wash" ~ -0.3,
          RVE_level == "Wash-Austria" ~ -0.1,
          RVE_level == "Austria-Austria" ~ -0.2,
          RVE_level == "Phuk-Phuk" ~ -0.2
        )
      ),
    RVE_decay = T
)

subtype_specific_model <- list(
  form = vax_strain_pairs %>%
      # In the subtype-specific RVE model, each subtype gets its own baseline repeat vaccination effect
      # the baseline effect decays over time with a shared rate
      # (controlled by parameter gamma and not reflected here)
      mutate(RVE_level = subtype) %>%
      # Setting some default values for synthetic value experiments
      mutate(
        default_synthetic_value = case_when(
          RVE_level == "H3N2" ~ 0,
          RVE_level == "H1N1" ~ -0.1,
          RVE_level == "B/Victoria" ~ -0.3,
          RVE_level == "B/Yamagata" ~ -0.2
        )
      ),
    RVE_decay = T
)

# Fully-shared  model
fully_shared_model <- list(
  # In the fully-shared model, there's a single baseline RVE
  # (therefore a single RVE level)
  # the baseline effect decays over time with a shared rate
  # (controlled by parameter gamma and not reflected here)
  form = vax_strain_pairs %>%
      mutate(RVE_level = "fully_shared_RVE") %>%
      # Setting some default values for synthetic value experiments
      mutate(
        default_synthetic_value = -0.2
      ),
    RVE_decay = T
)

# No-decay model
no_decay_model <- list(
  # Like the fully-shared model, but the single baseline RVE
  # does not decay for older prior vaccinations.
  form = vax_strain_pairs %>%
      mutate(RVE_level = "non_decay_RVE") %>%
      # Setting some default values for synthetic value experiments
      mutate(
        default_synthetic_value = -0.1
      ),
    RVE_decay = F
)

# Update-specific baseline model
update_specific_model <- list(
  # In the update-specific baseline model, there are separate baseline effects
  # when the strains in strain_pair are the same and when they are not
  # the baseline effect also decays over time with a share rate
  # (controlled by parameter gamma and not reflected here)
  form = vax_strain_pairs %>%
      separate(strain_pair, into = c("strain1", "strain2"), sep = '-', remove = F) %>%
      mutate(RVE_level = ifelse(strain1 == strain2, "same_strain", "different_strain")) %>%
      # Setting some default values for synthetic value experiments
      mutate(
        default_synthetic_value = case_when(
          RVE_level == "same_strain" ~ -0.2,
          RVE_level == "different_strain" ~ -0.1
        )
      ),
    RVE_decay = T
)

# Collect models in a list
models_list <- list(fully_flexible_RVE_model = fully_flexible_RVE_model, strain_pair_model = strain_pair_model,
                    subtype_specific_model = subtype_specific_model, fully_shared_model = fully_shared_model,
                    no_decay_model = no_decay_model, update_specific_model = update_specific_model)


compute_nonmonotonicity_condition <- function(model_results, input_data){

  # Only implemented for the fully shared RVE model
  stopifnot(all(model_results$RVE_model$form$RVE_level == "fully_shared_RVE"))

  data_and_predictions <- input_data %>%
    mutate(Hpre_group_index = get_Hpre_group_index(input_data)) %>%
    select(individual, year, strain, subtype, treatment, n_previous_vax, Hpre_group_index, titer_type, Hpre) %>%
    unique() %>%
    mutate(timepoint = 0) %>%
    left_join(
      model_results$posterior_summary_pars %>%
        filter(str_detect(parameter, "Hpre_continuous")) %>%
        mutate(Hpre_group_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
        select(Hpre_group_index, mean, lower, upper)
    )

  pairwise_comparisons_predicted <- run_pairwise_wilcoxon(
    response_data = data_and_predictions %>% rename(log2_titer = "mean"),
    response_var = "log2_titer", p.adjust.method = "holm") %>%
    mutate(Hpre_treatment_1_predicted = log2(GMT_treatment_1),
           Hpre_treatment_2_predicted = log2(GMT_treatment_2)) %>%
    rename(p_value_predicted = p_value) %>%
    select(-matches("GMT_"))

  pairwise_comparisons_obs <- run_pairwise_wilcoxon(
    response_data = data_and_predictions %>% rename(log2_titer = "Hpre"),
    response_var = "log2_titer", p.adjust.method = "holm") %>%
    mutate(Hpre_treatment_1_obs = log2(GMT_treatment_1),
           Hpre_treatment_2_obs = log2(GMT_treatment_2)) %>%
    rename(p_value_obs = p_value) %>%
    select(-matches("GMT_"))

  pairwise_comparisons <- pairwise_comparisons_predicted %>%
    left_join(pairwise_comparisons_obs)  %>%
    # Comparisons between first and second time vaccinees
      filter(treatment_1 %in% c("VV", "PVV", "PPVV"),
             treatment_2 %in% c("PV", "PPV", "PPPV"))

  ab_ceiling_effect <- model_results$posterior_summary_pars %>%
    filter(str_detect(parameter, "b_peak")) %>%
    select(strain, mean) %>%
    rename(b = mean)

  intercept <- model_results$posterior_summary_pars %>%
    filter(str_detect(parameter, "a\\[")) %>%
    select(strain, mean) %>%
    rename(a = mean)

  repeated_vaccination_effect <- model_results$posterior_summary_pars %>%
    filter(parameter == "r[1]") %>%
    pull(mean)

  stopifnot(length(repeated_vaccination_effect) == 1)

  # Check that second-time vaccinees are in "treatment_1" and first-time vaccinees in "treatment_2"
  stopifnot(all(str_count(pairwise_comparisons$treatment_1, "V") == 2))
  stopifnot(all(str_count(pairwise_comparisons$treatment_2, "V") == 1))

  return(
    pairwise_comparisons %>%
      left_join(ab_ceiling_effect) %>%
      left_join(intercept) %>%
      mutate(repeated_vaccination_effect) %>%
      mutate(non_monotonicity_threshold = (1 - b) * (Hpre_treatment_1_predicted - Hpre_treatment_2_predicted) / (b * (a / b - Hpre_treatment_1_predicted))) %>%
      mutate(condition_met = abs(repeated_vaccination_effect) < non_monotonicity_threshold) %>%
      arrange(year) %>%
      select(year, strain, Hpre_treatment_1_predicted, Hpre_treatment_2_predicted, b, a, repeated_vaccination_effect, non_monotonicity_threshold, condition_met) %>%
      rename(predicted_Hpre_second_times = Hpre_treatment_1_predicted,
             predicted_Hpre_first_times = Hpre_treatment_2_predicted)
  )
    
}