library(tidyverse)
library(rstan)
library(cowplot)
theme_set(theme_cowplot())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
source("R/base_functions.R")
source("R/bayesian_model_functions.R")

start_time <- Sys.time()

# Read fit specification directory from argument, look for fit_specs.R in that directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript fit_bayesian_model_real_data.R <fit_specs_dir> [loglik_unit] [holdout_index]")
}
fit_specs_dir <- args[1]
# Optional args with defaults
# Unit for LOO (1 = individual/strain/year, 2 = individual)
loglik_unit <- if (length(args) >= 2) as.integer(args[2]) else 1 
# Index of held-out unit (if 0, no unit is held out)
holdout_index <- if (length(args) >= 3) as.integer(args[3]) else 0


# Validate
if (is.na(loglik_unit) || !(loglik_unit %in% c(1, 2))) {
  stop("loglik_unit must be 1 or 2")
}
if (is.na(holdout_index) || holdout_index < 0) {
  stop("holdout_index must be an integer >= 0")
}

# Early exit if LOO output already exists (avoid running Stan)
if (holdout_index > 0) {
  LOO_subdir <- if (loglik_unit == 1) "LOO_Hpre_group" else "LOO_individual"
  LOO_dir_path <- file.path(fit_specs_dir, LOO_subdir)
  LOO_out_path <- file.path(LOO_dir_path, paste0("holdout_", holdout_index, ".csv"))
  
  # Check for temp file to avoid running simultaneous jobs for the same specs and holdout unit
  LOO_tmp_path <- paste0(LOO_out_path, ".tmp")
  
  if (file.exists(LOO_out_path) || file.exists(LOO_tmp_path)) {
    message(sprintf("LOO output already exists or job already running. Ending redundant job.", LOO_out_path))
    quit(save = "no")
  }
  
  # Create temp file
  dir.create(LOO_dir_path, recursive = TRUE, showWarnings = FALSE)
  file.create(LOO_tmp_path)
}

fit_specs_path <- file.path(fit_specs_dir, "fit_specs.R")
if (!file.exists(fit_specs_path)) {
  stop(paste("fit_specs.R not found in directory:", fit_specs_dir))
}
source(fit_specs_path)

# Read processed data
input_data <- read_csv("bayesian_fits_real_data/processed_data.csv")


if(!is.null(fit_specs$RVE_model)){
  include_RVE <- T
}else{
  # If fit_specs$RVE_model is NULL, fit model without RVEs
  include_RVE <- F
  # In this case, we'll pass an arbitrary RVE model to fit_specs
  fit_specs$RVE_model <- fully_shared_model
  # (This is necessary but won't matter, since include_RVE = F will be passed to stan)
}

# Assign RVE levels based on the specified model in fit_specs list
input_data <- input_data %>% assign_RVE_levels(RVE_model = fit_specs$RVE_model)


# We're disabling kinetics by default (so fitting to day 30 only)
include_kinetics_real_data <- 0
if(include_kinetics_real_data == F){
  input_data <- input_data %>%
    filter(timepoint == 30)
}

input_data <- input_data %>%
  mutate(strain = factor(strain, levels = unique(input_data$strain)))

save(input_data, file = paste0(fit_specs_dir, "/input_data.rds"))

input_list_real <- prepare_stan_input(
  data = input_data,
  censoring = 1, # Censoring is always true for real data
  include_kinetics = include_kinetics_real_data,
  include_RVE = include_RVE,
  include_upeak_shared = fit_specs$include_upeak_shared,
  include_upeak_subtype = fit_specs$include_upeak_subtype,
  include_upeak_year = fit_specs$include_upeak_year,
  include_k_peak_flu = fit_specs$include_k_peak_flu,
  include_k_peak_cov2 = fit_specs$include_k_peak_cov2,
  RVE_model = fit_specs$RVE_model,
  loglik_unit = loglik_unit,
  holdout_index = holdout_index,
  include_Hpre_priors = fit_specs$include_Hpre_priors,
  include_NAI_infections = fit_specs$include_NAI_infections,
  include_time_since_vax = fit_specs$include_time_since_vax,
  include_sex = fit_specs$include_sex,
  include_age = fit_specs$include_age,
  include_BMI = fit_specs$include_BMI,
  include_smoking = fit_specs$include_smoking,
  include_asthma = fit_specs$include_asthma
)

# The shell script for LOO counts the number of units based on processed_data.csv before to generate a job array.
# Some units may be filtered out from input_data, rendering some array indices irrelevant.
# This skips array indices in excess of the number of actual LOO units.
skip_missing_LOO_indices()

n_iterations <- 5000
n_chains <- 4

# Require these diagnostics for the log-posterior, otherwise repeat
# not required when performing LOO (holdout_index != 0)
ess_cutoff <- if (holdout_index != 0) 0 else n_chains * 100
rhat_cutoff <- if (holdout_index != 0) Inf else 1.01

repeat {
  model_fit <- stan(file = "R/linear_model.stan",
                    data = input_list_real, iter = n_iterations,
                    chains = n_chains)

  # Generate summary tibble
  model_results <- process_stan_fit(fit = model_fit, data = input_data,
                                    RVE_model = fit_specs$RVE_model)

  lp_bulk_ESS <- model_results$diagnostics %>% 
    filter(par == "lp__") %>%
    pull(Bulk_ESS)

  lp_rhat <- model_results$diagnostics %>% 
    filter(par == "lp__") %>%
    pull(Rhat)

  if((lp_bulk_ESS >= ess_cutoff) && round(lp_rhat, 2) <= rhat_cutoff){
    break
  }else{
      print(paste0("log-posterior ESS: ", lp_bulk_ESS, "; Rhat: ", round(lp_rhat, 3), ". Repeating"))
  }
}

end_time <- Sys.time()

# If no LOO was performed, export detailed results
if(holdout_index == 0){
  save(model_fit, file = paste0(fit_specs_dir, "/model_fit.rds"))
  save(model_results, file = paste0(fit_specs_dir, "/model_results.rds"))
}else{
  
  # If we held out a combination of individual, strain, year based on holdout_index,
  # identify what combination that was
  if(loglik_unit == 1){
    holdout_info <- tibble(Hpre_group = input_list_real$Hpre_group_index,
                           individual = input_list_real$individual,
                           measurement_replicate = input_list_real$measurement_replicate,
                           timepoint = input_data$timepoint,
                           treatment = input_data$treatment,
                           Hpre = input_list_real$Hpre,
                           Ht = input_list_real$Ht, 
                           strain = input_list_real$strain,
                           year = input_list_real$year) %>%
           mutate(observation_index = 1:n()) %>%
           unique() %>%
           filter(Hpre_group == holdout_index) %>%
           mutate(strain = levels(input_data$strain)[strain]) %>%
           select(Hpre_group, individual, timepoint, year, treatment, strain, everything())
  }else{
    stop("Post-processing not implemented if holding out an entire individual.")
    holdout_info <- tibble(individual = holdout_index)

  }
  
  # If in the future doing this for multiple timepoints, revise this function.
  stopifnot(length(unique(holdout_info$timepoint)) == 1)

  # If held-out unit has replicated measurements, keep only one
  # (This doesn't affect LOO LPD calculations, just simplifies computing additional metrics of predictive ability below)
  holdout_info <- holdout_info %>%
    mutate(has_replicate_measurement = any(measurement_replicate > 1)) %>%
    filter(measurement_replicate == 1) %>%
    select(-measurement_replicate, -timepoint)

  stopifnot(nrow(holdout_info) == 1)

  names(holdout_info) <- paste0("LOO_", names(holdout_info))
  
  LOO_results <- holdout_info %>%
    # Annotate with posterior mean, lower, upper LPD for this unit
    cross_join(
      model_results$posterior_summary_pars %>%
        filter(parameter %in% c("LOO_LPD", "LOO_LPD_pre", "LOO_LPD_post")) %>%
        select(parameter, mean, lower, upper) %>%
        pivot_wider(
          names_from = parameter,
          values_from = c(mean, lower, upper),
          names_glue = "{parameter}_{.value}"
        ) %>%
        # Annotate with diagnostics for the posterior distribution of the LOO fit
        cross_join(
          model_results$diagnostics %>%
            filter(par == "lp__") %>%
            select(-par) %>%
            rename(lp_valid = valid, lp_Rhat = Rhat,
                   lp_Bulk_ESS = Bulk_ESS, lp_Tail_ESS = Tail_ESS)
        ) 
    ) %>%
        mutate(start_time = start_time, end_time = end_time)

  # Other measures of predictive performance for post-vaccination titers
  predictive_performance_postvax <- tibble(
    Hpre_group_index = input_list_real$Hpre_group_index,
    titer_type = input_data$titer_type) %>%
    mutate(observation_index = 1:n()) %>%
    filter(Hpre_group_index == holdout_index) %>%
    # Take posterior sample of predictive values for the held-out post-vaccination titer
    left_join(model_results$predicted_value_samples %>%
                select(observation_index, H_hat)) %>%
    rename(LOO_H_hat = H_hat) %>%
    mutate(LOO_Ht = holdout_info$LOO_Ht) %>%
    # Censor those predictions
    mutate(censored_LOO_H_hat = log2(censor_titers(titer = 2^LOO_H_hat, assay = titer_type))) %>%
    summarise(
      # Mean residual 
      mean_residual_post = mean(LOO_Ht - LOO_H_hat),
      mean_censored_residual_post = mean(LOO_Ht - censored_LOO_H_hat),
      # Mean prediction 
      mean_prediction_post = mean(LOO_H_hat),
      # Mean censored prediction 
      mean_censored_prediction_post = mean(censored_LOO_H_hat),
      # Probability that censored prediction is equal to the held-out post-vaccination titers
      prob_LOO_postvax_equal_prediction = mean(abs(censored_LOO_H_hat - LOO_Ht) < 1e-6),
      # Probability that the held-out post-vac titer is within 1 unit of censored prediction.
      prob_LOO_postvax_within_1_unit = mean(abs(censored_LOO_H_hat - LOO_Ht) <= 1 + 1e-6)
    )
  
  # Same measures of predictive performance for pre-vaccination titers
  # Take sample of predicted post-vaccination titers for this Hpre group
  titer_type = unique(input_data$titer_type[input_list_real$Hpre_group_index == holdout_index])


  # Note: we compute predictive performance for pre-vaccination titers, but note that the observed
  # pre-vaccination titer of the held-out unit is included in the training data.
  predictive_performance_prevax <- model_results$parameter_samples %>%
    filter(str_detect(parameter, "Hpre_continuous")) %>%
    select(parameter, value) %>%
    mutate(Hpre_group = str_extract(parameter, "[0-9]+") %>% as.integer()) %>%
    filter(Hpre_group == holdout_index) %>%
    rename(Hpre_predicted = value) %>%
    select(Hpre_predicted) %>%
    mutate(titer_type = titer_type) %>%
    # Censor those predictions 
    mutate(censored_Hpre_predicted = log2(censor_titers(titer = 2^Hpre_predicted,
                                                        assay = titer_type))) %>%
    # Annotate observed pre-vaccination titer
    mutate(LOO_Hpre = LOO_results$LOO_Hpre) %>%
    # Compute predictive measures 
    summarise(
      mean_residual_pre = mean(LOO_Hpre - Hpre_predicted),
      mean_censored_residual_pre = mean(LOO_Hpre - censored_Hpre_predicted),
      mean_prediction_pre = mean(Hpre_predicted),
      mean_censored_prediction_pre = mean(censored_Hpre_predicted),
      prob_LOO_prevax_equal_prediction = mean(abs(censored_Hpre_predicted - LOO_Hpre) < 1e-6),
      prob_LOO_prevax_within_1_unit = mean(abs(censored_Hpre_predicted - LOO_Hpre) <= 1 + 1e-6)
    )
  
  LOO_results <- LOO_results %>%
    bind_cols(predictive_performance_postvax) %>%
    bind_cols(predictive_performance_prevax) %>%
    select(matches("LOO"), matches("mean_"), everything())

  # LOO-CV will double as a test of the function that calculates censored normals in Stan
  # (For convenience, we implement this only for samples with a single measurement)
  if(holdout_info$LOO_has_replicate_measurement == F){

    # Log2 of smallest and highest *observed* dilutions
    L = log2(2 * unique(input_data$undetectable_value[input_list_real$Hpre_group_index == holdout_index]))
    U = log2(as.integer(max(original_dilutions)))

    # Take posterior sample of predicted post-vaccination titers 
    H_hat <- rstan::extract(model_fit, pars = paste0("H_hat[", holdout_info$LOO_observation_index, "]"))[[1]] %>% as.numeric()
    
    # Take posterior sample of measurement error sds (sigma_FRNT or sigma_HAI)
    sigma <- rstan::extract(model_fit, pars = paste0("sigma_", titer_type))[[1]] %>% as.numeric()

    # Take posterior sample of LPDs
    stan_LPD_post <- rstan::extract(model_fit, pars = "LOO_LPD_post")[[1]] %>% as.numeric()

    stan_vs_R_LPD <- tibble(H_hat = H_hat) %>%
      # Take observed pre- and post-vaccination titer
      mutate(Ht = holdout_info$LOO_Ht,
             L = L, 
             U = U) %>%
      # Take posterior sample of measurement error sds (sigma_FRNT or sigma_HAI)
      mutate(sigma = sigma) %>%
      rowwise() %>%
      # Compute density using R function
      mutate(R_function_LPD_post = log_Pr_H(Ht = Ht, H_hat = H_hat, sigma = sigma, L = L, U = U)) %>%
      ungroup() %>%
      # Compare with Stan output
      mutate(stan_LPD_post = stan_LPD_post)  %>%
      mutate(diff_post = abs(R_function_LPD_post - stan_LPD_post))

    
    if(any(is.infinite(stan_vs_R_LPD$stan_LPD_post)) || any(stan_vs_R_LPD$diff_post > 1e-3)){
        warning(paste0("Underflow or discrepancy in LPD. Check the output file for details."))
        write_csv(
          stan_vs_R_LPD %>%
            filter(is.infinite(stan_vs_R_LPD$stan_LPD_post) | diff_post > 1e-3),
          file = file.path(LOO_dir_path, paste0("underflow_", holdout_index, ".csv")))
    }

  }

  # Write LOO_results
  dir.create(LOO_dir_path, recursive = TRUE, showWarnings = FALSE)
  write_csv(LOO_results, LOO_out_path) 

  # Remove temporary file that prevents other jobs for the same specs/holdout unit
  print(">>>>> Removing temporary lock file >>>>>")
  if (file.exists(LOO_tmp_path)) {
    file.remove(LOO_tmp_path)
  }
}
