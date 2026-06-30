library(tidyverse)
library(rstan)
library(cowplot)
theme_set(theme_cowplot())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
source("R/base_functions.R")
source("R/bayesian_model_functions.R")

start_time <- Sys.time()

# Read held-out unit from args 
args <- commandArgs(trailingOnly = TRUE)
# (if none passed, defaults to 0; no LOO-CV performed)
holdout_index <- if (length(args) > 0) as.integer(args[1]) else 0
if (is.na(holdout_index) || holdout_index < 0) {
  stop("holdout_index must be an integer >= 0")
}

# By default, any cross-validation will be done with individual/strain/year
# as unit type (1 = individual/strain/year, 2 = individual)
loglik_unit <-  1 

null_model_dir <- "results/null_model_fits_real_data/null_specific_model/no_individual_effects/"
dir.create(null_model_dir, recursive = TRUE, showWarnings = FALSE)

# Early exit if LOO output already exists (avoid running Stan)
if (holdout_index > 0) {
  LOO_subdir <- if (loglik_unit == 1) "LOO_Hpre_group" else "LOO_individual"
  LOO_dir_path <- file.path(null_model_dir, LOO_subdir)
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


# Read processed data
input_data <- read_csv("results/bayesian_fits_real_data/processed_data.csv")

# Downstream functions require a specification of RVE levels. 
# We'll use an arbitrary model, but this won't matter for the null model
RVE_model <- fully_shared_model

input_data <- input_data %>% assign_RVE_levels(RVE_model = RVE_model)


# Fitting to day 30 only)
input_data <- input_data %>%
    filter(timepoint == 30)


input_data <- input_data %>%
  mutate(strain = factor(strain, levels = unique(input_data$strain)))

input_list_real <- prepare_stan_input(
  data = input_data,
  censoring = 1, 
  include_kinetics = F,
  include_RVE = T, # Irrelevant for null model but a value required
  include_upeak_shared = F,
  include_upeak_subtype = F,
  include_upeak_year = F,
  include_k_peak_cov2 = F,
  RVE_model = RVE_model,
  loglik_unit = loglik_unit,
  holdout_index = holdout_index,
  include_Hpre_priors = F, # Irrelevant for null model but value required,
  include_NAI_infections = F # Irrelevant but value required
)

# Skip 
skip_missing_LOO_indices()

n_iterations <- 5000
n_chains <- 4

# Require these diagnostics for the log-posterior, otherwise repeat
# not required when performing LOO (holdout_index != 0)
ess_cutoff <- if (holdout_index != 0) 0 else n_chains * 100
rhat_cutoff <- if (holdout_index != 0) Inf else 1.01

repeat {
  model_fit <- stan(file = "R/null_model.stan",
                    data = input_list_real, iter = n_iterations,
                    chains = n_chains)

  diagnostics <- get_stan_diagnostics(model_fit)

  lp_bulk_ESS <- diagnostics %>% 
    filter(par == "lp__") %>%
    pull(Bulk_ESS)

  lp_rhat <- diagnostics %>% 
    filter(par == "lp__") %>%
    pull(Rhat)

  if((lp_bulk_ESS >= ess_cutoff) && round(lp_rhat, 2) <= rhat_cutoff){
    break
  }else{
      print(paste0("log-posterior ESS: ", lp_bulk_ESS, "; Rhat: ", round(lp_rhat, 3), ". Repeating"))
  }
}

end_time <- Sys.time()

# If no LOO was performed, export parameter estimates
if(holdout_index == 0){

    parameter_estimates <- rstan::extract(
        model_fit,
        pars = c(names(model_fit)[str_detect(names(model_fit), "null_mean")],
                 names(model_fit)[str_detect(names(model_fit), "sigma")], "lp__"),
        permuted = TRUE) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
        group_by(parameter) %>%
        summarise(
          mean = mean(value),
          lower = quantile(value, 0.025),
          upper = quantile(value, 0.975)
        )
    
    parameter_estimates <- parameter_estimates %>%
      left_join(diagnostics %>% rename(parameter = par))
    
    write_csv(parameter_estimates, file = paste0(null_model_dir, "/parameter_estimates.csv"))

    # Export, for each observation, the mean prediction and the censored mean residual across posterior samples.
    # Since the mean prediction is the same for all observations in each strain, first we compute the mean prediction by strain

    mean_prediction_post <- as.data.frame(model_fit) %>%
      as_tibble() %>%
      select(matches("null_mean_Ht")) %>%
      mutate(posterior_sample_index = 1:n()) %>%
      pivot_longer(cols = !any_of("posterior_sample_index"), names_to = "parameter") %>%
      mutate(strain_index = as.integer(str_extract(parameter, "[0-9]+"))) %>%
      mutate(strain = levels(input_data$strain)[strain_index]) %>%
      select(-strain_index) %>%
      mutate(parameter = str_remove(parameter, "\\[[0-9]+\\]")) %>%
      pivot_wider(names_from = parameter, values_from = value) %>%
      group_by(strain) %>%
      summarise(mean_prediction_post = mean(null_mean_Ht)) %>%
      ungroup()

    # Use that to compute the residuals
    model_residuals <- input_data %>%
      filter(measurement_replicate == 1) %>%
      # For consistency with LOO-CV, we take the first measurement if the same sample was measured twice
      select(individual, year, treatment, subtype, strain, titer_type, undetectable_value, Ht) %>%
      left_join(mean_prediction_post) %>%
      mutate(censored_mean_prediction_post = log2(censor_titers(2^mean_prediction_post, assay = titer_type)),
             # Rounding to prevent small numerical differences (though censoring already applied)
             censored_mean_residual_post = round(Ht - censored_mean_prediction_post))

    model_info <- retrieve_model_info(null_model_dir)

    model_residuals <- bind_cols(model_info, model_residuals)

    write_csv(model_residuals, paste0(null_model_dir, "/residuals.csv"))


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
  # (This doesn't affect LOO LPD calculations, just simplifies computing additional metrics of predictive ability)
  holdout_info <- holdout_info %>%
    mutate(has_replicate_measurement = any(measurement_replicate > 1)) %>%
    filter(measurement_replicate == 1) %>%
    select(-measurement_replicate, -timepoint)

  stopifnot(nrow(holdout_info) == 1)

  names(holdout_info) <- paste0("LOO_", names(holdout_info))
  
  titer_type = unique(input_data$titer_type[input_list_real$Hpre_group_index == holdout_index])

  LOO_results <- holdout_info %>%
    cross_join(
        rstan::extract(model_fit, pars = "LOO_LPD_post") %>%
          as_tibble()  %>%
          summarise(LOO_LPD_post_mean = mean(LOO_LPD_post),
                    LOO_LPD_post_lower = quantile(LOO_LPD_post, 0.025),
                    LOO_LPD_post_upper = quantile(LOO_LPD_post, 0.975)) %>%
          mutate(lp_Rhat = diagnostics %>% filter(par == "lp__") %>% pull(Rhat),
                 valid = diagnostics %>% filter(par == "lp__") %>% pull(valid),
                 Bulk_ESS = diagnostics %>% filter(par == "lp__") %>% pull(Bulk_ESS),
                 Tail_ESS = diagnostics %>% filter(par == "lp__") %>% pull(Tail_ESS),
                 start_time = start_time, 
                 end_time = end_time)

    )

    # Other metrics of predictive performance
    strain_index <- input_list_real$strain[input_list_real$Hpre_group_index == holdout_index][1]
    mean_par_names <- paste0("null_mean_", c("Ht", "Hpre"), "[", strain_index, "]")

    predictive_performance <- 
      rstan::extract(model_fit, pars = mean_par_names) %>%
          as_tibble() %>%
          rename_with(~ str_remove(.x, "\\[.*\\]$")) %>%
          mutate(LOO_Ht = holdout_info$LOO_Ht,
                 LOO_Hpre = holdout_info$LOO_Hpre) %>%
          mutate(titer_type = titer_type,
                 censored_null_mean_Ht = log2(censor_titers(titer = 2^null_mean_Ht,
                                              assay = titer_type)),
                 censored_null_mean_Hpre = log2(censor_titers(titer = 2^null_mean_Hpre,
                                                assay = titer_type))
          ) %>%
          summarise(
            mean_residual_post = mean(LOO_Ht - null_mean_Ht),
            mean_censored_residual_post = mean(LOO_Ht - censored_null_mean_Ht),
            mean_prediction_post = mean(null_mean_Ht),
            mean_censored_prediction_post = mean(censored_null_mean_Ht),
            prob_LOO_postvax_equal_prediction = mean(abs(censored_null_mean_Ht - LOO_Ht) < 1e-6),
            prob_LOO_postvax_within_1_unit = mean(abs(censored_null_mean_Ht - LOO_Ht) <= 1 + 1e-6),
            mean_residual_pre = mean(LOO_Hpre - null_mean_Hpre),
            mean_censored_residual_pre = mean(LOO_Hpre - censored_null_mean_Hpre),
            mean_prediction_pre = mean(null_mean_Hpre),
            mean_censored_prediction_pre = mean(censored_null_mean_Hpre),
            prob_LOO_prevax_equal_prediction = mean(abs(censored_null_mean_Hpre - LOO_Hpre) < 1e-6),
            prob_LOO_prevax_within_1_unit = mean(abs(censored_null_mean_Hpre - LOO_Hpre) <= 1 + 1e-6)
          )

  LOO_results <- LOO_results %>%
    bind_cols(predictive_performance) %>%
    select(matches("LOO"), matches("mean_"), everything())
    
  # Write LOO_results
  dir.create(LOO_dir_path, recursive = TRUE, showWarnings = FALSE)
  write_csv(LOO_results, LOO_out_path)
  
  # Remove temporary file that prevents other jobs for the same specs/holdout unit
  print(">>>>> Removing temporary lock file >>>>>")
  if (file.exists(LOO_tmp_path)) {
    file.remove(LOO_tmp_path)
  }
  
}