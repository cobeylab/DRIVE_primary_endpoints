library(tidyverse)
library(cowplot)
library(rstan)
n_cores <- parallelly::availableCores()
options(mc.cores = n_cores)
#rstan_options(auto_write = TRUE)

source("R/base_functions.R")
source("R/bayesian_model_functions.R")

args <- commandArgs(trailingOnly = TRUE)
#e.g. args <- c("results/synthetic_data_experiments/fully_flexible_RVE_model/default/", F)

experiment_dir <- args[1] # Directory where results will be exported

# By default, look for a file called 'setup.R' inside experiment_dir
setup_R_file <- file.path(experiment_dir, "setup.R")
source(setup_R_file)

# Current version of the stan model assumes same sigma for measurement error in pre- and post-vaccination titers
# include_prevax_error_synth_data must be TRUE or the Stan model will be mispecified relative to the synthetic data.
stopifnot(include_prevax_error_synth_data) 


detailed_replicate <- as.logical(args[2]) # If TRUE, exports detailed results (will be set to F when running replicate experiments)

if (length(args) < 2 || is.na(detailed_replicate)) {
  detailed_replicate <- FALSE
}

# Define output directory based on detailed_replicate
if (detailed_replicate) {
  output_dir <- file.path(experiment_dir, "detailed_replicate")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
} else {
  output_dir <- experiment_dir
}

# Use Sys.time to generate an id for this replicate of the experiment
replicate_id <- paste0(Sys.time(), "_", Sys.getpid()) %>%
  str_replace_all(" ", "_") %>%
  str_replace_all(":", "-")

summary_output_file <- file.path(output_dir, paste0("inference_results_", replicate_id, ".csv"))
synth_data_file <- file.path(experiment_dir, "detailed_replicate", "synth_data_object.RData")


# If fixed_synth_data is TRUE and detailed_replicate is F, use pre-existing fixed synthetic data instead of 
# simulating a new realization
if (fixed_synth_data && !detailed_replicate) {
  if (!file.exists(synth_data_file)) {
    warning("Experiment uses fixed synthetic data but no pre-existing fixed synthetic data was found. Run experiment with detailed_replicate set to TRUE to create one.")
    quit(status = 1)
  }
  load(synth_data_file) # loads 'synth_data_object'
  message(strrep("%", 40))
  message("Fitting to fixed synthetic data")
  message(strrep("%", 40))
} else {
  # Otherwise, simulate new realization of synthetic data...
  synth_data_object <- simulate_vaccine_responses(data_scaffold,
                                                  RVE_model = RVE_model,
                                                  params = params,
                                                  censoring = censoring,
                                                  include_kinetics = include_kinetics_synth_data,
                                                  include_prevax_measurement_error = include_prevax_error_synth_data,
                                                  include_NAI_infections = include_NAI_infections_synth_data,
                                                  scaffold_multiples = scaffold_multiples)

  # ...and export the data if detailed_replicate is TRUE
  # (This data will then be used when fixed_synth_data && !detailed_replicate)
  if (detailed_replicate) {
    # (Saving as an RData file to preserve factor levels)
    if (file.exists(synth_data_file)) {
      warning(paste("File", synth_data_file, "already exists. Not overwriting. Quitting."))
      quit(status = 1)
    } else {
      save(synth_data_object, file = synth_data_file)
    }
  }
}

simulated_data <- synth_data_object$data

# Visualize predicted pre vs. post-vaccination titers
synth_data_object$peak_plot

# Checks RVE levels were ordered correctly in the synthetic data
stopifnot(levels(simulated_data$RVE_level_Y1) == get_RVE_model_levels(RVE_model))
stopifnot(levels(simulated_data$RVE_level_Y2) == get_RVE_model_levels(RVE_model))
stopifnot(levels(simulated_data$RVE_level_Y3) == get_RVE_model_levels(RVE_model))

# Prepare data for Stan
input_list_synthetic <- prepare_stan_input(
  data = simulated_data,
  censoring = censoring,
  include_kinetics = include_kinetics_synth_data,
  include_RVE = include_RVE,
  include_upeak_shared = include_upeak_shared,
  include_upeak_subtype = include_upeak_subtype,
  include_upeak_year = include_upeak_year,
  include_k_peak_flu = include_k_peak_flu_synth_data,
  include_k_peak_cov2 = include_k_peak_cov2,
  RVE_model = RVE_model,
  include_Hpre_priors = F, # Disabled by default in synthetic data (assumptions for undetectable pops not implemented in the simulation)
  include_NAI_infections = include_NAI_infections_synth_data,
  include_time_since_vax = include_time_since_vax_synth_data,
  include_sex = include_sex_synth_data,
  include_age = include_age_synth_data,
  include_BMI = include_BMI_synth_data,
  include_smoking = include_smoking_synth_data,
  include_asthma = include_asthma_synth_data
)

# Check level assignment was done correctly for titer type
stopifnot(all(input_list_synthetic$titer_type[str_detect(simulated_data$titer_type, "HAI")] == 1))
stopifnot(all(input_list_synthetic$titer_type[str_detect(simulated_data$titer_type, "FRNT")] == 2))

message(strrep("%", 40))
message(glue::glue("USING {n_cores} CORES"))
message(strrep("%", 40))

# Fit model using Stan and measure duration
stan_start_time <- Sys.time()
fit_synthetic <- stan(
  file = "R/linear_model.stan",
  data = input_list_synthetic,
  iter = synth_experiment_iterations,
  chains = synth_experiment_chains,
  control = list(adapt_delta = 0.9)
)
stan_end_time <- Sys.time()
stan_fit_duration <- as.numeric(difftime(stan_end_time, stan_start_time, units = "secs"))

# Process Stan output
model_results_synthetic <- process_stan_fit(fit = fit_synthetic, data = simulated_data,
                                            RVE_model = RVE_model,
                                            # No need for detailed samples of predicted quantities
                                            predicted_values_sample_size = 0) 

# Annotate model results with true parameter values
model_results_synthetic$posterior_summary_pars <- model_results_synthetic$posterior_summary_pars %>%
  mutate(parameter_label = str_extract(parameter, "^[^\\[]+"))

true_values_df <- bind_rows(
  # Strain-specific parameters
  model_results_synthetic$posterior_summary_pars %>%
    filter(!is.na(strain)) %>%
    mutate(parameter_label = str_extract(parameter, "^[^\\[]+")) %>%
    left_join(params$strain_specific %>%
                select(any_of(names(strain_specific_params))) %>%
                unique() %>%
                pivot_longer(cols = !matches('strain'),
                             names_to = "parameter_label", values_to = "true_value"),
              by =  join_by(strain, parameter_label)
            ) %>%
    select(parameter, true_value),
  # RVE parameters
  model_results_synthetic$posterior_summary_pars %>%
    filter(!is.na(RVE_level)) %>%
    mutate(parameter_label = str_extract(parameter, "^[^\\[]+")) %>%
    left_join(params$RVE %>%
                unique() %>%
                pivot_longer(cols = !matches('RVE_level'),
                             names_to = "parameter_label", values_to = "true_value"),
              by = join_by(RVE_level, parameter_label)) %>%
    select(parameter, true_value),
  # Universal parameters
  params$universal %>%
    unique() %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter", values_to = "true_value"),
  # Hpre means and sds (latent pre-vaccination titer distribution parameters)
  synth_data_object$Hpre_means_and_sds %>% select(parameter, true_value)
) 

# Merge posterior summary with true values
model_results_synthetic$posterior_summary_pars <- left_join(model_results_synthetic$posterior_summary_pars,
  true_values_df, by = "parameter")

results <- model_results_synthetic$posterior_summary_pars %>%
  select(parameter, true_value, mean, lower, upper) %>%
  left_join(model_results_synthetic$diagnostics %>%
              rename(parameter = par)) %>%
  mutate(
    replicate_id = replicate_id,
    censoring = censoring,
    include_prevax_measurement_error = include_prevax_error_synth_data,
    include_kinetics = include_kinetics_synth_data,
    scaffold_multiples = scaffold_multiples,
    n_iterations = synth_experiment_iterations,
    n_chains = synth_experiment_chains,
    n_cores = n_cores,
    stan_fit_duration_sec = stan_fit_duration
  ) %>%
  select(replicate_id, censoring, include_prevax_measurement_error,
         include_kinetics, scaffold_multiples, n_iterations, n_chains, n_cores,
         stan_fit_duration_sec, everything())

# Write results to a replicate-specific file
if (file.exists(summary_output_file)) {
  # Redundant with earlier check, but to be safe...
  warning(paste("File", summary_output_file, "already exists. Not overwriting. Quitting."))
  quit(status = 1)
} else {
  write_csv(results, summary_output_file)
}

if (detailed_replicate) {
  save(model_results_synthetic, file = file.path(output_dir, "model_results_synthetic.RData"))

  # Plot true vs estimated parameter values
  g1 <- plot_parameter_estimates(model_results_synthetic) + ylab("") +
    ggtitle(paste0("NOTE: Censoring is ", as.character(censoring),
                   "; used ", scaffold_multiples, " multiples of the scaffold"))
  ggsave(filename = file.path(output_dir, "parameter_estimates.pdf"), plot = g1, width = 8, height = 6)

  g2 <- plot_obs_vs_predicted_titers(data = simulated_data,
                                     model_results = model_results_synthetic,
                                     censor_predictions = censoring)
  ggsave(filename = file.path(output_dir, "obs_vs_predicted_titers.pdf"), plot = g2, width = 8, height = 6)

 

  g3 <- plot_Ht_vs_Hpre_predictions(data = simulated_data,
                                    t = 30,
                                    use_oficial_timepoint = T, data_time_tolerance = NULL,
                                    params_summary = format_params_for_plotting(model_results_synthetic),
                                    params_sample = format_params_for_plotting(model_results_synthetic, posterior_summary = F),
                                    plot_data = T)

  ggsave(filename = file.path(output_dir, "Ht_vs_Hpre_predictions.pdf"), plot = g3, width = 8, height = 6)

  ggsave(file.path(output_dir, "traceplot_posterior.pdf"), plot = traceplot(fit_synthetic, "lp__"), width = 8, height = 6)
  ggsave(file.path(output_dir, "traceplot_sigma.pdf"), plot = traceplot(fit_synthetic, "sigma"), width = 8, height = 6)
  if (include_upeak_shared) {
    ggsave(file.path(output_dir, "traceplot_sigma_upeak_shared.pdf"), plot = traceplot(fit_synthetic, "sigma_upeak_shared"), width = 8, height = 6)
  }
  ggsave(file.path(output_dir, "traceplot_a.pdf"), plot = traceplot(fit_synthetic, "a"), width = 8, height = 6)
  ggsave(file.path(output_dir, "traceplot_b_peak.pdf"), plot = traceplot(fit_synthetic, "b_peak"), width = 8, height = 6)
  ggsave(file.path(output_dir, "traceplot_r.pdf"), plot = traceplot(fit_synthetic, "r"), width = 8, height = 6)
  ggsave(file.path(output_dir, "traceplot_Hpre_mean.pdf"), plot = traceplot(fit_synthetic, "Hpre_mean"), width = 8, height = 6)
  ggsave(file.path(output_dir, "traceplot_Hpre_sd.pdf"), plot = traceplot(fit_synthetic, "Hpre_sd"), width = 8, height = 6) 
}
