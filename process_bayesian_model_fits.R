library(tidyverse)
library(cowplot)
library(rstan)
theme_set(theme_cowplot())
source("R/bayesian_model_functions.R")

fit_dirs <- commandArgs(trailingOnly = TRUE)[1]

# If no directory provided...
if(is.na(fit_dirs)){
    # Recursively find all directories with a model_results.rds file
    fit_dirs <- list.files(
    path = "results/bayesian_fits_real_data",
    pattern = "model_results\\.rds$",
    recursive = TRUE,
    full.names = TRUE) %>%
    dirname() %>%
    unique()
}

for (fit_specs_dir in fit_dirs) {

    # e.g., fit_specs_dir <- "bayesian_fits_real_data//update_specific_model/"

    # Create "plots" directory within fit_specs_dir if it doesn't exist
    plots_dir <- file.path(fit_specs_dir, "plots")
    if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
    }

    print(paste0("%%%%%%%%%% PROCESSING ", fit_specs_dir, " %%%%%%%%%%"))

    load(paste0(fit_specs_dir, "/model_results.rds"))
    load(paste0(fit_specs_dir, "/input_data.rds"))
    source(paste0(fit_specs_dir, "/fit_specs.R"))

    model_info <- retrieve_model_info(fit_specs_dir)

    # Export model residuals (for model_comparison_full_data.R)
    model_residuals <- combine_data_and_predictions(
        data = input_data,
        model_results =  model_results, 
        censor_predictions = F) %>%
        rename(mean_prediction_post = H_hat_mean)

    censored_predictions <- combine_data_and_predictions(
        data = input_data,
        model_results =  model_results,
        censor_predictions = T) %>%
        select(measurement_replicate, individual, timepoint, subtype, strain, year, H_hat_mean) %>%
        rename(censored_mean_prediction_post = H_hat_mean)
    
    model_residuals <- left_join(model_residuals, censored_predictions) %>%
        # For consistency with the performance metrics in LOO-CV analysis, we use the first replicate measurement in samples with 
        # replicated measurements 
        filter(measurement_replicate == 1) %>%
        select(individual, year, treatment, subtype, strain, titer_type,
               undetectable_value, Ht, mean_prediction_post, censored_mean_prediction_post) %>%
        # Rounding to prevent small numerical differences (though censoring already applied)
        mutate(censored_mean_residual_post = round(Ht - censored_mean_prediction_post))

    model_residuals <- bind_cols(model_info, model_residuals)

    write_csv(model_residuals, paste0(fit_specs_dir, "/residuals.csv"))
    
    # Looking at inference diagnostics
    diagnostics <- model_results$diagnostics %>%
        filter(!str_detect(par, "Hpre_"), !str_detect(par, "H_hat"), !str_detect(par, "u_peak")) %>%
        arrange(par)

    write_csv(diagnostics, paste0(plots_dir, "/inference_diagnostics.csv"))

    # Only plot RVE estimates if model included RVE
    if(model_results$posterior_summary_pars %>% filter(str_detect(parameter, "r\\[")) %>% nrow() > 0){
        RVE_estimates_pl <- plot_RVE_estimates(model_results)
        save_plot(paste0(plots_dir, "/RVE_estimates.pdf"),
            RVE_estimates_pl,
            base_height = if (str_detect(fit_specs_dir, "fully_flexible_RVE")) 5 else 2.3,
            base_width = if (str_detect(fit_specs_dir, "fully_flexible_RVE")) 18 else 2.8)

        save_plot(paste0(plots_dir, "/RVE_estimates.svg"),
            RVE_estimates_pl,
            base_height = if (str_detect(fit_specs_dir, "fully_flexible_RVE")) 5 else 2.3,
            base_width = if (str_detect(fit_specs_dir, "fully_flexible_RVE")) 18 else 2.8)

    }

    # Looking at parameter summary statistics
    posterior_summary_pars <- model_results$posterior_summary_pars %>%
        filter(!str_detect(parameter, "Hpre_"), !str_detect(parameter, "H_hat"), !str_detect(parameter, "u_peak"),
               !str_detect(parameter, "log_lik")) %>%
        arrange(parameter)

    write_csv(posterior_summary_pars, paste0(plots_dir, "/parameter_posterior_summary.csv"))


    parameter_estimates_pl <- plot_parameter_estimates(model_results, include_Hpre_means_and_sd = F)
    save_plot(paste0(plots_dir, "/parameter_estimates.pdf"),
        parameter_estimates_pl,
        base_height = 10,
        base_width = 15
    )

    Hpre_means_and_sds <- plot_Hpre_means_and_sds(model_results, input_data)
    save_plot(paste0(plots_dir, "/Hpre_means_and_sds.pdf"),
        Hpre_means_and_sds,
        base_height = 7,
        base_width = 18
    )

    obs_vs_predicted_titers <- plot_obs_vs_predicted_titers(
        data = input_data,
        model_results = model_results,
        censor_predictions = T)

    save_plot(paste0(plots_dir, "/obs_vs_predicted_titers.pdf"),
        obs_vs_predicted_titers,
        base_height = supp_fig_height_full - 3,
        base_width = supp_fig_witdh_full,
    )

    obs_vs_predicted_prevax_titers <- plot_obs_vs_predicted_prevax_titers(
        data = input_data,
        model_results = model_results,
        censor_predictions = T
    )
    
    save_plot(paste0(plots_dir, "/obs_vs_predicted_prevax_titers.pdf"),
        obs_vs_predicted_prevax_titers,
        base_height = 8,
        base_width = 9
    )

    pre_vs_post_titers <- plot_Ht_vs_Hpre_predictions(
        data = input_data,
        t = 30,
        use_oficial_timepoint = T,
        params_summary = format_params_for_plotting(model_results),
        params_sample = format_params_for_plotting(model_results, posterior_summary = F),
        plot_data = T,
        data_time_tolerance = NULL,
        include_NAI_infections = fit_specs$include_NAI_infections)

    save_plot(paste0(plots_dir, "/pre_vs_post_titers.pdf"),
        pre_vs_post_titers,
        base_height = 3.6,
        base_width = 4.1
    )

    save_plot(paste0(plots_dir, "/pre_vs_post_titers.svg"),
        pre_vs_post_titers,
        base_height = 3.6,
        base_width = 4.1
    )
    
    latent_pre_vs_obs_post_titers <- plot_Ht_vs_Hpre_predictions(
        data = input_data,
        t = 30,
        use_oficial_timepoint = T,
        params_summary = format_params_for_plotting(model_results),
        params_sample = format_params_for_plotting(model_results, posterior_summary = F),
        plot_data = T,
        data_time_tolerance = NULL,
        use_latent_prevax_titers = T,
        include_NAI_infections = fit_specs$include_NAI_infections
    )

    save_plot(paste0(plots_dir, "/latent_pre_vs_obs_post_titers.pdf"),
        latent_pre_vs_obs_post_titers,
        base_height = supp_fig_height_full / 2,
        base_width = supp_fig_witdh_full
    )

    # Figure with distribution of model predictions by timepoint

    predictions_by_timepoint_Y2 <- plot_model_predictions_by_timepoint(
        data = input_data, model_results = model_results, year = 2, censor_predictions = T, transpose = T)

    predictions_by_timepoint_Y3 <- plot_model_predictions_by_timepoint(
        data = input_data, model_results = model_results, year = 3, censor_predictions = T, transpose = T)

    predictions_by_timepoint_Y4 <- plot_model_predictions_by_timepoint(
        data = input_data, model_results = model_results, year = 4, censor_predictions = T, transpose = T)


    predictions_by_timepoint <- plot_grid(
        predictions_by_timepoint_Y2 + ggtitle("Year 2") +
            theme(plot.title = element_text(size = default_figure_font_size), legend.position = "none"),
        predictions_by_timepoint_Y3 + ggtitle("Year 3") +
            theme(plot.title = element_text(size = default_figure_font_size), legend.position = "none"),
        predictions_by_timepoint_Y4 + ggtitle("Year 4") +
            theme(plot.title = element_text(size = default_figure_font_size),, legend.position = "none",
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = small_font_size)),
        ncol = 1
    )

    save_plot(paste0(plots_dir, "/predictions_by_timepoint.pdf"),
        predictions_by_timepoint,
        base_height = supp_fig_height_full - 2,
        base_width = supp_fig_witdh_full
    )

    # Smaller version for main text figure panel
    predictions_by_timepoint_main_text <- plot_model_predictions_by_timepoint(
        data = input_data, model_results = model_results, year = 4, censor_predictions = T, subtype = c("H1N1", "B/Victoria")) +        
        theme(strip.text.x = element_blank(), legend.box.spacing = unit(0, "pt")) +
        scale_x_discrete(labels = count_prior_vax) +
        ylab("Titers 30 days after vaccination in year 4") +
        xlab("Prior vaccinations")


    save_plot(paste0(plots_dir, "/predictions_by_timepoint_main_text.pdf"),
        predictions_by_timepoint_main_text,
        base_height = supp_fig_height_full - 6,
        base_width = main_text_width_full / 3
    )

    save_plot(paste0(plots_dir, "/predictions_by_timepoint_main_text.svg"),
        predictions_by_timepoint_main_text,
        base_height = supp_fig_height_full - 6,
        base_width = main_text_width_full / 3
    )




    variance_partition <- plot_grid(
        plot_variance_partition(data = input_data, model_results = model_results,  by_individual_effect_type = F) + ggtitle("Individual effects combined into sigma2_u"),
        plot_variance_partition(data = input_data, model_results = model_results, by_individual_effect_type = T) + ggtitle("Individual effects by type"),
        nrow = 2
    )
    
    save_plot(paste0(plots_dir, "/variance_partition.pdf"),
        variance_partition,
        base_height = 7,
        base_width = 15
    )

    model_residuals_censored <- plot_model_residuals(input_data, model_results, censor_predictions = T)

    save_plot(paste0(plots_dir, "/model_residuals_censored.pdf"),
        model_residuals_censored,
        base_height = 4.5,
        base_width = 6.5
    )

    model_residuals_uncensored <- plot_model_residuals(input_data, model_results, censor_predictions = F)

    save_plot(paste0(plots_dir, "/model_residuals_uncensored.pdf"),
        model_residuals_uncensored,
        base_height = 4.5,
        base_width = 6.5
    )

    if(file.exists(paste0(fit_specs_dir, "/model_fit.rds"))){
        load(paste0(fit_specs_dir, "/model_fit.rds"))
        save_plot(
            paste0(plots_dir, "/traceplot.pdf"),
            traceplot(model_fit, diagnostics$par) +
                theme(legend.position = 'top',
                      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)),
            base_height = 10,
            base_width = 13)
    }

  # For the fully-shared RVE models, whether conditions for higher titers in 2nd vs 1st time vaccinees are met
  if(all(model_results$RVE_model$form$RVE_level == "fully_shared_RVE")){
    monotonicity_contition <- compute_nonmonotonicity_condition(model_results, input_data)
    write_csv(monotonicity_contition, paste0(plots_dir, "/non_monotonicity_condition.csv"))
  }

}
