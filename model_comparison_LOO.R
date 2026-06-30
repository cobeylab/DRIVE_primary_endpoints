library(tidyverse)
library(cowplot)
library(ggpubr)
library(ggh4x)
theme_set(theme_cowplot())
source("R/base_functions.R")
source("R/bayesian_model_functions.R")

fit_directories <- list.dirs(getwd(), recursive = TRUE, full.names = TRUE) %>%
    keep(~ basename(.) == "LOO_Hpre_group")

fit_directories <- fit_directories %>%
    keep(~ file.exists(file.path(., "combined_results.csv")))

# Ignore some older exploratory fits
fit_directories <- fit_directories[!str_detect(fit_directories, "OLDER")]
# And the sensitivity analysis that includes an effect of time since prior vaccination.
fit_directories <- fit_directories[!str_detect(fit_directories, "TIME_SINCE_VAX")]

combined_results <- lapply(as.list(fit_directories),
    function(dir){
        results <- retrieve_model_info(dirname(dir)) %>%
            bind_cols(as_tibble(read.csv(paste0(dir, "/combined_results.csv"))))

            names(results)[names(results) == "Bulk_ESS"] <- "lp_Bulk_ESS"
            names(results)[names(results) == "Tail_ESS"] <- "lp_Tail_ESS"
        return(results)
    }
) %>%
    bind_rows() %>%
    mutate(
        start_time = as.POSIXct(start_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        end_time = as.POSIXct(end_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        duration_hours = as.numeric(difftime(end_time, start_time, units = "hours"))
    ) %>%
    mutate(subtype = subtype_labeller(LOO_strain),
           subtype = factor(subtype, levels = subtype_levels),
           n_previous_vax = count_prior_vax(LOO_treatment))

# The LOO-CV code computes among other things the mean value of the censored predictions for each unit across posterior samples
# A more interesting quantity is the censored value of the mean of uncensored predictions. We compute that here
# (I.e., we compute the censored mean prediction in addition to the mean censored prediction)
combined_results <- combined_results %>%
    mutate(assay = ifelse(subtype == "H3N2", "FRNT", "HAI"),
           censored_mean_prediction_post = log2(censor_titers(2^mean_prediction_post, assay = assay)),
           censored_mean_prediction_pre = log2(censor_titers(2^mean_prediction_pre, assay = assay)))

# Similarly, we define the censored mean residual, based on the censored mean prediction
# (as opposed to the "mean censored residual" based on the mean censored prediction)
combined_results <- combined_results %>%
    # Rounding to prevent small numerical differences (although censoring already applied)
    mutate(censored_mean_residual_post = round(LOO_Ht - censored_mean_prediction_post),
           censored_mean_residual_pre = round(LOO_Hpre - censored_mean_prediction_pre))

    
combined_results <- combined_results %>%
    assign_model_labels()

# How long the job arrays take
group_vars <- c("model", "individual_effects", "Hpre_priors", "NAI_infections",
                combined_results %>%
                    select(matches("include")) %>%
                    names())

combined_results %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(n_jobs = n(),
              total_duration = sum(duration_hours),
              max_duration = max(duration_hours),
              percentile_95 = quantile(duration_hours, 0.95, na.rm = T),
              array_timespan = as.numeric(difftime(max(end_time), min(start_time), units = "hours")))

#  ========== Comparing fixed effects models  ===========

fixed_effects_comparison_tibble <- combined_results  %>%
    filter(individual_effects == "none",
           # The default comparisons do not include NAI-inferred infections...
           !NAI_infections,
           # and consider only models with all these covariates included
           # (plus the null model)
           str_detect(model, "null") |
             (include_sex & include_age & include_BMI & include_smoking & include_asthma & include_k_peak_flu)
           ) 
    
# Summary plots showing absolute error and proportion of observations <= 2-fold or identical to prediction
fixed_effects_comparison <- plot_fixed_effects_comparison(fixed_effects_comparison_tibble, by_subtype = F)

fixed_effects_comparison_by_subtype <- plot_fixed_effects_comparison(fixed_effects_comparison_tibble, by_subtype = T)

save_plot(
    "results/model_comparison/LOO_CV_fixed_effects_abs_error.pdf",
    fixed_effects_comparison$abs_error_plot +
        theme(axis.text.x = element_text(size = small_font_size)),
    base_height = supp_fig_height_full / 2.5,
    base_width = supp_fig_witdh_full
)

save_plot(
    "results/model_comparison/LOO_CV_fixed_effects_percentages.pdf",
    fixed_effects_comparison$percents_plot,
    base_height = 5,
    base_width = supp_fig_witdh_full
)


#  ========== Comparing individual effects on a restricted set of fixed-effects forms  ===========
individual_effects_comparison_tibble <- combined_results %>%
        filter(model %in% c("fully_shared_model", "subtype_specific_model"),
               !NAI_infections,
               str_detect(model, "null") |
                 (include_sex & include_age & include_BMI & include_smoking & include_asthma & include_k_peak_flu)
            )


individual_effects_comparison <- plot_individual_effects_comparison(individual_effects_comparison_tibble, by_subtype = F)

individual_effects_comparison_by_subtype <- plot_individual_effects_comparison(individual_effects_comparison_tibble, by_subtype = T)


individual_effects_comparison$abs_error_plot
individual_effects_comparison$percents_plot

save_plot(
    "results/model_comparison/LOO_CV_individual_effects_percentages.pdf",
    individual_effects_comparison$percents_plot,
    base_height = supp_fig_height_full / 1.5,
    base_width = supp_fig_witdh_full
)

individual_effects_comparison_by_subtype$abs_error_plot
individual_effects_comparison_by_subtype$percents_plot


# ========= For best-performing model determined after fixed-effects and individual-effects selection ======
# =========== compare pre-vax titers + vax history alone vs . all default covariates

covariate_set_comparison_tibble <- combined_results %>%
    filter(model == "fully_shared_model", Hpre_priors == F, !NAI_infections, individual_effects == "none",
           # For this comparison, exclude a model with prior infection effects removed but other default
           # covariates kept
           !(include_age & !include_k_peak_flu))

covariate_set_comparison <- plot_covariate_set_comparison(covariate_set_comparison_tibble, by_subtype = F)
covariate_set_comparison_by_subtype <- plot_covariate_set_comparison(covariate_set_comparison_tibble, by_subtype = T)   

covariate_set_comparison$abs_error_plot
covariate_set_comparison$percents_plot

save_plot(
    "results/model_comparison/LOO_CV_covariate_set_percentages.pdf",
    covariate_set_comparison$percents_plot,
    base_width = supp_fig_witdh_full / 2,
    base_height = supp_fig_height_full / 2
)

ggsave(
  filename = "results/model_comparison/LOO_CV_covariate_set_percentages.png",
  plot = covariate_set_comparison$percents_plot,
  width = supp_fig_witdh_full / 2,
  height = supp_fig_height_full / 2,
  units = "in",
  dpi = 500
)