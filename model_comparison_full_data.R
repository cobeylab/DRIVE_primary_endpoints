library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
source("R/bayesian_model_functions.R")

dir.create("results/model_comparison/", showWarnings = F)

# Recursively find all directories with a model_fit.rds file
fit_dirs <- list.files(
    path = "results/bayesian_fits_real_data",
    pattern = "residuals.csv",
    recursive = TRUE,
    full.names = TRUE) %>%
    dirname() %>%
    unique()

fit_dirs <- c(fit_dirs, "results/null_model_fits_real_data/null_specific_model/no_individual_effects/")

# Models with an effect of time since vaccination were fitted as sensitivity analysis, not included in these comparisons 
fit_dirs <- fit_dirs[!str_detect(fit_dirs, "TIME_SINCE_VAX")] 

# Read residuals.csv files in these directories, bind rows
model_residuals <- map_dfr(fit_dirs, ~ read_csv(file.path(.x, "residuals.csv")))

model_residuals <- model_residuals %>%
    assign_model_labels() %>%
    mutate(subtype = factor(subtype, levels = subtype_levels))

# For the main comparison, exclude models fitted with NAI-inferred infections
main_comparison_tibble <- model_residuals %>%
    filter(!NAI_infections, individual_effects == "none") %>%
    # and models without the full set of "default" covariates (except for the null model)
    filter(str_detect(model, "null") | (include_k_peak_flu & include_age & include_sex &
            include_BMI & include_smoking & include_asthma))

residual_summary <- summarize_model_residuals(model_residuals = main_comparison_tibble, by_subtype = F)
residual_summary_by_subtype <- summarize_model_residuals(model_residuals = main_comparison_tibble, by_subtype = T)


summary_plots <- plot_fixed_effects_comparison(main_comparison_tibble, by_subtype = F)
summary_plots_by_subtype <- plot_fixed_effects_comparison(main_comparison_tibble, by_subtype = T)

summary_plots$abs_error_plot
summary_plots$percents_plot

summary_plots_by_subtype$abs_error_plot
summary_plots_by_subtype$percents_plot

    
main_comparison_tibble %>%
    group_by(model, subtype, model_type, Hpre_priors, model_label) %>%
    summarise(RMSE = sqrt(mean(censored_mean_residual_post^2))) %>%
    ggplot(aes(x = model_label, y = RMSE, color = Hpre_priors, group = Hpre_priors)) +
    geom_point(size = 4) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap("subtype")    


# Compare model fitted with/without NAI-inferred infections
NAI_infections_comparison_tibble <- model_residuals %>%
    filter(model == "fully_shared_model", Hpre_priors == F, individual_effects == "none") %>%
    # (Between models with all default covariates)
    filter(include_k_peak_flu & include_age & include_sex & include_BMI & include_smoking & include_asthma)

plot_grid(
   summarize_model_residuals(NAI_infections_comparison_tibble, by_subtype = F) %>%
    make_percent_within_2fold_pl(xvar = "NAI_infections"),
   summarize_model_residuals(NAI_infections_comparison_tibble, by_subtype = F) %>%
    make_percent_identical_pl(xvar = "NAI_infections"),
   nrow = 2
)


# Compare individual effects 
individual_effects_comparison_tibble <- model_residuals %>%
        filter(model %in% c("fully_shared_model", "subtype_specific_model"),
               !NAI_infections) %>%
        filter(include_k_peak_flu & include_age & include_sex & include_BMI & include_smoking & include_asthma)


individual_effects_comparison <- plot_individual_effects_comparison(individual_effects_comparison_tibble, by_subtype = F)

individual_effects_comparison_by_subtype <- plot_individual_effects_comparison(individual_effects_comparison_tibble, by_subtype = T)


individual_effects_comparison$abs_error_plot
individual_effects_comparison$percents_plot

individual_effects_comparison_by_subtype$abs_error_plot
individual_effects_comparison_by_subtype$percents_plot

# Compare full covariate set vs. sig. effects only vs. pre-vaccination titers + vaccination history only
covariate_set_comparison_tibble <- model_residuals %>%
    filter(model == "fully_shared_model", Hpre_priors == F, !NAI_infections, individual_effects == "none")

covariate_set_comparison <- plot_covariate_set_comparison(covariate_set_comparison_tibble, by_subtype = F)
covariate_set_comparison_by_subtype <- plot_covariate_set_comparison(covariate_set_comparison_tibble, by_subtype = T)   

covariate_set_comparison$abs_error_plot
covariate_set_comparison$percents_plot

save_plot(
    "results/model_comparison/full_data_covariate_set_percentages.pdf",
    covariate_set_comparison$percents_plot,
    base_width = supp_fig_witdh_full * 2 / 3,
    base_height = supp_fig_height_full / 2
)

ggsave(
  filename = "results/model_comparison/full_data_covariate_set_percentages.png",
  plot = covariate_set_comparison$percents_plot,
  width = supp_fig_witdh_full * 2 / 3,
  height = supp_fig_height_full / 2,
  units = "in",
  dpi = 500
)

covariate_set_comparison_by_subtype$abs_error_plot
covariate_set_comparison_by_subtype$percents_plot
