library(tidyverse)
library(cowplot)
library(ggtext)
theme_set(theme_cowplot())
source("R/base_functions.R")

add_param_type <- function(df) {
  annotated_df <- df %>%
        mutate(param_type = case_when(
            str_detect(parameter, "r\\[") ~ "Prior vaccination effects",
            str_detect(parameter, "b_peak") ~ "Strain-specific antibody\nceiling effect",
            str_detect(parameter, "sigma") ~ "Standard deviations of measurement\nerror and individual effects",
            # We coded infection effects with k instead of beta, but showing together
            # with other covariates in the manuscript.
            str_detect(parameter, "k_peak") ~ "Covariate effects", 
            str_detect(parameter, "a\\[") ~ "Strain-specific\nintercept",
            str_detect(parameter, "beta_") ~ "Covariate effects",
            T ~ parameter
        ))

    param_type_values <- unique(annotated_df$param_type)

    param_type_levels <- param_type_values[c(
        which(str_detect(param_type_values, "intercept")),
        which(str_detect(param_type_values, "ceiling")),
        which(str_detect(param_type_values, "vaccination")),
        which(str_detect(param_type_values, "Covariate")),
        which(str_detect(param_type_values, "error"))
    )]

    param_type_levels <- c(param_type_levels,
        setdiff(param_type_values, param_type_levels))

    annotated_df <- annotated_df %>%
        mutate(param_type = factor(param_type, levels = param_type_levels))

    return(annotated_df)
}

param_labeller <- function(parameters){
    tibble(parameters) %>%
        mutate(param_label = case_when(
            str_detect(parameters, "a\\[") ~ "",
            str_detect(parameters, "r\\[") ~ "",
            str_detect(parameters, "b_peak") ~ "",
            str_detect(parameters, "age3140") ~ "Age 31-40",
            str_detect(parameters, "age4150") ~ "Age 41-50",
            str_detect(parameters, "asthma") ~ "Asthma",
            str_detect(parameters, "BMI_over") ~ "BMI ≥ 25",
            str_detect(parameters, "BMI_under") ~ "BMI < 18.5",
            str_detect(parameters, "beta_M") ~ "Male",            
            str_detect(parameters, "sigma_FRNT") ~ "σ<sub>FRNT</sub>",
            str_detect(parameters, "sigma_HAI") ~ "σ<sub>HAI</sub>",
            str_detect(parameters, "sigma_upeak_shared") ~ "σ<sub>global</sub>",            
            str_detect(parameters, "sigma_upeak_year") ~ "σ<sub>year</sub>",             
            str_detect(parameters, "sigma_upeak_subtype") ~ "σ<sub>subtype</sub>",             
            str_detect(parameters, "beta_smoking") ~ "Smoking",            
            str_detect(parameters, "k_peak_flu") ~ "Subtype-matched<br>influenza infection",        
            T ~ parameters
        )) %>%
        pull(param_label)
}

# Get all experiment directories within synthetic_data_experiments/ that contain combined_results.csv
base_dir <- "synthetic_data_experiments"
experiment_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
# Experiment directories are one level deeper (i.e., subdirectories of subdirectories of base_dir)
experiment_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE) %>%
    map(list.dirs, recursive = FALSE, full.names = TRUE) %>%
    unlist()

experiment_dirs <- experiment_dirs[file.exists(file.path(experiment_dirs, "combined_results.csv"))]

for (experiment_dir in experiment_dirs) {
    # Create plots directory if it doesn't exist
    plots_dir <- file.path(experiment_dir, "plots")
    if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

    inference_results <- read_csv(file.path(experiment_dir, "combined_results.csv"))

    n_pars_per_replicate <- inference_results %>%
        group_by(replicate_id) %>%
        summarise(S = length(unique(parameter))) %>%
        pull(S) %>%
        unique()

    if(length(n_pars_per_replicate) > 1){
        stop(paste("Not all replicates have the same number of parameters in", experiment_dir))
    }

    if(any(duplicated(inference_results %>% select(replicate_id, parameter)))) {
        stop(paste("Some combinations of replicate_id and parameter occur more than once in", experiment_dir))
    }

    inference_results <- inference_results %>%
        mutate(error = mean - true_value)

    # Plot 1: Histogram of log-posterior estimate across replicates
    p1 <- inference_results %>%
        filter(str_detect(parameter, "lp__")) %>%
        ggplot(aes(x = mean)) +
        geom_histogram() +
        labs(x = "Estimate (mean)", y = "Count", title = "Distribution of the log-posterior across replicates")
    ggsave(file.path(plots_dir, "lp_mean_hist.pdf"), p1, width = 7, height = 5)

    # Plot 2: Histogram of bulk_ess for "__lp" across replicates
    p2 <- inference_results %>%
        filter(str_detect(parameter, "lp__")) %>%
        ggplot(aes(x = Bulk_ESS)) +
        geom_histogram() +
        labs(x = "Bulk ESS", y = "Count", title = "Distribution of lp__ Bulk ESS Across Replicates")
    ggsave(file.path(plots_dir, "lp_bulk_ess_hist.pdf"), p2, width = 7, height = 5)

    # Plot 3: Parameter estimates vs. lp__ Bulk ESS
    lp_ess <- inference_results %>%
        filter(parameter == "lp__") %>%
        select(replicate_id, Bulk_ESS) %>%
        rename(lp_bulk_ess = Bulk_ESS)

    p3 <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               parameter != "lp__",
               !str_detect(parameter, "LOO_LPD")) %>%
        left_join(lp_ess, by = "replicate_id") %>%
        ggplot(aes(x = lp_bulk_ess, y = mean)) +
        geom_point(alpha = 0.5) +
        facet_wrap(~parameter, scales = "free") +
        geom_smooth(method = 'lm') +
        labs(x = "Bulk ESS of lp__", y = "Parameter Estimate (mean)", title = "Parameter estimates vs. lp__ Bulk ESS")
    ggsave(file.path(plots_dir, "param_vs_lp_bulk_ess.pdf"), p3, width = 10, height = 7)

    # Plot 4: Histogram of errors by parameter
    p4 <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        ggplot(aes(x = error)) +
        geom_histogram() +
        facet_wrap('parameter', scales = 'free') +
        geom_vline(xintercept = 0)
    ggsave(file.path(plots_dir, "error_hist_by_param.pdf"), p4, width = 10, height = 7)

    # Plot 4b: Histogram of Bulk_ESS by parameter
    p4b <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        ggplot(aes(x = Bulk_ESS)) +
        geom_histogram() +
        facet_wrap('parameter', scales = 'free') +
        labs(x = "Bulk ESS", y = "Count", title = "Distribution of Bulk ESS by Parameter") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.margin = margin(10, 10, 40, 10))
    ggsave(file.path(plots_dir, "bulk_ess_hist_by_param.pdf"), p4b, width = 14, height = 9)

    # Plot 4c: Histogram of Tail_ESS by parameter
    p4c <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        ggplot(aes(x = Tail_ESS)) +
        geom_histogram() +
        facet_wrap('parameter', scales = 'free') +
        labs(x = "Tail ESS", y = "Count", title = "Distribution of Tail ESS by Parameter") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.margin = margin(10, 10, 40, 10))
    ggsave(file.path(plots_dir, "tail_ess_hist_by_param.pdf"), p4c, width = 14, height = 9)

    # Plot 4d: Histogram of Rhat by parameter
    p4d <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        ggplot(aes(x = Rhat)) +
        geom_histogram() +
        geom_vline(xintercept = 1.01, linetype = "dashed", color = "red") +
        facet_wrap('parameter', scales = 'free') +
        labs(x = "Rhat", y = "Count", title = "Distribution of Rhat by Parameter") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.margin = margin(10, 10, 40, 10))
    ggsave(file.path(plots_dir, "rhat_hist_by_param.pdf"), p4d, width = 14, height = 9)


    # Plot 5: Boxplot of estimates and true values
    p5 <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        add_param_type() %>%
        ggplot(aes(x = parameter)) +
        geom_point(aes(y = mean), color = "black",,
                    position = position_jitter(width = 0.2, height = 0),
                    alpha = 0.2, size = 1, shape = 1) +
        geom_boxplot(aes(y = mean), fill = "#DDDDDD", alpha = 0.5, outlier.shape = NA,
                     width = 0.5, box.linewidth = boxplot_line_width) +
        geom_point(aes(y = true_value), color = "red", size = 2) +
        facet_wrap(~param_type, scales = "free", nrow = 3) +
        baseline_figure_settings +
        theme(axis.text.x = element_markdown(size = default_figure_font_size, angle = 45, hjust = 1)) +
        labs(y = "Parameter value", x = "") +
        scale_x_discrete(labels = param_labeller)


    ggsave(file.path(plots_dir, "boxplot_estimates_true.svg"),
            p5, width = supp_fig_witdh_full,
            height = supp_fig_height_full - 2)

    # Plot 6: Coverage barplot
    p6 <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        add_param_type() %>%
        mutate(covered = if_else(true_value >= lower & true_value <= upper, "Covered", "Not Covered")) %>%
        group_by(parameter, param_type, covered) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(parameter, param_type) %>%
        mutate(percent = 100 * n / sum(n)) %>%
        ggplot(aes(x = parameter, y = percent, fill = covered)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 5, linetype = "dashed", color = "black") +
        facet_wrap(~param_type, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
        labs(y = "% Replicates", x = "Parameter", fill = "True value covered by CrI")
    ggsave(file.path(plots_dir, "coverage_barplot.pdf"), p6, width = 10, height = 9)

    # Plot 7: Error boxplot
    p7 <- inference_results %>%
        filter(!str_detect(parameter, "Hpre_"),
               !str_detect(parameter, "log_lik"),
               !str_detect(parameter, "lp__"),
               !str_detect(parameter, "LOO_LPD")) %>%
        add_param_type() %>%
        ggplot(aes(x = parameter, y = error)) +
        geom_point(color = "black", alpha = 0.2, size = 2, shape = 1,
                   position = position_jitter(width = 0.2, height = 0)) +
        geom_boxplot(fill = "#DDDDDD", alpha = 0.5, outlier.shape = NA, width = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        facet_wrap(~param_type, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Error (mean - true)", x = "Parameter")
    ggsave(file.path(plots_dir, "error_boxplot.pdf"), p7, width = 15, height = 10)
}
