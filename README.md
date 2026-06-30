# DRIVE

Code for Vieira et al, Repeated annual vaccination attenuates influenza vaccine responses in a randomized placebo-controlled trial.


## Dependencies ##


## 1. Analysis of the clinical trial endpoints (Fig. 1)

1.1. Run `infer_infections_via_NAI.R` to export lists of potential infections detected via rises in NAI titers. The next scripts assume these lists have been produced. The results will be written to `results/NAI_inferred_infections`.

1.2. Run `endpoint_analyses.R` to peform the analysis of primary and some secondary endpoints of the clinical trial, shown in Fig. 1D and S2-S14. The results will be written to `results/endpoint_analyses/`.

## 2. Bayesian hierarchical model of influenza vaccine responses (Fig. 2)

2.1. Run `process_data_for_bayesian_model.R` to prepare the input data for inference. This will create a `processed_data.csv` file in `results/bayesian_fits_real_data/` and a `data_scaffold.csv` in `results/synthetic_data_experiments/`. The latter is used as a basis to generate synthetic data for model validation.

2.2. `run_synth_data_experiment.sh` is the master script for running synthetic data experiments on a SLURM-based computing cluster. Given an experiment directory, it calls `fit_bayesian_model_synthetic_data.R` to run the desired number of replicate experiments, each consisting of an independent realization of a synthetic dataset followed by model fitting via Stan. Each experiment directory in `results/synthetic_data_experiments/` contains a `setup.R` file with the configurations for that particular experiment. 

2.3 `combine_synth_experiment_replicates.sh` combines results from independent realizations of a synthetic data experiment into a single `combined_results.csv` file exported to the experiment directory.

2.4 `process_synth_data_experiments.R` looks in `results/synthetic_data_experiments/` for experiment directories that have `combined_results.csv` files and, for each experiment, plots summaries of the model inference across replicate synthetic datasets.

2.5 `fit_bayesian_model_real_data.R` fits the Bayesian hierarchical model to the real (as opposed to synthetic) data, given a directory in `results/bayesian_fits_real_data` containing a specification file `fit_specs.R`. The specification file controls which version of the model is fitted, what covariates and individual effects are are included, and whether informative priors are put on the population distribution of latent pre-vaccination titers. `run_bayesian_fits_real_data.sh` calls `fit_bayesian_model_real_data.R` to perform fits to the entire data set in a SLURM-based cluster. `fit_null_model_real_data.R` fits the null model (in which the predicted post-vaccination titer for all participants is simply a strain-specific population average; no specification file needed).

2.6 `process_bayesian_model_fits.R` summarizes the inference results for each model configuration fitted to the observed data, producing plots like those shown in Fig 2B-D and Figs. S17, S20 and S24.

2.7 `run_LOO.sh` performs leave-one-out cross-validation on a SLURM-based cluster given a directory with a fit specification file, automatically calling either `fit_bayesian_model_real_data.R` or `fit_null_model_real_data.R` as appropriate. `combine_LOO_output_files.sh` combines the results across held-out units to produce a single `combined_results.csv` file with cross-validation results for a given model configuration.

2.8 `model_comparison_full_data.R` and `model_comparison_LOO.R` plot performance comparisons for models fitted to the full data (Fig. S19) or in leave-one-out cross-validation (Figs. S15, S16, S18 and S21).







