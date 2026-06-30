# DRIVE

Code for Vieira et al, Repeated annual vaccination attenuates influenza vaccine responses in a randomized placebo-controlled trial.


## Dependencies ##


## 1. Analysis of the clinical trial endpoints (Fig. 1)

1.1. Run `infer_infections_via_NAI.R` to export lists of potential infections detected via rises in NAI titers. The next scripts assume these lists have been produced. The results will be written to `results/NAI_inferred_infections`.

1.2. Run `endpoint_analyses.R` to peform the analysis of primary and some secondary endpoints of the clinical trial, shown in Fig. 1D and S2-S14. The results will be written to `results/endpoint_analyses/`.

## 2. Bayesian hierarchical model of influenza vaccine responses

2.1. Run `process_data_for_bayesian_model.R` to prepare the input data for inference. This will create a `processed_data.csv` file in `results/bayesian_fits_real_data/` and a `data_scaffold.csv` in `results/synthetic_data_experiments/`. The latter is used as a basis to generate synthetic data for model validation.

2.2. `run_synth_data_experiment.sh` is the master script for running synthetic data experiments on a SLURM-based computing cluster. Given an experiment directory, it calls `fit_bayesian_model_synthetic_data.R` to run the desired number of replicate experiments, each consisting of an independent realization of a synthetic dataset followed by model fitting via Stan. Each experiment directory in `results/synthetic_data_experiments/` contains a `setup.R` file with the configurations for that particular experiment. 

2.3 `combine_synth_experiment_replicates.sh` combines results from independent realizations of a synthetic data experiment into a single `combined_results.csv` file exported to the experiment directory.

2.4 `process_synth_data_experiments.R` looks in `results/synthetic_data_experiments/` for experiment directories that have `combined_results.csv` files and, for each experiment, plots summaries of the model inference across replicate synthetic datasets.







