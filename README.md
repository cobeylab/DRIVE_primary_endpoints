# DRIVE

Code for Vieira et al, Repeated annual vaccination attenuates influenza vaccine responses in a randomized placebo-controlled trial.


## Dependencies ##


## 1. Analysis of the clinical trial endpoints (Fig. 1)

1.1. Run `infer_infections_via_NAI.R` to export lists of potential infections detected via rises in NAI titers. The next scripts assume these lists have been produced. The results will be written to `results/NAI_inferred_infections`.

1.2. Run `endpoint_analyses.R` to peform the analysis of primary and some secondary endpoints of the clinical trial, shown in Fig. 1 and S2-S14. The results will be written to `results/endpoint_analyses/`.

## 2. Bayesian hierarchical model of influenza vaccine responses

2.1. Run `process_data_for_bayesian_model.R` to prepare the input data for inference. This will create a `processed_data.csv` file in `results/bayesian_fits_real_data/` and a `data_scaffold.csv` in `results/synthetic_data_experiments/`. The latter is used as a basis to generate synthetic data for model validation.

2.2.


