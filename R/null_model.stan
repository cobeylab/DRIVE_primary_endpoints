functions {
  // Generalized censored likelihood function
  // Log-probability of observing a censored titer Ht given a mean of H_hat and variance sigma
  // Used for calculating the likelihood of a post-vaccination titer given a predicted value.
  real log_Pr_H(
    real Ht,
    real H_hat,
    real sigma,
    real L,
    real U
  ) {
    // Special case: Ht = L-1 → P(H_hat < L)
    if (fabs(Ht - (L - 1)) < 1e-9) {
      return normal_lcdf(L | H_hat, sigma);
    }
    // Special case: Ht = U → P(H_hat ≥ U)
    else if (fabs(Ht - U) < 1e-9) {
      return log1m_exp(normal_lcdf(U | H_hat, sigma));
    }
    // General case: P(Ht ≤ H_hat < Ht + 1)
    else {
      return log_diff_exp(
        normal_lcdf(Ht + 1 | H_hat, sigma),
        normal_lcdf(Ht | H_hat, sigma)
      );
    }
  }

  // Function to calculate the likelihood for pre-vaccination titers
  // Because pre-vaccination titers are represented in wide-format (repeatedly for each post-vaccination titer value)
  // the pre-vaccination titer for a given person/strain/year will be represented multiple times in the data
  // Here we make sure that the pre-vaccination titer for a given person/strain/year contributes only once regardless of how many
  // post-infection time points were included.
  // (It can contribute twice if there have been two replicated measurements of the same person/strain/year, which happens for FRNT)
  // Supports LOO via holdout_action:
  //   0 = include all units
  //   1 = exclude held-out unit
  //   2 = include only held-out unit
  vector compute_Hpre_loglik(
    int loglik_unit,
    int holdout_index,
    int holdout_action,
    int N_total,
    int I,
    int[] individual,
    int N_Hpre_groups,
    int[] Hpre_group_index,
    vector Hpre,
    vector undetectable_value,
    int max_dilution,
    int censoring,
    int[] titer_type,
    vector sigma_Hpre, 
    vector null_mean_Hpre,
    int[] measurement_replicate,
    int[] strain
  ) {
    // Initialize vector of log-likelihoods for each Hpre_group or for each individual
    // depending on loglik_unit
    // If loglik_unit = 1, this vector has length N_Hpre_groups, if = 2, length = I
    int length_output = (loglik_unit == 1 ? N_Hpre_groups : I);

    vector[length_output] lp = rep_vector(0, length_output);

    // For each group g of individual/strain/year observations...
    for (g in 1:N_Hpre_groups) {
      // locate first observations in this group that are or aren't a repeated measurement
      int i1 = 0;
      int i2 = 0;
      for (i in 1:N_total) {
        if (Hpre_group_index[i] == g) {
          if (measurement_replicate[i] == 1 && i1 == 0) i1 = i;
          if (measurement_replicate[i] == 2 && i2 == 0) i2 = i;
        }
      }
      int rep_measurement_indices[2];
      rep_measurement_indices[1] = i1;
      rep_measurement_indices[2] = i2;

      // Find the individual index at i1 (will be the same in i2)
      int individual_idx = individual[i1];

      // define storing_index. If loglik_unit = 1, it's g, otherwise it's individual_idx
      int storing_index = (loglik_unit == 1 ? g : individual_idx);

      // Skip based on holdout_action
      if ((holdout_action == 1 && storing_index == holdout_index) ||
          (holdout_action == 2 && storing_index != holdout_index)) {
        continue;
      }

      // Iterate across i1 and i2 to calculate the likelihood of 1 or 2 pre-vaccination measurements for this Hpre group
      for (k in 1:2) {
        // Skipping zero handles cases when there was a single measurement
        if (rep_measurement_indices[k] != 0) {
          int rep_idx = rep_measurement_indices[k];
          real L = log(2 * undetectable_value[rep_idx]) / log(2);
          real U = log(max_dilution) / log(2);
          if (censoring == 1) {
            lp[storing_index] += log_Pr_H(Hpre[rep_idx], null_mean_Hpre[strain[rep_idx]], sigma_Hpre[strain[rep_idx]], L, U);
          } else {
            lp[storing_index] += normal_lpdf(Hpre[rep_idx] | null_mean_Hpre[strain[rep_idx]], sigma_Hpre[strain[rep_idx]]);
          }
        }
      }
    }
    return lp;
  }

  // Function to calculate the likelihood for post-vaccination titers
  // Supports LOO via holdout_action (same semantics as above)
  vector compute_postvax_loglik(
    int loglik_unit,
    int holdout_index,
    int holdout_action,
    int N_total,
    int I,
    int[] individual,
    int N_Hpre_groups,
    int[] Hpre_group_index,
    int[] titer_type,
    vector Ht,
    vector undetectable_value,
    int max_dilution,
    int censoring,
    vector null_mean_Ht,
    vector sigma_Ht,
    int[] strain
  ) {

    // Determine output length based on loglik_unit
    int length_output = (loglik_unit == 1 ? N_Hpre_groups : I);
    vector[length_output] lp = rep_vector(0, length_output);

    // Here we get to do a straigthforward loop across all observations,
    // With multiple post-vaccination time points
    // (or repeated measurements of the same post-vax timepoint)
    // represented in long format
    for (n in 1:N_total) {

      // Identify Hpre group (needed to retrieve latent pre-vax titer)
      int g = Hpre_group_index[n];
      int storing_index = (loglik_unit == 1 ? g : individual[n]);

      // Skip based on holdout_action
      if ((holdout_action == 1 && storing_index == holdout_index) ||
          (holdout_action == 2 && storing_index != holdout_index)) {
        continue;
      }

      // Define L (the log2 of the smallest real dilution) as 
      // (note that the smallest real dilution is 2 * undetectable_value)
      real L = log(2 * undetectable_value[n]) / log(2);
      real U = log(max_dilution) / log(2);

      if (censoring == 1) {
        lp[storing_index] += log_Pr_H(Ht[n], null_mean_Ht[strain[n]], sigma_Ht[strain[n]], L, U);
      } else {
        lp[storing_index] += normal_lpdf(Ht[n] | null_mean_Ht[strain[n]], sigma_Ht[strain[n]]);
      }
    }
    return lp;
  }
}

data {
  int<lower=0> N_total;
  int<lower=1> strain[N_total];  // Strain index for each observation
  int<lower=1> K;  // Number of unique strains
  int<lower=1, upper=2> measurement_replicate[N_total];
  vector[N_total] Ht;
  vector[N_total] Hpre;
  int<lower=1> I;
  int<lower=1> individual[N_total];
  int<lower=1> N_Hpre_groups;
  int<lower=1, upper=N_Hpre_groups> Hpre_group_index[N_total];
  vector[N_total] undetectable_value;
  int<lower=1> max_dilution;
  int<lower=0, upper=1> censoring;
  int<lower=1, upper=2> titer_type[N_total];
  int<lower=1, upper=2> loglik_unit;
  int<lower=0> holdout_index;
}

parameters {
  vector<lower = -10, upper = 20>[K] null_mean_Hpre;
  vector<lower = -10, upper = 20>[K] null_mean_Ht;
  vector<lower=0>[K] sigma_Hpre; // strain-specific sigma for Hpre
  vector<lower=0>[K] sigma_Ht;   // strain-specific sigma for Ht
}


model {
  // Null model: all titers are modeled as coming from a single mean (null_mean_Hpre or null_mean_Ht)
  target += sum(compute_Hpre_loglik(
    loglik_unit,
    holdout_index,
    0, // For pre-vaccination titers, we include the held-out unit in the training data
    N_total,
    I,
    individual,
    N_Hpre_groups,
    Hpre_group_index,
    Hpre,
    undetectable_value,
    max_dilution,
    censoring,
    titer_type,
    sigma_Hpre,
    null_mean_Hpre,
    measurement_replicate,
    strain
  ));

  target += sum(compute_postvax_loglik(
    loglik_unit,
    holdout_index,
    1, // exclude post-vaccination titer of the held-out unit from the likelihood
    N_total,
    I,
    individual,
    N_Hpre_groups,
    Hpre_group_index,
    titer_type,
    Ht,
    undetectable_value,
    max_dilution,
    censoring,
    null_mean_Ht,
    sigma_Ht,
    strain
  ));
}

generated quantities {

  // Log predictive density of the held-out post vaccination titer(s)
    real LOO_LPD_post;

    if (holdout_index > 0) {

      vector[loglik_unit == 1 ? N_Hpre_groups : I] loglik_post = compute_postvax_loglik(
        loglik_unit,
        holdout_index,
        2, // include only held-out unit
        N_total,
        I,
        individual,
        N_Hpre_groups,
        Hpre_group_index,
        titer_type,
        Ht,
        undetectable_value,
        max_dilution,
        censoring,
        null_mean_Ht,
        sigma_Ht,
        strain
      );

      // The helper returns a vector that is 0 except for held-out unit; sum() collapses it to a single scalar.
      LOO_LPD_post = sum(loglik_post);

    } else {
      // No holdout requested; skip LOO_LPD computation
      LOO_LPD_post = 0;
    }
}
