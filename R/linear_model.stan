functions {
  real calculate_Hpeak(
    real Hpre,
    real RVE,
    real a,
    real b_peak,
    real k_peak_flu,
    real k_peak_cov2, 
    real u_peak_shared,
    real u_peak_subtype,
    real u_peak_year,
    int infection_before_sample_subtype_matched,
    int recent_infection_before_sample_cov2, 
    int sex,
    real beta_M,
    int age_group,
    real beta_age3140,
    real beta_age4150,
    int BMI_group, 
    real beta_BMI_under,
    real beta_BMI_over,
    int smoking,
    real beta_smoking,
    int asthma,
    real beta_asthma,
    int time_since_vax_group,
    real beta_time_under28,
    real beta_time_over42
  ) {
    real age_effect = (age_group == 2 ? beta_age3140 : (age_group == 3 ? beta_age4150 : 0));
    real bmi_effect = (BMI_group == 2 ? beta_BMI_under : (BMI_group == 3 ? beta_BMI_over : 0));
    real time_effect = (time_since_vax_group == 2 ? beta_time_under28 : (time_since_vax_group == 3 ? beta_time_over42 : 0));
    real Hpeak = (a * (1 + RVE)) +
                     (1 - b_peak * (1 + RVE)) * Hpre +
                     k_peak_flu * infection_before_sample_subtype_matched +
                     k_peak_cov2 * recent_infection_before_sample_cov2 +
                     u_peak_shared + u_peak_subtype + u_peak_year +
                     (sex == 1 ? beta_M : 0) +
                     age_effect +
                     bmi_effect +
                     smoking * beta_smoking +
                     asthma * beta_asthma +
                     time_effect;
    return Hpeak;
  }

  real calculate_HLT(
    real Hpeak,
    real Hpre,
    real f
  ) {
    return Hpre + f * (Hpeak - Hpre);
  }

  real calculate_H_hat(
    real Hpre,
    real Hpeak,
    real HLT,
    real t,
    real tau,
    real omega
  ) {
    if (t <= tau) {
      return Hpre + (Hpeak - Hpre) * (t / tau);
    } else {
      return HLT + (Hpeak - HLT) * exp(-omega * (t - tau));
    }
  }

  real calculate_f(
    real f0,
    real r_LT,
    real b_LT,
    real Hpre,
    real n_previous_vax
  ) {
    real logit_f = logit(f0) + r_LT * n_previous_vax + b_LT * Hpre;
    return inv_logit(logit_f);
  }

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
    if (fabs(Ht - (L - 1)) < 1e-5) {
      return normal_lcdf(L | H_hat, sigma);
    }
    // Special case: Ht = U → P(H_hat ≥ U)
    else if (fabs(Ht - U) < 1e-5) {
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

  // Function to calculate repeat vaccination effect
  real calculate_RVE(
    real r_Y1, real r_Y2, real r_Y3,
    int vax_Y1, int vax_Y2, int vax_Y3,
    real gamma,
    int year
  ) {
    // RVE = sum of vax_Yk * r_Yk * (1 + gamma)^(year - k - 1) for k < year of current observation
    // cases k >= year are handled by r_Y being set to zero.
    // so current and future years don't count toward RVEs.
    return vax_Y1 * r_Y1 * pow(1 + gamma, year - 1 - 1) +
           vax_Y2 * r_Y2 * pow(1 + gamma, year - 2 - 1) +
           vax_Y3 * r_Y3 * pow(1 + gamma, year - 3 - 1);
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
    vector L,
    real U,
    int censoring,
    int[] titer_type,
    real sigma_HAI,
    real sigma_FRNT,
    vector Hpre_continuous,
    int[] measurement_replicate
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
          real sigma = titer_type[rep_idx] == 1 ? sigma_HAI : sigma_FRNT;
          if (censoring == 1) {
            lp[storing_index] += log_Pr_H(Hpre[rep_idx], Hpre_continuous[g], sigma, L[rep_idx], U);
          } else {
            lp[storing_index] += normal_lpdf(Hpre[rep_idx] | Hpre_continuous[g], sigma);
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
    int[] strain,
    int[] ind_subtype_combination,
    int[] ind_year_combination,
    int[] infection_before_sample_subtype_matched,
    int[] recent_infection_before_sample_cov2,
    int[] sex,
    int[] age_group,
    int[] BMI_group,
    int[] smoking,
    int[] asthma,
    int[] time_since_vax_group,
    int[] titer_type,
    int[] vax_Y1,
    int[] vax_Y2,
    int[] vax_Y3,
    int[] year,
    int[] RVE_level_Y1,
    int[] RVE_level_Y2,
    int[] RVE_level_Y3,
    vector Ht,
    vector t,
    vector L,
    real U,
    int censoring,
    vector n_previous_vax,
    vector Hpre_continuous,
    vector a,
    vector b_peak,
    real k_peak_flu,
    vector k_peak_cov2,
    vector u_peak_shared,
    vector u_peak_subtype,
    vector u_peak_year,
    real beta_M,
    real beta_age3140,
    real beta_age4150,
    real beta_BMI_under,
    real beta_BMI_over,
    real beta_smoking,
    real beta_asthma,
    real beta_time_under28,
    real beta_time_over42,
    real sigma_HAI,
    real sigma_FRNT,
    real f0,
    real r_LT,
    real b_LT,
    real tau,
    real omega,
    vector r,
    real gamma
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

      // Determine which sigma applies 
      real sigma = titer_type[n] == 1 ? sigma_HAI : sigma_FRNT;

      // Compute the total repeat vaccination effect
      real RVE = calculate_RVE(
        r[RVE_level_Y1[n]],
        r[RVE_level_Y2[n]],
        r[RVE_level_Y3[n]],
        vax_Y1[n],
        vax_Y2[n],
        vax_Y3[n],
        gamma,
        year[n]
      );

      // Compute the peak titer
      real Hpeak = calculate_Hpeak(
        // Find the sampled continuous pre-vaccination titer associated with this observation
        // based on the individual-strain-year combination the observation is in, i.e. Hpre_group_index[n]
        Hpre_continuous[g], RVE, a[strain[n]], b_peak[strain[n]],
        k_peak_flu, k_peak_cov2[sex[n]],
        u_peak_shared[individual[n]],
        u_peak_subtype[ind_subtype_combination[n]],
        u_peak_year[ind_year_combination[n]],
        infection_before_sample_subtype_matched[n],
        recent_infection_before_sample_cov2[n],
        sex[n],
        beta_M,
        age_group[n],
        beta_age3140,
        beta_age4150,
        BMI_group[n],
        beta_BMI_under,
        beta_BMI_over,
        smoking[n],
        beta_smoking,
        asthma[n],
        beta_asthma,
        time_since_vax_group[n],
        beta_time_under28,
        beta_time_over42
      );

      // Long-term kinetics
      real f = calculate_f(f0, r_LT, b_LT, Hpre_continuous[Hpre_group_index[n]], n_previous_vax[n]);
      real HLT = calculate_HLT(Hpeak, Hpre_continuous[Hpre_group_index[n]], f);

      // Final predicted post-vaccination titer for this observation
      real H_hat = calculate_H_hat(Hpre_continuous[Hpre_group_index[n]], Hpeak, HLT, t[n], tau, omega);

      if (censoring == 1) {
        lp[storing_index] += log_Pr_H(Ht[n], H_hat, sigma, L[n], U);
      } else {
        lp[storing_index] += normal_lpdf(Ht[n] | H_hat, sigma);
      }
    }
    return lp;
  }
}

data {
  int<lower=0> N_total;  // Total number of observations
  int<lower=1, upper=2> measurement_replicate[N_total]; // Indicates whether each observation is the first or possibly second repeated measurement of the same person/strain/year
  vector[N_total] Ht;  // Observed log-transformed antibody titers at time t
  vector[N_total] Hpre;  // Log-transformed pre-vaccination antibody titers
  int<lower=1> strain[N_total];  // Strain index for each observation
  int<lower=1> K;  // Number of unique strains
  vector[N_total] n_previous_vax;  // Number of prior vaccinations for each observation
  int<lower=1> I;  // Number of unique individuals
  int<lower=1> individual[N_total];  // Individual IDs
  vector[N_total] t;  // Days since yearly vaccination
  vector[N_total] undetectable_value;  // Value assigned to undetectable titer (vector, handles FRNT vs HAI)
  int<lower=1> max_dilution;  // Maximum dilution level (integer)
  int<lower=0, upper=1> censoring;  // Boolean switch for censoring (1 = apply censoring, 0 = use raw values)
  int<lower=0, upper=1> infection_before_sample_subtype_matched[N_total];  // Indicator for prior infection before sample date
  int<lower=0, upper=1> recent_infection_before_sample_cov2[N_total];
  int<lower=1> N_ind_sub_combinations;  // Number of unique individual-subtype combinations
  int<lower=1> ind_subtype_combination[N_total];  // Individual-subtype IDs
  int<lower=0, upper=1> include_upeak_subtype;  // Switch for including subtype-specific individual effects
  int<lower=0, upper=1> include_upeak_shared; // Switch for including shared individual effects
  int<lower=0, upper=1> include_upeak_year; // Switch for including shared year-specific individual effects
  int<lower=0, upper=1> include_kinetics; //
  int<lower=0, upper=1> include_RVE; // Switch including repeat vaccination effects
  int<lower=1> N_Hpre_groups; // Number of unique individual-strain-year groups (we'll call them Hpre groups because all obs in each group have same pre-vac titer)
  int<lower=1, upper=N_Hpre_groups> Hpre_group_index[N_total]; // Mapping from observations to individual-strain-year combination
  int<lower=1> N_strain_npriorvax_groups; // Number of unique strain/n_priorvax combinations
  int<lower=1, upper=N_strain_npriorvax_groups> strain_npriorvax_group_index[N_total]; // Index linking each observation to a strain/npriorvax group.
  int<lower=1, upper=2> sex[N_total]; // Sex for each observation (1=male, 2=female)
  int<lower=1> N_ind_year_combinations;  // Number of unique individual-year combinations
  int<lower=1, upper=N_ind_year_combinations> ind_year_combination[N_total];  // Individual-year IDs
  int<lower=1, upper=3> age_group[N_total]; // Age group for each observation (1="18-30", 2="31-40", 3="41-50")
  int<lower=1> N_RVE_effects; // Number of repeat vaccination effect parameters.
  int<lower=1, upper=N_RVE_effects+1> RVE_level_Y1[N_total];
  int<lower=1, upper=N_RVE_effects+1> RVE_level_Y2[N_total];
  int<lower=1, upper=N_RVE_effects+1> RVE_level_Y3[N_total];
  int<lower=0, upper=1> include_k_peak_flu; // Switch for including k_peak_flu parameter
  int<lower=0, upper=1> include_k_peak_cov2; // Switch for including k_peak_cov2 parameter
  int<lower=0, upper=1> smoking[N_total]; // Smoking status (logical)
  int<lower=0, upper=1> asthma[N_total];  // Asthma status (logical)
  int<lower=1, upper=2> titer_type[N_total]; // 1=HAI, 2=FRNT for each observation
  int<lower=1, upper=3> BMI_group[N_total]; // BMI group mapping (1=healthy,2=under,3=over)
  // Vaccination indicators for each year
  int<lower=0, upper=1> vax_Y1[N_total];
  int<lower=0, upper=1> vax_Y2[N_total];
  int<lower=0, upper=1> vax_Y3[N_total];
  int<lower=1> year[N_total]; // Add year as data for each observation
  int<lower=0, upper=1> RVE_decay; // Switch for including RVE decay parameter
  int<lower=1, upper=2> loglik_unit; // How to factorize log-lik. 1 = Hpre group, 2 = individual
  int<lower=0> holdout_index; // Index of held-out unit: 0=no holdout; else 1..N_Hpre_groups if loglik_unit==1, or 1..I
  vector[N_strain_npriorvax_groups] prior_mean_Hpre_mean; // Prior means for Hpre_mean
  vector[N_strain_npriorvax_groups] prior_mean_Hpre_sd; // Prior means for Hpre_sd 
  vector[N_strain_npriorvax_groups] prior_sd_Hpre_mean; // Prior SD for Hpre_mean
  vector[N_strain_npriorvax_groups] prior_sd_Hpre_sd; // Prior SD for Hpre_sd
  int<lower=0, upper=1> include_Hpre_priors; // Switch for assigning priors to Hpre_mean and Hpre_sd
  int<lower=0, upper=1> include_time_since_vax; // Switch for including time-since-vaccination effects
  int<lower=1, upper=3> time_since_vax_group[N_total]; // Time-since-vaccination group (1=reference, 2=under28, 3=over42)
  int<lower=0, upper=1> include_sex; // Switch for including sex effect
  int<lower=0, upper=1> include_age; // Switch for including age group effects
  int<lower=0, upper=1> include_BMI; // Switch for including BMI group effects
  int<lower=0, upper=1> include_smoking; // Switch for including smoking effect
  int<lower=0, upper=1> include_asthma; // Switch for including asthma effect
}

transformed data {
  
  // Calculate L and U (log2 values for the lowest and highest measured titers)
  vector[N_total] L;
  real U;
  for (n in 1:N_total) {
    // Define L (the log2 of the smallest real dilution) as 
    // (note that the smallest real dilution is 2 * undetectable_value)
    L[n] = log(2 * undetectable_value[n]) / log(2);
  }
  U = log(max_dilution) / log(2);
}

parameters {
  vector<lower=0>[K] a;  // Positive intercepts for each strain
  vector<lower=0, upper =1>[K] b_peak;  // Baseline antibody ceiling effect for each strain
  vector<lower=-1,upper=1>[include_RVE ? N_RVE_effects : 0] r;  // Repeat vaccination effects
  real<lower=0> sigma_HAI;  // Residual SD for HAI
  real<lower=0> sigma_FRNT; // Residual SD for FRNT
  real<lower=-5, upper=5> k_peak_flu[include_k_peak_flu ? 1 : 0];  // Effect of prior infection on peak titer
  real<lower=-5, upper=5> k_peak_cov2[include_k_peak_cov2 ? 2 : 0]; // Sex-specific effect of prior cov2 infection on peak titer
  real beta_M[include_sex ? 1 : 0]; // Effect of male sex on Hpeak
  real beta_age3140[include_age ? 1 : 0]; // Effect of age group 31-40 on Hpeak
  real beta_age4150[include_age ? 1 : 0]; // Effect of age group 41-50 on Hpeak
  real beta_smoking[include_smoking ? 1 : 0]; // Effect of smoking on Hpeak
  real beta_asthma[include_asthma ? 1 : 0];  // Effect of asthma on Hpeak
  real beta_BMI_under[include_BMI ? 1 : 0]; // Effect of underweight BMI relative to healthy
  real beta_BMI_over[include_BMI ? 1 : 0];  // Effect of overweight BMI relative to healthy
  real<lower=-1, upper=1> gamma[RVE_decay ? 1 : 0]; // Decay parameter for RVE if enabled
  real beta_time_under28[include_time_since_vax ? 1 : 0]; // Effect of time-since-vax group 2 (under 28 days) on Hpeak
  real beta_time_over42[include_time_since_vax ? 1 : 0];  // Effect of time-since-vax group 3 (over 42 days) on Hpeak
  // ===== The parameters below will only be included if include_upeak_shared == 1 =====
  // The ternary statements (condition ? value_if_true: value_if_false) make it so that 
  // the parameters have length 0 (and thus won't be estimated) if include_upeak_subtype == 0
  // (https://www.martinmodrak.cz/2018/04/24/optional-parameters-data-in-stan/) 
  real<lower=0> sigma_upeak_shared[include_upeak_shared ? 1 : 0];  // SD of individual effects on peak titer
  vector[include_upeak_shared ? I : 0] u_peak_shared;  // Individual random effects on peak titer
  // ===== The parameters below will only be included if include_upeak_subtype == 1 =====
  real<lower=0> sigma_upeak_subtype[include_upeak_subtype ? 1 : 0];  // SD of subtype-specific individual effects on peak titer
  vector[include_upeak_subtype ? N_ind_sub_combinations : 0] u_peak_subtype;  // Subtype-specific individual random effects on peak titer
  // ===== The parameters below will only be included if include_upeak_year == 1 =====
  real<lower=0> sigma_upeak_year[include_upeak_year ? 1 : 0];  // SD of individual-year effects on peak titer
  vector[include_upeak_year ? N_ind_year_combinations : 0] u_peak_year;  // Individual-year random effects on peak titer
  // ===== The parameters below will only be estimated if include_kinetics is TRUE in the data block ===
  // if it isn't, we'll set tau, b_LT, r_LT and omega to 0 and f0 to 1 in all calculations
  // (titers immediately rise to the peak value and stay there forever)
  real<lower=0, upper=1> f0[include_kinetics ? 1 : 0];  // Shared baseline long-term fraction
  real<lower=-10, upper=10> r_LT[include_kinetics ? 1 : 0];  // Shared repeat vaccination effect on long-term fraction
  real<lower=-1, upper=1> b_LT[include_kinetics ? 1 : 0];  // Shared effect of Hpre on long-term fraction
  real<lower=0, upper=90> tau[include_kinetics ? 1 : 0];  // Shared baseline tau
  real<lower=0, upper=0.01> omega[include_kinetics ? 1 : 0];  // Baseline omega
  // ========
  vector<lower = -10, upper = 20>[N_Hpre_groups] Hpre_continuous; // Continuous pre-vaccination titer for each individual-strain-year combination
  vector[N_strain_npriorvax_groups] Hpre_mean; // Prior mean for Hpre_continuous
  vector<lower=0>[N_strain_npriorvax_groups] Hpre_sd; // Prior sd for Hpre_continuous
}

transformed parameters {
  // Defines generalized (transformed) versions of u_peak_subtype and sigma_upeak_subtype
  // to handle instances where include_upeak_subtype == 1 or == 0

  vector[I] u_peak_shared_transformed;
  real<lower=0> sigma_upeak_shared_transformed;

  if (include_upeak_shared == 1) {
    u_peak_shared_transformed = u_peak_shared;
    sigma_upeak_shared_transformed = sigma_upeak_shared[1];
  } else {
    for (i in 1:I) {
      u_peak_shared_transformed[i] = 0;
    }
    sigma_upeak_shared_transformed = 0;
  }

  vector[N_ind_sub_combinations] u_peak_subtype_transformed;
  real<lower=0> sigma_upeak_subtype_transformed;

  if (include_upeak_subtype == 1) {
    u_peak_subtype_transformed = u_peak_subtype;
    sigma_upeak_subtype_transformed = sigma_upeak_subtype[1];
  } else {
    for( n in 1:N_ind_sub_combinations ){
      u_peak_subtype_transformed[n] = 0;
    }
    sigma_upeak_subtype_transformed = 0;
  }

  vector[N_ind_year_combinations] u_peak_year_transformed;
  real<lower=0> sigma_upeak_year_transformed;

  if (include_upeak_year == 1) {
    u_peak_year_transformed = u_peak_year;
    sigma_upeak_year_transformed = sigma_upeak_year[1];
  } else {
    for (n in 1:N_ind_year_combinations) {
      u_peak_year_transformed[n] = 0;
    }
    sigma_upeak_year_transformed = 0;
  }

  real f0_transformed;
  real tau_transformed;
  real b_LT_transformed;
  real r_LT_transformed;
  real omega_transformed;

  // Define the transformed parameters for kinetics
  if (include_kinetics == 1){
    tau_transformed = tau[1];
    f0_transformed = f0[1];
    b_LT_transformed = b_LT[1];
    r_LT_transformed = r_LT[1];
    omega_transformed = omega[1];
  } else {
    tau_transformed = 0;
    f0_transformed = 1;
    b_LT_transformed = 0;
    r_LT_transformed = 0;
    omega_transformed = 0;
  }
  real k_peak_flu_transformed;
  k_peak_flu_transformed = include_k_peak_flu ? k_peak_flu[1] : 0;

  vector[2] k_peak_cov2_transformed;
  if (include_k_peak_cov2 == 1) {
    k_peak_cov2_transformed[1] = k_peak_cov2[1];
    k_peak_cov2_transformed[2] = k_peak_cov2[2];
  } else {
    k_peak_cov2_transformed = rep_vector(0, 2);
  }
  real gamma_transformed;
  if (RVE_decay == 1) {
    gamma_transformed = gamma[1];
  } else {
    gamma_transformed = 0;
  }

  real beta_time_under28_transformed;
  real beta_time_over42_transformed;
  if (include_time_since_vax == 1) {
    beta_time_under28_transformed = beta_time_under28[1];
    beta_time_over42_transformed = beta_time_over42[1];
  } else {
    beta_time_under28_transformed = 0;
    beta_time_over42_transformed = 0;
  }

  real beta_M_transformed;
  beta_M_transformed = include_sex ? beta_M[1] : 0;

  real beta_age3140_transformed;
  real beta_age4150_transformed;
  beta_age3140_transformed = include_age ? beta_age3140[1] : 0;
  beta_age4150_transformed = include_age ? beta_age4150[1] : 0;

  real beta_BMI_under_transformed;
  real beta_BMI_over_transformed;
  beta_BMI_under_transformed = include_BMI ? beta_BMI_under[1] : 0;
  beta_BMI_over_transformed = include_BMI ? beta_BMI_over[1] : 0;

  real beta_smoking_transformed;
  beta_smoking_transformed = include_smoking ? beta_smoking[1] : 0;

  real beta_asthma_transformed;
  beta_asthma_transformed = include_asthma ? beta_asthma[1] : 0;

  // Transformed repeat vaccination effect vector
  // We add a value at the end that will prevent vaccination in current and future years
  // from counting toward RVEs in the current year of each observation
  // For instance, if year = 2, RVE_level_Y2 and RVE_level_Y3 will have an integer value N_RVE_effects + 1
  // So, for observation n, r_transformed[RVE_level_Y2[n]] = r_transformed[RVE_level_Y3[n]] r_transformed[N_RVE_effects + 1] = 0
  // (vaccination status in year 3 and year 2 itself won't count toward repeat vaccination effects in year 2)
  vector[N_RVE_effects + 1] r_transformed;

  // Initialize to zeros; fill from r only if include_RVE == 1.
  for (i in 1:(N_RVE_effects + 1)) {
    r_transformed[i] = 0;
  }
  if (include_RVE == 1) {
    for (i in 1:N_RVE_effects) {
      r_transformed[i] = r[i];
    }
  }
  
}

model {
  
  if(include_kinetics){
      // Prior for tau
      tau ~ uniform(0, 90);

      // Prior for omega
      //omega ~ uniform(0,0.05);

      // Prior for f0
      //f0 ~ beta(8, 8);  

      // Prior for r_LT
      //r_LT ~ normal(0, 1);

  }

  // Conditionally assign priors to Hpre_mean and Hpre_sd
  if (include_Hpre_priors == 1) {
    // Normal prior for Hpre_mean, centered on provided values
    Hpre_mean ~ normal(prior_mean_Hpre_mean, prior_sd_Hpre_mean);

    // Normal prior for Hpre_sd, centered on provided values
    Hpre_sd ~ normal(prior_mean_Hpre_sd, prior_sd_Hpre_sd);
  }

  // Prior for individual random effects
  if (include_upeak_shared == 1) {
    u_peak_shared ~ normal(0, sigma_upeak_shared[1]);
  }
  if (include_upeak_subtype == 1) {
    // Prior for subtype-specific individual random effects
    u_peak_subtype ~ normal(0, sigma_upeak_subtype[1]);
  }
  if (include_upeak_year == 1) {
    u_peak_year ~ normal(0, sigma_upeak_year[1]);
  }

  // ======== Prior for latent pre-vaccination titer ========
  // For each group g of individual/strain/year observations
  for (g in 1:N_Hpre_groups) {
    // Find the first observation in the group
    // (if more than one post-vaccination time point was included, there will be multiple post-vaccination observations, but
    // they will have the same latent pre-vaccination titer. We find the first observation in the group to retrieve its
    // strain/n. prior vax combination, which in turns determines the mean and sd of the latent pre-vaccination titer)
    int i = 1;
    while (Hpre_group_index[i] != g) i += 1;
    int prior_idx = strain_npriorvax_group_index[i];
    // Sample latent pre-vaccination titer 
    Hpre_continuous[g] ~ normal(Hpre_mean[prior_idx], Hpre_sd[prior_idx]);
  }

  // ======== Likelihood contribution of observed pre-vaccination titers ========
  // compute_Hpre_loglik returns a vector by Hpre_group or individual (for later cross validation)
  // we sum to get a total contribution
  target += sum(compute_Hpre_loglik(
    loglik_unit,
    holdout_index,
    0, // For pre-vaccination titers, we include the held-out unit in the fitting 
    N_total,
    I,
    individual,
    N_Hpre_groups,
    Hpre_group_index,
    Hpre,
    L,
    U,
    censoring,
    titer_type,
    sigma_HAI,
    sigma_FRNT,
    Hpre_continuous,
    measurement_replicate
  ));

  // ======== Likelihood for post-vaccination titers ========
  // compute_postvax_loglik returns a vector by Hpre_group or individual (for later cross validation)
  // we sum to get a total contribution
  target += sum(compute_postvax_loglik(
    loglik_unit,
    holdout_index,
    1, // For post-vaccination titers, we exclude the held-out unit from the fitting.
    N_total,
    I,
    individual,
    N_Hpre_groups,
    Hpre_group_index,
    strain,
    ind_subtype_combination,
    ind_year_combination,
    infection_before_sample_subtype_matched,
    recent_infection_before_sample_cov2,
    sex,
    age_group,
    BMI_group,
    smoking,
    asthma,
    time_since_vax_group,
    titer_type,
    vax_Y1,
    vax_Y2,
    vax_Y3,
    year,
    RVE_level_Y1,
    RVE_level_Y2,
    RVE_level_Y3,
    Ht,
    t,
    L,
    U,
    censoring,
    n_previous_vax,
    Hpre_continuous,
    a,
    b_peak,
    k_peak_flu_transformed,
    k_peak_cov2_transformed,
    u_peak_shared_transformed,
    u_peak_subtype_transformed,
    u_peak_year_transformed,
    beta_M_transformed,
    beta_age3140_transformed,
    beta_age4150_transformed,
    beta_BMI_under_transformed,
    beta_BMI_over_transformed,
    beta_smoking_transformed,
    beta_asthma_transformed,
    beta_time_under28_transformed,
    beta_time_over42_transformed,
    sigma_HAI,
    sigma_FRNT,
    f0_transformed,
    r_LT_transformed,
    b_LT_transformed,
    tau_transformed,
    omega_transformed,
    r_transformed,
    gamma_transformed
  ));
}

generated quantities {

  // Predicted titers, including individual effects
  vector[N_total] H_hat;

  // Predicted titers given only fixed effects
  vector[N_total] H_hat_fixed;

  // LOO log-likelihood (log predictive density) for the post-vaccination titer of the held-out unit
  real LOO_LPD_post;

  // Export predicted titer H_hat for all observations
  for (n in 1:N_total) {
    real f = calculate_f(f0_transformed, r_LT_transformed, b_LT_transformed, Hpre_continuous[Hpre_group_index[n]], n_previous_vax[n]);
    real RVE = calculate_RVE(
      r_transformed[RVE_level_Y1[n]],
      r_transformed[RVE_level_Y2[n]],
      r_transformed[RVE_level_Y3[n]],
      vax_Y1[n],
      vax_Y2[n],
      vax_Y3[n],
      gamma_transformed,
      year[n]
    );
    real Hpeak = calculate_Hpeak(
      Hpre_continuous[Hpre_group_index[n]], RVE, a[strain[n]], b_peak[strain[n]],
      k_peak_flu_transformed, k_peak_cov2_transformed[sex[n]],
      u_peak_shared_transformed[individual[n]],
      u_peak_subtype_transformed[ind_subtype_combination[n]],
      u_peak_year_transformed[ind_year_combination[n]],
      infection_before_sample_subtype_matched[n],
      recent_infection_before_sample_cov2[n],
      sex[n], beta_M_transformed,
      age_group[n], beta_age3140_transformed, beta_age4150_transformed,
      BMI_group[n], beta_BMI_under_transformed, beta_BMI_over_transformed,
      smoking[n], beta_smoking_transformed,
      asthma[n], beta_asthma_transformed,
      time_since_vax_group[n], beta_time_under28_transformed, beta_time_over42_transformed);
    // Fixed-effect predictions obtained by setting individual effects to zero
    real Hpeak_fixed = calculate_Hpeak(
      Hpre_continuous[Hpre_group_index[n]], RVE, a[strain[n]], b_peak[strain[n]],
      k_peak_flu_transformed, k_peak_cov2_transformed[sex[n]],
      0, 0, 0,
      infection_before_sample_subtype_matched[n],
      recent_infection_before_sample_cov2[n],
      sex[n], beta_M_transformed,
      age_group[n], beta_age3140_transformed, beta_age4150_transformed,
      BMI_group[n], beta_BMI_under_transformed, beta_BMI_over_transformed,
      smoking[n], beta_smoking_transformed,
      asthma[n], beta_asthma_transformed,
      time_since_vax_group[n], beta_time_under28_transformed, beta_time_over42_transformed);
    real HLT = calculate_HLT(Hpeak, Hpre_continuous[Hpre_group_index[n]], f);
    real HLT_fixed = calculate_HLT(Hpeak_fixed, Hpre_continuous[Hpre_group_index[n]], f); // Long term titer given fixed-effects Hpeak

    H_hat[n] = calculate_H_hat(Hpre_continuous[Hpre_group_index[n]], Hpeak, HLT, t[n], tau_transformed, omega_transformed);

    H_hat_fixed[n] = calculate_H_hat(Hpre_continuous[Hpre_group_index[n]], Hpeak_fixed, HLT_fixed, t[n], tau_transformed, omega_transformed);
  }

  // Compute predictive densitity (log-likelihood) only for the held-out unit
  if (holdout_index > 0) {
    
    // We compute LPD only for the post-vaccination titers
    vector[loglik_unit == 1 ? N_Hpre_groups : I] loglik_post = compute_postvax_loglik(
      loglik_unit,
      holdout_index,
      2, // include only held-out unit
      N_total,
      I,
      individual,
      N_Hpre_groups,
      Hpre_group_index,
      strain,
      ind_subtype_combination,
      ind_year_combination,
      infection_before_sample_subtype_matched,
      recent_infection_before_sample_cov2,
      sex,
      age_group,
      BMI_group,
      smoking,
      asthma,
      time_since_vax_group,
      titer_type,
      vax_Y1,
      vax_Y2,
      vax_Y3,
      year,
      RVE_level_Y1,
      RVE_level_Y2,
      RVE_level_Y3,
      Ht,
      t,
      L,
      U,
      censoring,
      n_previous_vax,
      Hpre_continuous,
      a,
      b_peak,
      k_peak_flu_transformed,
      k_peak_cov2_transformed,
      u_peak_shared_transformed,
      u_peak_subtype_transformed,
      u_peak_year_transformed,
      beta_M_transformed,
      beta_age3140_transformed,
      beta_age4150_transformed,
      beta_BMI_under_transformed,
      beta_BMI_over_transformed,
      beta_smoking_transformed,
      beta_asthma_transformed,
      beta_time_under28_transformed,
      beta_time_over42_transformed,
      sigma_HAI,
      sigma_FRNT,
      f0_transformed,
      r_LT_transformed,
      b_LT_transformed,
      tau_transformed,
      omega_transformed,
      r_transformed,
      gamma_transformed
    );

    // The helpers returns a vector that is 0 except for held-out unit; sum() collapses it to a single scalar.
    LOO_LPD_post = sum(loglik_post);

  } else {
    // No holdout requested; skip LOO_LPD computation
    LOO_LPD_post = 0;
  }
}
