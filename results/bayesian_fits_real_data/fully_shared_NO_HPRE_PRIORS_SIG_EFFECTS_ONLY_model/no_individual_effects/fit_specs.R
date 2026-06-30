fit_specs <- list(
    RVE_model = fully_shared_model,
    include_upeak_shared = 0,
    include_upeak_subtype = 0,
    include_upeak_year = 0,
    include_k_peak_flu = 0, # Not significant in fit with default covariate set
    include_k_peak_cov2 = 0,
    include_Hpre_priors = 0,
    include_NAI_infections = 0,
    include_time_since_vax = F,
    include_sex = T,
    include_age = T,
    include_BMI = F,
    include_smoking = F,
    include_asthma = T
)