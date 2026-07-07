library(RColorBrewer)

original_dilutions <- c(5, 10, 20, 40, 80, 160, 320, 640, 1280) # Original dilutions, including the LOD for HAI (5) and FRNT (10)


# ==== Functions ====

# Snap a Date vector to the first day (Monday) of its ISO week, preserving NAs.
# Intended for plotting only, not for calculations.
snap_to_iso_week_start <- function(d) {
  iso <- ISOweek::ISOweek(d)
  ISOweek::ISOweek2date(ifelse(is.na(iso), NA_character_, paste0(iso, "-1")))
}

annotate_calendar_year <- function(data){
  calendar_year_link <- 
    bind_rows(
      tibble(drive = 1, year = 1:5, calendar_year = 2020:2024),
      tibble(drive = 2, year = 1:4, calendar_year = 2021:2024)
    ) %>%
    mutate(next_year = as.integer(str_sub(calendar_year, 3, 4)) + 1,
           calendar_year = paste0(calendar_year, "-", next_year)) %>%
    select(-next_year)

  stopifnot(all(c("drive", "year") %in% names(data)))

  return(
    data %>%
      left_join(calendar_year_link, by = c("drive", "year"))
  )

}

strain_labeller <- function(strain){
  case_when(
    is.na(strain) ~ NA_character_,
    strain == "A/Hong Kong/45/2019" ~ "HK (H3)",
    strain == "A/Cambodia/e0826360/2020" ~ "Camb (H3)",
    strain ==  "A/Darwin/9/2021" | strain == "A/Darwin/6/2021" ~ "Darwin (H3)",
    strain == "A/Hawaii/70/2019" ~ "Hawaii (H1)",
    strain == "A/Wisconsin/588/2019" ~ "Wis588 (H1)",
    strain == "A/Wisconsin/67/2022" ~ "Wis67 (H1)",
    strain == "B/Washington/02/2019" ~ "Wash (B/Vic)",
    strain == "B/Austria/1359417/2021" ~ "Austria (B/Vic)",
    strain == "B/Phuket/3073/2013" ~ "Phuk (B/Yam)"
  ) %>%
  factor(levels = c("HK (H3)", "Camb (H3)", "Darwin (H3)",
                    "Hawaii (H1)", "Wis588 (H1)", "Wis67 (H1)",
                    "Wash (B/Vic)", "Austria (B/Vic)", "Phuk (B/Yam)"))
}

subtype_labeller <- function(strain){
  case_when(
    strain %in% c("A/Hong Kong/45/2019", "A/Cambodia/e0826360/2020", "A/Darwin/9/2021", "A/Darwin/6/2021", "A/Minnesota/41/2019", "A/Tasmania/503/2020") ~ "H3N2",
    strain %in% c("A/Hawaii/66/2019", "A/Hawaii/70/2019", "A/Wisconsin/588/2019", "A/Victoria/2570/2019", "A/Wisconsin/67/2022", "A/Wisconsin/50/2022") ~ "H1N1",
    strain %in% c("B/Phuket/3073/2013") ~ "B/Yamagata",
    strain %in% c("B/Washington/02/2019", "B/Austria/1359417/2021") ~ "B/Victoria"
  )
}

vax_pair_labeller <- function(strain_pair){
  case_when(
    is.na(strain_pair) ~ NA_character_,
    strain_pair %in% vax_strain_pairs$strain_pair ~
      vax_strain_pairs$pair_label[match(strain_pair, vax_strain_pairs$strain_pair)],
    TRUE ~ as.character(strain_pair)
  ) %>%
  factor(levels = unique(vax_strain_pairs$pair_label))
}
  
get_age_groups <- function(data){
  data %>%
    mutate(age_group = case_when(
      age >= 18 & age <= 30 ~ '18-30',
      age > 30 & age <= 40 ~ '31-40',
      age > 40 & age <= 50 ~ '41-50'
    )) %>%
    mutate(age_group = factor(age_group, levels = c('18-30','31-40','41-50')))
}

flag_vaccine_strains <- function(data){
  data %>%
    mutate(is_vaccine_strain =
             (paste(year, strain) %in%
                paste(vaccine_strains$year, vaccine_strains$strain)))
}

p_value_labeller <- function(tibble_with_p_values){
  tibble_with_p_values %>%
    mutate(comparison_with_placebo = 
             (!str_detect(treatment_1, "V") | ! str_detect(treatment_2, "V"))) %>%
    mutate(plot_label = case_when(
      # Comparisons with placebo will be indicated by line type, so we give empty label
      comparison_with_placebo & p_value > 0.05 ~ "",
      # For comparisons between vaccinated groups, label if they are significant 
      !comparison_with_placebo & p_value < 0.0001 ~ "***",
      !comparison_with_placebo & p_value < 0.01 ~ "**",
      !comparison_with_placebo & p_value < 0.05 ~ "*"
    )) %>%
    filter(!is.na(plot_label))
  
}

relabel_timepoint <- function(timepoint){
  timepoint_levels <- paste0("Day ", sort(unique(timepoint)))
  timepoint <- paste0("Day ", timepoint)
  return(factor(timepoint, levels = timepoint_levels))
}

relabel_year <- function(year){
  year_levels <- paste0("Year ", sort(unique(year)))
  year <- paste0("Year ", year)
  return(factor(year, levels = year_levels))
}

get_continuous_timepoints <- function(data){
  data %>%
    mutate(continuous_timepoint = (year - 1) * 365 + timepoint)
}

get_continuous_timepoint_labels <- function(continuous_timepoint){
  tibble(continuous_timepoint = continuous_timepoint) %>%
    mutate(year = floor(continuous_timepoint / 365) + 1,
           timepoint = continuous_timepoint - (year - 1) * 365) %>%
    mutate(label = paste0('Y', year, '\nd', timepoint)) %>%
    pull(label)
}

get_log_base <- function(response_var){
  if(response_var == 'titer'){
    log_base <- 2
  }else{
    stopifnot(response_var == 'luminex')
    log_base <- 10
  }
  return(log_base)
}

# Checks if a reported titer value is among a hard-coded set of dilutions, or in the set of their parwise geometric means
check_titer_values <- function(titer){
  non_na_values <- titer[!is.na(titer)]
  for(x in non_na_values){
    closest_expected_GMT <- possible_pairwise_GMTs[abs(x - possible_pairwise_GMTs) == min(abs(x - possible_pairwise_GMTs))]
    stopifnot(abs(x - closest_expected_GMT) < 1e-4)
  }
}

censor_titers <- function(titer, assay){

  # Titer: a numeric vector of uncensored titer values
  # Assay: a character vector of the same length, with values either "HAI" or "FRNT"

  stopifnot(length(titer) == length(assay))  

  converted_titers <- c()
  for(i in 1:length(titer)){

    value <- titer[i]

    if(str_detect(assay[i], "FRNT")){
      # This handles the fact that a titer of, say, 3 will be coded as 5 for HAI but 10 for FRNT.
      dilutions <- original_dilutions[-1]
    }else{
      stopifnot(str_detect(assay[i], "HAI"))
      dilutions <- original_dilutions
    }
    
    # If value differs from one the original dilutions by less than 1e-8, call it that dilution
    # this handles like 2^log(10, base =2), which R says is less than 10 by 1.776357e-15
    smallest_difference = min(abs(value - dilutions))
    if(smallest_difference < 1e-8){
      discrete_value <- dilutions[abs(value-dilutions) == smallest_difference]
      stopifnot(length(discrete_value) == 1)
    }else{
      if(value < min(dilutions)){
        discrete_value <- min(dilutions)
      }else{
        discrete_value <- max(dilutions[value >= dilutions])
      }
    }
    
    converted_titers <- c(converted_titers, discrete_value)
  }
  return(converted_titers)
}

# In this test, we assume the assay is HAI so that 3 gets coded as 5
stopifnot(censor_titers(c(3, 10, 15, 20, 35, Inf), assay = rep("HAI", 6)) == c(5, 10, 10, 20, 20, 1280))

# Here we assume FRNT so 3 is coded as 10
stopifnot(censor_titers(c(3, 10, 15, 20, 35, Inf), assay = rep("FRNT", 6)) == c(10, 10, 10, 20, 20, 1280))

stopifnot(censor_titers(2^log(40, base = 2), assay = "HAI") == 40)

annotate_years_with_vaccine_updates <- function(data){
  stopifnot((5 %in% data$year) == F) # Needs to be updated when we get to year 5 data
  
  data %>%
    mutate(vaccine_strain_updated = case_when(
      year == 2 & subtype %in% c('H3N2','H1N1') ~ T,
      year == 3 & subtype %in% c('H3N2','B/Victoria') ~ T,
      year == 4 & subtype == "H1N1" ~ T,
      T ~ F
    ))
}

count_prior_vax <- function(treatment){
  str_count(treatment, 'V') - 1 + (!str_detect(treatment, "V"))
}
stopifnot(count_prior_vax(c("P", "V", "PV", "VV", "PPP", "PPV", "VVV")) == c(0, 0, 0, 1, 0, 0, 2))

annotate_with_n_prior_vax <- function(data){
  data %>% 
    mutate(first_time_vax = str_count(treatment, 'V') == 1,
           n_previous_vax = count_prior_vax(treatment),
           # This variable is used to model a repeat vaccination effect only in seasons when the vaccine
           # was repeated from previous year
           n_previous_vax_if_repeated_vaccine = case_when(
             # In years when vaccine was updated, assign everyone 0 previous vaccinations
             vaccine_strain_updated ~ 0,
             # In years when vaccine was repeated, count number of past vaccinations given person's treatment
             !vaccine_strain_updated ~ n_previous_vax
           ))
}

annotate_with_timepoint_days <- function(data, timepoint_days){
  data %>% 
    left_join(timepoint_days) %>%
    select(pID, year, matches('treatment'), matches("timepoint"), matches("ndays"), everything())
}

annotate_with_previous_infections <- function(data, positive_pcr_tests, NAI_inferred_infections){

  infection_dates <- positive_pcr_tests %>% mutate(detection_method = "PCR")

  if(!is.null(NAI_inferred_infections)){
    infection_dates <- bind_rows(
      infection_dates,
      NAI_inferred_infections %>%
        select(-matches("in_low_circ_int")) %>%
        mutate(detection_method = "NAI")
    )
  }
  
  past_infection_summary_by_timepoint <- data %>%
    # Record date of the most recent infections prior 
    # to each timepoint, using exact specimen collection and timepoint dates
    select(pID, year, year_vdate, matches("timepoint")) %>%
    unique() %>%
    left_join(infection_dates, relationship = "many-to-many") %>%
    filter(!is.na(cdate),
           cdate < timepoint_date) %>%
    group_by(pID, year, timepoint, timepoint_date, year_vdate, result, detection_method) %>%
    summarise(most_recent_infection = max(cdate)) %>%
    mutate(result = str_remove(result, "_positive")) %>%
    # Create variables to indicate whether an infection happened before this 
    # year's vaccination date
    # And variables indicating whether an infection happened between the year's
    # vaccination date and the current timepoint.
    mutate(infection_since_year_vdate = most_recent_infection >= year_vdate,
           infection_before_year_vdate = most_recent_infection < year_vdate,
           # True for everyone in this tibble...
           infection_before_sample = T,
           # Recent infection = less than 365 days
           recent_infection_before_sample = as.numeric(
            timepoint_date - most_recent_infection,
            units = "days") < 365)
            
  # Left join results into the input data, separately for PCR and NAI in wide format
  PCR_wide <- past_infection_summary_by_timepoint %>%
    filter(detection_method == "PCR") %>%
    pivot_wider(names_from = result,
                values_from = matches("infection")) %>%
    rename_with(~paste0(., "_PCR"), matches("infection")) %>%
    select(-detection_method)

  annotated_data <- data %>%
    left_join(PCR_wide) 

  if(any(past_infection_summary_by_timepoint$detection_method == "NAI")){
    NAI_wide <- past_infection_summary_by_timepoint %>% 
        filter(detection_method == "NAI") %>%
        pivot_wider(names_from = result,
                    values_from = matches("infection")) %>%
        rename_with(~paste0(., "_NAI"), matches("infection")) %>%
        select(-detection_method)
    annotated_data <- annotated_data %>%
      left_join(NAI_wide) 
  }

  annotated_data <- annotated_data %>%
    # For all infection columns, replace NA with F
    mutate(across(.cols = matches('infection_since') |
                    matches('infection_before'),
                  .fns = ~replace_na(.x, F)))
  
  # If annotated_data doesn't include any infections of a particular
  # virus/subtype, the following loop will make that explicit (otherwise the
  # corresponding columns will simply be missing)

  for(dmethod in unique(past_infection_summary_by_timepoint$detection_method)){
    detection_results <- infection_dates %>%
      filter(detection_method == !!dmethod) %>%
      pull(result) %>%
      unique() %>%
      str_remove("_positive")

    for(result in detection_results){
      column_names <- paste0(c("most_recent_infection", 
                               "infection_since_year_vdate", 
                            "infection_before_year_vdate", 
                            "infection_before_sample", 
                            "recent_infection_before_sample"), "_", result, "_", dmethod)
      
      # For any columns that are missing from annotated data, add that column
      # populated with NAs (if it's most_recent_infection) or F (otherwise)
      for(column_name in column_names){
          if(! column_name %in% names(annotated_data)){
            if(str_detect(column_name, "most_recent_infection")){
              annotated_data[[column_name]] <- NA
            }else{
              annotated_data[[column_name]] <- F
            }
          }
        }
    }
  }

  # Add a special column indicating any prior infection with flu (by PCR only)
  # (used in specimen selection)
  annotated_data <- annotated_data %>%
    mutate(
      infection_before_sample_anyflu_PCR = 
        infection_before_sample_flu_A_PCR | infection_before_sample_flu_B_PCR |
        infection_before_sample_H3N2_PCR | infection_before_sample_H1N1_PCR)

  
  # If data has a subtype column (e.g., titer data), create infection summary
  # variables specific to the subtype in each row
  # For instance, for a titer measured against H3N2,
  # most_recent_infection_subtype_matched is the date of the most recent
  # infection with H3N2, though we keep most_recent_infection columnn for other
  # subtypes and SARS-CoV-2 in wide format
  if("subtype" %in% names(data)){

    # General function to pivot infection status columns to long format
    get_subtype_matched_columns <- function(annotated_data, data, column_name) {
      matched_cols <- grep(paste0("^", column_name), names(annotated_data), value = TRUE)

      # Exclude columns matching "flu_A", "flu_B", "cov2"
      matched_cols <- matched_cols[!str_detect(matched_cols, "flu_A|flu_B|cov2|anyflu")]

      annotated_data %>%
        select(pID, year, all_of(matched_cols), matches('timepoint')) %>%
        unique() %>%
        pivot_longer(
          cols = all_of(matched_cols),
          names_prefix = paste0(column_name, "_"),
          names_to = "subtype_and_method",
          values_to = paste0(column_name, "_subtype_matched")
        ) %>%
        mutate(method = case_when(
          str_detect(subtype_and_method, "PCR") ~ "PCR",
          str_detect(subtype_and_method, "NAI") ~ "NAI"
        )) %>%
        mutate(subtype = str_remove(subtype_and_method, paste0("_", method))) %>%
        select(-subtype_and_method) %>%
        pivot_wider(names_from = method, values_from = matches(column_name),
                    names_prefix = paste0(column_name, "_subtype_matched_")) %>%
        mutate(across(.cols = matches(column_name),
                  .fns = ~replace_na(.x, F))) %>%
        filter(subtype %in% unique(data$subtype)) %>%
        mutate(subtype = factor(subtype, levels = subtype_levels))
    }

    long_format_infection_before_vax <- get_subtype_matched_columns(annotated_data, data, "infection_before_year_vdate")
    long_format_infection_since_vax <- get_subtype_matched_columns(annotated_data, data, "infection_since_year_vdate")
    long_format_infection_before_sample <- get_subtype_matched_columns(annotated_data, data, "infection_before_sample")
    long_format_recent_infection_before_sample <- get_subtype_matched_columns(annotated_data, data, "recent_infection_before_sample")

    annotated_data <- annotated_data %>%
      left_join(long_format_infection_before_vax) %>%
      left_join(long_format_infection_since_vax) %>%
      left_join(long_format_infection_before_sample) %>%
      left_join(long_format_recent_infection_before_sample) %>%
      mutate(across(matches('subtype_matched'), ~replace_na(.x, F)))

  }

  return(annotated_data)
}

annotate_with_age_and_sex <- function(data, age_and_sex){
  left_join(data, age_and_sex) %>%
    mutate(age = age_in_year_1 + year - 1) %>%
    get_age_groups() %>%
    select(-age_in_year_1)
}

# Takes output of read_luminex keeps only vaccine
# strains, with names adjusted to match FRNT/HAI data
filter_luminex_to_vax_strains <- function(luminex){
  vax_strain_names_luminex <- c(
    "HI19rHA-A/Hawaii/70/2019",
    "VIC19rHA-A/Victoria/2570/2019",
    "MN19rHA-A/Minnesota/41/2019",
    "Tas20rHA-A/Tasmania/503/2020",
    "DA21rHA-A/Darwin/9/2021/H3N2",
    "WA19rHA-B/Washington/02/2019",
    "AT21rHA-B/Austria/1359417/2021",
    "PH13rHA-B/Phuket/3073/2013"
  )
  
  luminex %>%
    filter(strain %in% vax_strain_names_luminex) %>%
    mutate(strain = str_extract(strain, "[A-B]\\/[^\\s]*")) %>%
    mutate(strain = str_remove(strain, "\\/H3N2")) %>%
    mutate(strain = case_when(
      # These 2 strains are identical at aa level, so we'll code A/Victoria as A/Wisconsin
      # To match FNRT/HAI titer data
      strain == "A/Victoria/2570/2019" ~ "A/Wisconsin/588/2019",
      T ~ strain
    )) %>%
    flag_vaccine_strains() %>%
    filter(is_vaccine_strain) %>%
    select(-is_vaccine_strain)

}

annotate_with_cov2_vaccination <- function(data, cov2_vdates){

  most_recent_sc2_at_each_tp <- data %>%
    select(pID, year, timepoint, timepoint_date) %>%
    unique() %>%
    left_join(
      cov2_vdates %>%
        pivot_longer(cols = matches("cov2_vdate"), names_to = "dose", values_to = "sc2_vdate",
                     names_prefix = "cov2_vdate_d"),
        relationship = "many-to-many"
    ) %>%
    group_by(pID, year, timepoint, timepoint_date) %>%
    # Will return -Inf and warning if no prior SC vaccination. Suppress warning, replace -Inf with NA
    summarise(most_recent_sc2_vdate = suppressWarnings(max(sc2_vdate[sc2_vdate <= timepoint_date], na.rm = T))) %>%
    mutate(most_recent_sc2_vdate = case_when(
      is.infinite(most_recent_sc2_vdate) ~ NA_Date_,
      TRUE ~ most_recent_sc2_vdate
    )) %>%
    ungroup()

    annotated_data <- data %>%
      left_join(most_recent_sc2_at_each_tp)

    stopifnot(nrow(annotated_data) == nrow(data))
  return(annotated_data)

}

compute_fold_changes <- function(titers){
    titers %>%
      group_by(across(any_of(c('pID', 'year', 'strain')))) %>%
      mutate(fold_change = titer / titer[timepoint == 0]) %>%
      ungroup() %>%
      mutate(fold_change = case_when(
        timepoint == 0 ~ NA,
        timepoint != 0 ~ fold_change
      ))
}

get_fraction_by_group <- function(input_data, response_var){

  if(response_var == "fourfold_rise_or_greater"){
    input_data <- input_data %>% filter(timepoint != 0)
  }else{
    # Not strictly necessary, but review if using other variables...
    stopifnot(response_var == "titer_40_or_greater")
  }
  fraction_var_name <- paste0("fraction_", response_var)
  
  input_data %>%
    group_by(treatment, strain, year, timepoint, titer_type) %>%
    summarise(!!fraction_var_name := sum(.data[[response_var]])/n()) %>%
    ungroup()
    
}

# Plot "day 30" titers as a function of the actual interval from day 0
plot_d30_titers_vs_actual_interval <- function(titer_vaccine_responses, year){
  titer_vaccine_responses %>%
    filter(timepoint == 30) %>%
    filter(year == !!year) %>%
    ggplot(aes(x= ndays_since_year_vax,
               y = log(titer, base = 2))) +
    geom_point(aes(color = treatment)) +
    facet_wrap("subtype") +
    geom_smooth(se = T, alpha = 0.1,
                method = "lm",
                aes(color = treatment)) +
    scale_y_continuous(breaks = log(original_dilutions, base = 2),
                       labels = ~2^.x, limits = range(log(original_dilutions, base = 2))) +
    xlab("Actual days since vaccination") +
    ylab('"Day 30" titer') +
    theme(panel.border = element_rect(color = 'grey50'),
          legend.position = "top")
  
}

annotate_and_order_strains_and_subtypes <- function(data, strain_levels){
  if(is.null(strain_levels)){
    strain_levels <- unique(data$strain)
  }
  data %>%
    mutate(strain = factor(strain, levels = strain_levels)) %>%
    mutate(subtype = subtype_labeller(strain),
           subtype = factor(subtype, levels = subtype_levels))
}

annotate_with_treatment <- function(data, treatment_assignment){
  left_join(data,
            treatment_assignment %>% select(pID, drive, matches("treatment_")) %>%
              pivot_longer(cols = matches('treatment_'), names_to = "year",
                           values_to = "treatment") %>%
              mutate(year = as.integer(str_remove(year, "treatment_Y")))) %>%
    select(pID, year, treatment, everything())

}

# When fitting models, selects which batch to use for each person/timepoint/strain
select_titer_batches <- function(titers){
  
  # For each titer type (HAI vs FRNT), for each person/year/strain, if there are multiple
  # measurements for the same timepoint (multiple batches), keep the batch the
  # that includes the same person for the most years.
  
  selected_titers <- titers %>%
    group_by(pID, titer_type, batch) %>%
    mutate(n_years_in_batch = length(unique(year))) %>%
    group_by(pID, year, timepoint, strain) %>%
    filter(n_years_in_batch == max(n_years_in_batch)) %>%
    ungroup()
  
  # Check that we retained a single measurement per individual/year/timepoint/
  # strain/titer type (HAI vs FRNT)
  stopifnot(
    selected_titers %>%
      group_by(pID, year, timepoint, strain, titer_type) %>%
      count() %>%
      pull(n) %>%
      unique() == 1
  )
  
  return(selected_titers)
  
}

# Runs pairwise Wilcoxon tests for comparing absolute post-vaccination titers/luminex
run_pairwise_wilcoxon <- function(response_data, response_var, p.adjust.method){

  base_function <- function(response_data, strain, year, timepoint){
    data_subset <- response_data %>%
      filter(year == !!year, timepoint == !!timepoint, strain == !!strain)

    if(nrow(data_subset) > 0 && length(unique(data_subset$treatment)) > 1) {

      # Compute all pairwise tests
      pw_wilcox_test <- pairwise.wilcox.test(x = data_subset[[response_var]], g = data_subset$treatment, alternative = 'two.sided',
                                             p.adjust.method = p.adjust.method)
    
      pw_wilcox_test <- as_tibble(pw_wilcox_test[[3]], rownames = 'treatment_1')
    
      pw_wilcox_test <- pw_wilcox_test %>%
        pivot_longer(cols = colnames(pw_wilcox_test)[colnames(pw_wilcox_test) != 'treatment_1'],
                     names_to = 'treatment_2', values_to = 'p_value') %>%
        mutate(strain = strain, year = year, timepoint = timepoint) %>%
        # These filters keep only the lower part of the comparison matrix
        filter(treatment_1 != treatment_2) %>%
        filter(!is.na(p_value)) %>%
        select(year, timepoint, strain, everything())  

      # Check that the correct number of pairwise comparisons was made
      expected_comparisons <- choose(length(unique(data_subset$treatment)), 2)
      stopifnot(nrow(pw_wilcox_test) == expected_comparisons)

      # Annotate table with GMT of each group for each pairwise comparison
      # Table of GMTs by group/year/timepoints
      if(response_var == "log2_titer"){
        GMT_table <- data_subset %>%
          group_by(year, timepoint, strain, treatment) %>%
          summarise(mean_log2_titer = mean(log2_titer)) %>%
          mutate(GMT = 2^mean_log2_titer) %>%
          select(-mean_log2_titer)
      
        pw_wilcox_test <- left_join(
          pw_wilcox_test, GMT_table,
          by = c("year", "timepoint", "strain", "treatment_1" = "treatment")) %>%
          rename(GMT_treatment_1 = GMT) %>%
          left_join(GMT_table,
          by = c("year", "timepoint", "strain", "treatment_2" = "treatment")) %>%
          rename(GMT_treatment_2 = GMT) 
      }
    }else{
      pw_wilcox_test <- NULL
    }
    return(pw_wilcox_test)
  }
  strain_year_timepoint_combinations <- response_data %>% select(strain, year, timepoint) %>% unique()
  
  pw_wilcox_results <- bind_rows(mapply(base_function,
                                        strain = strain_year_timepoint_combinations$strain,
                                        year = strain_year_timepoint_combinations$year,
                                        timepoint  = strain_year_timepoint_combinations$timepoint,
                                        MoreArgs = list(response_data = response_data), SIMPLIFY = F))
  
  return(pw_wilcox_results)

}

# To compare fraction fourfold rises or fraction of people at or above 40
# Very similar to run_pairwise_wilcoxon, but not worth generalizing.
run_pairwise_fisher_tests <- function(titers, response_var, p.adjust.method = "holm"){

  input_data <- titers

  if(response_var == "fourfold_rise_or_greater"){
    input_data <- titers %>% filter(timepoint != 0)
  }

  base_function <- function(input_data, strain, year, timepoint){
    data_subset <- input_data %>%
      filter(year == !!year, timepoint == !!timepoint, strain == !!strain)

    pw_fisher_test <- pairwise.fisher.test(x = data_subset[[response_var]],
                                           g = data_subset$treatment,
                                           alternative = 'two.sided',
                                           p.adjust.method = p.adjust.method,
                                           workspace = 2000000)
    
    pw_fisher_test <- as_tibble(pw_fisher_test$p.value, rownames = 'treatment_1')
    
    pw_fisher_test %>%
      pivot_longer(cols = colnames(pw_fisher_test)[colnames(pw_fisher_test) != 'treatment_1'],
                   names_to = 'treatment_2', values_to = 'p_value') %>%
      mutate(strain = strain, year = year, timepoint = timepoint) %>%
      select(year, timepoint, strain, everything())
    
  }
  strain_year_timepoint_combinations <- input_data %>% select(strain, year, timepoint) %>% unique()
  
  pw_fisher_results <- bind_rows(mapply(base_function,
                                        strain = strain_year_timepoint_combinations$strain,
                                        year = strain_year_timepoint_combinations$year,
                                        timepoint  = strain_year_timepoint_combinations$timepoint,
                                        MoreArgs = list(input_data = input_data),
                                        SIMPLIFY = F))
  return(pw_fisher_results)
  
}

plot_postvax_response <- function(response_data, response_var,
                                  year, vaccine_strains_only = T,
                                  prior_infection_timing = "before_sample",
                                  prior_infection_virus = "subtype_matched",
                                  include_NAI_infections = F,
                                  include_d0 = T,
                                  p.adjust.method = "holm",
                                  show_points = T,
                                  show_pairwise_comparisons = T,
                                  transpose = F){

  if(str_detect(response_var, "multiplex")){
    stopifnot(length(unique(response_data$measurement_type)) == 1)
  }

  if(transpose == T){
    facet_control <- facet_grid(timepoint ~ subtype)
  }else{
    facet_control <- facet_grid(subtype ~ timepoint)
  }
  
  # Run pairwise Wilcoxon tests internally
  pw_wilcox_test <- run_pairwise_wilcoxon(response_data = response_data, response_var = response_var, p.adjust.method = p.adjust.method)
  
  highlight_var <- paste("infection", prior_infection_timing,
                         prior_infection_virus, sep = "_")

  highlight_title <- paste0(
    recode_values(prior_infection_virus,
               "H3N2" ~ "H3N2",
               "H1N1" ~ "H1N1",
               "B/Victoria" ~ "B/Victoria",
               "B/Yamagata" ~ "B/Yamagata",
               "subtype_matched" ~ "Subtype-matched",
               "anyflu" ~ "Any influenza",
               "flu_A" ~ "Influenza A",
               "flu_B" ~ "Influenza B",
               "cov2" ~ "SARS-CoV-2"),
    " infection ",
    recode_values(prior_infection_timing,
               "before_year_vdate" ~ "before\nmost recent vaccination",
               "since_year_vdate" ~ "since\nmost recent vaccination",
               "before_sample" ~ "before current time point")
  )
  
  plot_data <- response_data %>%
    filter(year == !!year)

  if(include_NAI_infections){
    if(prior_infection_virus == "cov2"){
      stop("Can't highlight SC2 infections detected by NAI")
    }
    plot_data <- plot_data %>%
      mutate(highlight_var = !!sym(paste0(highlight_var, "_PCR")) | !!sym(paste0(highlight_var, "_NAI")))
  }else{
    plot_data <- plot_data %>%
      mutate(highlight_var = !!sym(paste0(highlight_var, "_PCR")))
  }
  

  if(!include_d0){
    plot_data <- plot_data %>%
      filter(timepoint != 0)
  }

  if(nrow(plot_data) == 0){
    return(NULL)
  }else{
    if(vaccine_strains_only){
      plot_data <- plot_data %>%
      flag_vaccine_strains() %>%
      filter(is_vaccine_strain)
    }

    if(response_var == "multiplex"){
      ylabel <- paste0("Multiplex – ", unique(response_data$measurement_type))
    }else{
          ylabel <- recode_values(response_var,
                         "log2_titer" ~ "Titer",
                         "log10_luminex" ~ "Luminex reading (log10)",
                         "luminex" ~ "Luminex reading",
                         "subclass_proportion" ~ "Luminex subclass proportion")
    }

    significance_symbol_shift <- recode_values(response_var,
                                            "log2_titer" ~ 0.5,
                                            "log10_luminex" ~ 0.1,
                                            "luminex" ~ (max(plot_data[[response_var]]) -
                                            min(plot_data[[response_var]]))  / 20,
                                            "subclass_proportion" ~ 0.01,
                                            "multiplex" ~
                                            (max(plot_data[[response_var]]) -
                                            min(plot_data[[response_var]]))  / 20)
  
    relabel_titer_type <- function(data, response_var){
    
      if("titer_type" %in% names(data)){
        data %>%
          mutate(titer_type = case_when(
            str_detect(titer_type, 'FRNT') ~ 'FRNT titers',
            str_detect(titer_type, 'HAI') ~ 'HAI titers'
          ))
      }else{
        data %>%
          mutate(titer_type = str_to_title(str_remove(response_var, "log10_")))
      }
    }
  
    base_pl <- plot_data %>%
      relabel_titer_type(response_var = response_var) %>%
      mutate(treatment = factor(treatment, levels = treatment_levels)) %>%
      mutate(timepoint = relabel_timepoint(timepoint)) %>%
      ggplot(aes(x = treatment, y = .data[[response_var]])) +
      geom_boxplot(outlier.alpha = 0, aes(fill = treatment), box.linewidth = boxplot_line_width, median.linewidth =  boxplot_median_width)

    if(show_points){
      base_pl <- base_pl +
        geom_point(
          position = position_jitter(height = 0, width = 0.2),
          aes(size = highlight_var, color = highlight_var,
              alpha = highlight_var),
          stroke = boxplot_point_stroke,
          shape = boxplot_point_shape)
    }

    base_pl <- base_pl + 
      #scale_shape_manual(values = c(16,4)) +
      scale_size_manual(name = highlight_title,
                        values = c(boxplot_point_size, 1.5))+
      scale_color_manual(name = highlight_title,
                         values = c('black', 'red')) +
      scale_alpha_manual(name = highlight_title,
                         values = c(boxplot_point_alpha, 1)) +
      facet_control + 
      xlab('Intervention') + 
      ylab(ylabel) +
      baseline_figure_settings +
      theme(legend.position = 'top') +
      guides(fill = "none")
  
    base_pl <- base_pl +
      scale_fill_manual(values = treatment_colors)
  
    if(response_var == 'log2_titer'){
      y_axis_breaks = log(original_dilutions, base = 2)
      base_pl <- base_pl + 
        scale_y_continuous(labels = function(x){2^x}, breaks = y_axis_breaks)
    }
  
    if (!is.null(pw_wilcox_test)) {
      p_value_data <- pw_wilcox_test %>%
        p_value_labeller() %>%
        filter(year == !!year)
        
      if(!include_d0){
        p_value_data <- p_value_data %>%
        filter(timepoint != 0)
      }  
        
      if (vaccine_strains_only) {
        p_value_data <- p_value_data %>%
          flag_vaccine_strains() %>%
          filter(is_vaccine_strain)
        
        pw_wilcox_test <- pw_wilcox_test
      }
    
      # Only need to add bars if there were significant tests 
      if (nrow(p_value_data) > 0 & show_pairwise_comparisons) {
        pw_heights <- p_value_data %>%
          select(treatment_1, treatment_2) %>%
          unique() %>%
          mutate(n = 1:n()) %>%
          mutate(height = max(response_data[[response_var]]) + significance_symbol_shift * n)
      
        p_value_data <- left_join(p_value_data, pw_heights) %>%
          mutate(timepoint = relabel_timepoint(timepoint)) 
      
        p_value_data <- left_join(p_value_data, response_data %>% select(subtype, strain, any_of("titer_type")) %>%
                                    unique()) %>%
          relabel_titer_type(response_var = response_var)
      
        pl <- base_pl +
          geom_segment(data = p_value_data,
                       aes(y = height, yend = height, x = treatment_1, xend = treatment_2,
                           linetype = comparison_with_placebo), show.legend = F) +
          geom_text(data = p_value_data,
                    aes(x = treatment_2, y = height, label = plot_label),
                    position = position_nudge(x = 0.1, y = 0.05),
                    size = default_figure_font_size, size.unit = "pt")
      } else {
        pl <- base_pl
      }
    } else {
      pl <- base_pl
    }
  return(pl)
  }
}

# Plots fractions by group
# response_var = "fourfold_rise_or_greater" or "titer_40_or_greater"
plot_fraction_by_group <- function(input_data, response_var, year,
                                   vaccine_strains_only = T){
  
  if(vaccine_strains_only){
    input_data <- input_data %>%
      flag_vaccine_strains() %>%
      filter(is_vaccine_strain)
  }
  
  
  if(response_var == "fourfold_rise_or_greater"){
    ylabel <- expression("Fraction of participants with " >= "4-fold rise")
  }else{
    stopifnot(response_var == "titer_40_or_greater")
    ylabel <- expression("Fraction of participants with titers " >= "40")
  }
  
  plot_data <- get_fraction_by_group(input_data %>% filter(timepoint != 0, timepoint != 365),
                                     response_var)
  fraction_var_name <- names(plot_data %>% select(matches("fraction_")))
  
  stopifnot(length(fraction_var_name) == 1)
  
  plot_data <- plot_data %>%
    annotate_and_order_strains_and_subtypes(strain_levels) %>%
    filter(year == !!year)
  
  pw_fisher_tests <- run_pairwise_fisher_tests(input_data,
                                               response_var = response_var)
  
  pw_fisher_tests <- pw_fisher_tests %>%
    semi_join(plot_data %>%
                select(year, timepoint, strain) %>%
                unique())

  p_value_data <- pw_fisher_tests %>%
    p_value_labeller() %>%
    annotate_and_order_strains_and_subtypes(strain_levels) %>%
    mutate(timepoint = relabel_timepoint(timepoint)) %>%
    mutate(comparison_with_placebo = factor(comparison_with_placebo, levels = c(F, T)))
  
  base_pl <- plot_data %>%
    mutate(treatment = factor(treatment, levels = treatment_levels)) %>%
    mutate(timepoint = relabel_timepoint(timepoint)) %>%
    ggplot(aes(x = treatment, y = .data[[fraction_var_name]])) +
    geom_col(aes(fill = treatment), color = 'gray40', linewidth = boxplot_line_width) +
    facet_grid(subtype ~ timepoint) +
    labs(
      y = ylabel, 
      x = "Intervention"
    ) +
    baseline_figure_settings +
    theme(legend.position = 'None') +
    scale_fill_manual(values = treatment_colors)

  if(nrow(p_value_data) > 0){
    significance_symbol_shift <- 0.04
  
      pw_heights <- p_value_data %>%
        select(treatment_1, treatment_2) %>%
        unique() %>%
        mutate(n = 1:n()) %>%
        mutate(height = max(plot_data[[fraction_var_name]]) + significance_symbol_shift * n)
  
      p_value_data <- p_value_data %>% left_join(pw_heights)
    
    if(all(c(T, F) %in% unique(p_value_data$comparison_with_placebo))){
      linetype_values <- c(1, 2)
    }else{
      if(all(p_value_data$comparison_with_placebo == T)){
        linetype_values <- 2
      }
    }

    pl <- base_pl +
      geom_segment(data = p_value_data,
                   aes(y = height, yend = height, x = treatment_1, xend = treatment_2,
                       linetype = comparison_with_placebo)) +
      geom_text(data = p_value_data,
                aes(x = treatment_2, y = height, label = plot_label),
                position = position_nudge(x = -0.2, y = -0.01),
                size = default_figure_font_size, size.unit = "pt") +
      scale_linetype_manual(values = linetype_values)
  }else{
    pl <- base_pl
  }
  return(pl)

}

plot_luminex_vs_titer_correlations <- function(titers, luminex){
  joint_data <- left_join(luminex, titers) %>%
    mutate(timepoint = factor(paste0('Day ', timepoint), levels = paste0('Day ', sort(unique(titers$timepoint)))))
    
  
  correlation_coefficient <- joint_data %>%
    group_by(timepoint, subtype) %>%
    mutate(n_obs = sum(!is.na(log2_titer)&!is.na(log10_luminex))) %>%
    filter(n_obs > 0) %>%
    summarise(cor_coef = cor.test(log2_titer, log10_luminex)$estimate,
              cor_CI = paste(round(cor.test(log2_titer, log10_luminex)$conf.int,2), collapse = '-')) %>%
    ungroup() %>%
    mutate(label = paste0("r = ",round(cor_coef,2), "(", cor_CI, ")"))
  
    
  joint_data %>%
    ggplot(aes(x = log2_titer, y = log10_luminex)) +
    geom_point() +
    facet_grid(subtype~timepoint, scales = 'free_y') +
    geom_text(data = correlation_coefficient, aes(x = 4.5, y = 4.8, label = label),
              size = 3) +
    #geom_smooth(method = 'lm') +
    xlab("Titer\n(FRNT for H3N2, HAI for the others)") +
    ylab("Luminex reading (log10)") +
    scale_x_continuous(breaks = log(original_dilutions, base = 2),
                       labels = function(x){2^x}) +
    theme(panel.border = element_rect(colour = 'gray90'),
          axis.text.x = element_text(angle = -45, vjust = 0))
    
  
}


# Similar to plot_postvax_response, but makes main text figure showing years 3 and 4 together
plot_years_3_4_postvax <- function(response_data, response_var, p.adjust.method = "holm"){

  # Looking only at vaccine strains
  response_data <- response_data %>%
    flag_vaccine_strains() %>%
    filter(is_vaccine_strain)
  
  input_data <-  response_data %>%
    filter(year == 3 & timepoint %in% c(0, 30, 182, 365) |
           year == 4 & timepoint %in% c(0, 30, 182)) %>%
    mutate(n_prior_vax = case_when(
      !str_detect(treatment, "V") ~ "Placebo",
      str_detect(treatment, "V") ~ paste0(ifelse(count_prior_vax(treatment) > 0,
                                                 count_prior_vax(treatment), "No"),
                                          ifelse(count_prior_vax(treatment) > 1,
                                                 " prior vaccinations", 
                                                 " prior vaccination"))
    )) %>%
    mutate(n_prior_vax = factor(n_prior_vax,
           levels = c("Placebo", "No prior vaccination", "1 prior vaccination",
                      "2 prior vaccinations", "3 prior vaccinations")))

    colors <- input_data %>%
      select(treatment, n_prior_vax) %>% 
      unique() %>%
      left_join(tibble(treatment = names(treatment_colors),
                       color = treatment_colors))

    color_vector <- colors$color
    names(color_vector) <- colors$n_prior_vax

    pw_wilcox_test <- run_pairwise_wilcoxon(response_data = input_data, response_var = response_var, p.adjust.method = p.adjust.method)
    p_value_annotation <- p_value_labeller(pw_wilcox_test)

    medians <- input_data %>%
      group_by(treatment, year, timepoint, subtype, n_prior_vax) %>%
      summarise(median = median(.data[[response_var]]))

    base_pl <- input_data %>% 
      ggplot(aes(x = factor(timepoint), y = .data[[response_var]], fill = factor(n_prior_vax))) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.alpha = 0,
                   box.linewidth = boxplot_line_width, median.linewidth =  boxplot_median_width) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
            alpha = main_text_alpha, size = main_text_point_size, stroke = boxplot_point_stroke,
            shape = boxplot_point_shape, show.legend = FALSE) +
      xlab("Days since vaccination") +
      ylab(response_var) +
      labs(fill = "") +
      baseline_figure_settings + 
      theme(legend.position = "top", legend.justification = "left",
            legend.direction = "horizontal", legend.text = element_text(size = small_font_size)) +
      scale_fill_manual(values = color_vector) +
      #geom_label(data = medians, aes(y = median, label = treatment),
      #     position = position_dodge(width = 0.8), size = 2, show.legend = F) +
      facet_grid(subtype ~ year, scales = "free_x", 
         labeller = labeller(year = function(x) paste0("Year ", x))) +
      scale_y_continuous(labels = function(x){2^x},
                         limits = c(NA, log2(2560)),
                       breaks = log(original_dilutions[c(1, 3, 5, 7, 9)], base = 2))
  
  return(base_pl)
}

plot_percent_with_infections_relative_to_vdate <- function(timepoint_days, treatment_assignment,  positive_pcr_tests, NAI_inferred_infections,
                                                           infection_time = "since_year_vdate",
                                                           p.adjust.method){

  valid_infection_times <- c("since_year_vdate", "before_year_vdate")
  if (!infection_time %in% valid_infection_times) {
    stop(paste0("`infection_time` must be one of: ", paste(valid_infection_times, collapse = ", ")))
  }

  col_prefix <- paste0("infection_", infection_time)
  pcr_col    <- paste0(col_prefix, "_PCR")
  nai_col    <- paste0(col_prefix, "_NAI")
  any_col    <- paste0(col_prefix, "_any_method")

  infection_since_vdate_data <- timepoint_days %>%
    filter(timepoint == 365, !is.na(timepoint_date)) %>%
    annotate_with_treatment(treatment_assignment) %>%
    filter(drive == 1) %>%
    annotate_with_previous_infections(positive_pcr_tests = positive_pcr_tests, NAI_inferred_infections) %>%
    mutate(
      !!pcr_col := .data[[paste0(col_prefix, "_flu_A_PCR")]] | .data[[paste0(col_prefix, "_H1N1_PCR")]] |
        .data[[paste0(col_prefix, "_H3N2_PCR")]] | .data[[paste0(col_prefix, "_flu_B_PCR")]] | .data[[paste0(col_prefix, "_B/Victoria_PCR")]],
      !!nai_col := .data[[paste0(col_prefix, "_H3N2_NAI")]] | .data[[paste0(col_prefix, "_H1N1_NAI")]]
    ) %>%
    mutate(!!any_col := .data[[pcr_col]] | .data[[nai_col]]) %>%
    select(pID, year, treatment, timepoint, all_of(c(pcr_col, nai_col, any_col))) %>%
    pivot_longer(cols = starts_with(col_prefix),
                 names_to = "method",
                 values_to = col_prefix,
                 names_prefix = paste0(col_prefix, "_")) %>%
    mutate(method = case_when(
      method == "any_method" ~ "PCR or NAI",
      T ~ method
    ))


  # This function was writen for titer data but can be used here by renaming "method" as "strain"
  pw_fisher_tests_since_vdate <- run_pairwise_fisher_tests(
    titers = infection_since_vdate_data %>% rename(strain = method),
    response_var = col_prefix,
    p.adjust.method = p.adjust.method) %>%
    rename(method = strain) %>%
    # When both treatments had 0 infections, p_value will be NA
    filter(!is.na(p_value)) %>%
    p_value_labeller()

  ylabel <- if (infection_time == "before_year_vdate") "Fraction infected BEFORE vaccination" else "Fraction infected AFTER vaccination"


  fractions_by_group <- infection_since_vdate_data %>%
    group_by(year, timepoint, treatment, method, .data[[col_prefix]]) %>%
    count() %>%
    group_by(year, timepoint, treatment, method) %>%
    mutate(total = sum(n),
           fraction = n / total,
           label = paste0(round(fraction * 100, 1), "% (", n, "/", total, ")")) %>%
    ungroup()  %>%
    filter(.data[[col_prefix]])

  since_vdate_pl <- fractions_by_group %>%
    ggplot(aes(x = treatment, y = fraction)) +
    facet_grid(method ~ year, scales = "free_x") +
    geom_col() +
    geom_text(aes(label = label), vjust = -0.5, size = 3)  +
    xlab("Intervention") +
    ylab(ylabel) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Method\n", breaks = NULL)) +
    scale_x_discrete(sec.axis = sec_axis(transform = identity, name = "Year\n", breaks = NULL))

  if(nrow(pw_fisher_tests_since_vdate) > 0){
    pw_heights <- pw_fisher_tests_since_vdate %>%
      select(treatment_1, treatment_2) %>%
      unique() %>%
      mutate(n = 1:n()) %>%
      mutate(height = max(fractions_by_group$fraction) + 0.02 + 0.02 * n)

    p_value_annotation <- pw_fisher_tests_since_vdate %>%
      left_join(pw_heights)
      
    since_vdate_pl <- since_vdate_pl +
      geom_segment(data = p_value_annotation,
                   aes(y = height, yend = height, x = treatment_1, xend = treatment_2),
                   position = position_nudge(x = 0, y = -0.01)) +
      geom_text(data = p_value_annotation,
                aes(x = treatment_2, y = height, label = plot_label),
                position = position_nudge(x = -0.2, y = -0.01))
  }

  return(since_vdate_pl)

}

plot_prevaccination_titer_correlations <- function(titers){

  # For the same person, compute correlation in pre-vaccination titers across subtypes 
  prevax_titers <- titers %>%
    flag_vaccine_strains() %>%
    filter(is_vaccine_strain) %>%
    # Look only at pre-vaccination titers for vaccinated people
    filter(str_detect(treatment, "V")) %>%
    filter(timepoint == 0) %>%
    select(pID, year, treatment, subtype, log2_titer) %>%
    mutate(n_previous_vax = count_prior_vax(treatment)) %>%
    select(-treatment)

  wide_format_tibble <- prevax_titers %>%
    pivot_wider(names_from = subtype, values_from = log2_titer)

  # Build every unordered pair of subtypes that is actually present in the data
  present_subtypes <- intersect(subtype_levels, names(wide_format_tibble))
  subtype_pairs <- combn(present_subtypes, 2, simplify = FALSE)

  # Reshape to one row per (person, year, subtype pair): x = first subtype titer,
  # y = second subtype titer. This lets a single facet_grid handle all pairs.
  pairwise_long <- map_dfr(subtype_pairs, function(pair) {
    wide_format_tibble %>%
      transmute(
        pID, year, n_previous_vax,
        x_subtype = pair[1],
        y_subtype = pair[2],
        pair_label = paste0(pair[1], " vs ", pair[2]),
        x = .data[[pair[1]]],
        y = .data[[pair[2]]]
      )
  }) %>%
    filter(!is.na(x), !is.na(y)) %>%
    mutate(pair_label = factor(pair_label,
                               levels = map_chr(subtype_pairs,
                                                ~ paste0(.x[1], " vs ", .x[2]))))

  # Check that, for each individual and year, the correct number of rows is present
  stopifnot(pairwise_long %>% group_by(pID, year) %>% count() %>% pull(n) %>% unique() == length(subtype_pairs))

  scatter_pl <- ggplot(pairwise_long, aes(x = x, y = y, color = factor(n_previous_vax))) +
    geom_jitter(width = 0.15, height = 0.15,
                alpha = boxplot_point_alpha, size = 0.6, shape = boxplot_point_shape, stroke = boxplot_point_stroke) +
    geom_smooth(method = "lm", se = F, linewidth = 0.7) +
    facet_grid(pair_label ~ year, labeller = labeller(year = function(y) paste0("Year ", y))) +
    scale_color_manual(values = year_4_nprior_vax_colors) +
    labs(x = "Pre-vaccination titer (subtype 1)",
         y = "Pre-vaccination titer (subtype 2)",
         color = "N. previous vaccinations") +
    baseline_figure_settings +
    theme(panel.border = element_rect(color = "gray90", fill = NA))

  # Pooled scatter: all subtype pairs and years together, one panel per n_previous_vax 
  # One pooled Pearson correlation (95% CI) per panel, placed in the top-left corner.
  pooled_cor_labels <- pairwise_long %>%
    group_by(n_previous_vax) %>%
    summarise(
      ct    = list(cor.test(x, y)),
      label = sprintf("r = %.2f (%.2f, %.2f)",
                      ct[[1]]$estimate, ct[[1]]$conf.int[1], ct[[1]]$conf.int[2]),
      x     = min(pairwise_long$x, na.rm = TRUE),
      y     = max(pairwise_long$y, na.rm = TRUE),
      .groups = "drop"
    )

  pooled_scatter_pl <- ggplot(pairwise_long, aes(x = x, y = y, color = factor(n_previous_vax))) +
    geom_jitter(width = 0.15, height = 0.15,
                alpha = boxplot_point_alpha, size = 0.6,
                shape = boxplot_point_shape, stroke = boxplot_point_stroke) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
    geom_text(data = pooled_cor_labels, aes(label = label),
              hjust = 0, vjust = 1, size = small_font_size, size.unit = "pt",
              color = "black") +
    facet_wrap(~ n_previous_vax, nrow = 1,
               labeller = labeller(n_previous_vax = function(n) paste0(n, " prior vax."))) +
    scale_color_manual(values = year_4_nprior_vax_colors, guide = "none") +
    labs(x = "Pre-vaccination titer (subtype 1)",
         y = "Pre-vaccination titer (subtype 2)") +
    baseline_figure_settings +
    theme(panel.border = element_rect(color = "gray90", fill = NA))

  # Correlation heatmap: Pearson r (95% CI) per subtype pair
  # Compute Pearson correlation + 95% CI within each (year, n_previous_vax, pair).
  # cor.test needs >= 3 finite pairs to return a CI, so guard small strata.
  pairwise_cor <- pairwise_long %>%
    group_by(year, n_previous_vax, x_subtype, y_subtype) %>%
    summarise(
      n = n(),
      ct = if (n() >= 3 && sd(x) > 0 && sd(y) > 0) list(cor.test(x, y)) else list(NULL),
      .groups = "drop"
    ) %>%
    mutate(
      r  = map_dbl(ct, ~ if (is.null(.x)) NA_real_ else unname(.x$estimate)),
      lo = map_dbl(ct, ~ if (is.null(.x)) NA_real_ else .x$conf.int[1]),
      hi = map_dbl(ct, ~ if (is.null(.x)) NA_real_ else .x$conf.int[2]),
      p  = map_dbl(ct, ~ if (is.null(.x)) NA_real_ else .x$p.value)
    ) %>%
    select(-ct) %>%
    # Keep only significant correlations (no multiple-testing correction); blank the rest
    mutate(
      significant = !is.na(p) & p < 0.05,
      r     = if_else(significant, r, NA_real_),
      label = if_else(significant, sprintf("%.2f\n(%.2f, %.2f)", r, lo, hi), NA_character_)
    )

  # Show only the lower triangle: each subtype pair appears once. combn() orders
  # x_subtype before y_subtype in present_subtypes, so with x normal and y reversed
  # these land below the diagonal.
  heatmap_data <- pairwise_cor %>%
    mutate(
      x_subtype = factor(x_subtype, levels = present_subtypes),
      y_subtype = factor(y_subtype, levels = rev(present_subtypes))
    )

  heatmap_pl <- ggplot(heatmap_data, aes(x = x_subtype, y = y_subtype, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 2) +
    facet_grid(n_previous_vax ~ year,
               labeller = labeller(year = function(y) paste0("Year ", y),
                                   n_previous_vax = function(n) paste0(n, " prior vax."))) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         midpoint = 0, limits = c(-1, 1), na.value = "gray90") +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    baseline_figure_settings +
    # Square tiles via aspect.ratio (not coord_equal, which pads panels with margin)
    theme(panel.border = element_rect(color = "gray90", fill = NA),
          axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1)

  return(list(scatter = scatter_pl, pooled_scatter = pooled_scatter_pl, heatmap = heatmap_pl))

}


# ==== Other objects ====

treatment_levels <- c("P", "V", "PP", "PV", "VV",
                      "PPP", "PPV", "PVV", "VVV",
                      "PPPP", "PPPV", "PPVV", "PVVV", "VVVV",
                      "PPPPV", "PPPVV", "PPVVV", "PVVVV", "VVVVV")

subtype_levels <- c('H1N1', 'H3N2', 'B/Victoria', 'B/Yamagata')

original_dilutions <- c(5, 10, 20, 40, 80, 160, 320, 640, 1280) # Original dilutions, including the LOD for HAI (5) and FRNT (10)

possible_pairwise_GMTs <- expand_grid(value1 = original_dilutions, value2 = original_dilutions) %>%
  mutate(GMT = sqrt(value1 * value2)) %>%
  select(GMT) %>%
  unique() %>%
  pull()

vaccine_strains <- tibble(
  year = rep(1:4, times = c(6, 6, 4, 5)),
  strain = c(
    "A/Hong Kong/45/2019", "A/Minnesota/41/2019", "A/Hawaii/66/2019", 
    "A/Hawaii/70/2019", "B/Phuket/3073/2013", "B/Washington/02/2019",
    "A/Cambodia/e0826360/2020", "A/Tasmania/503/2020", "A/Wisconsin/588/2019", 
    "A/Victoria/2570/2019", "B/Phuket/3073/2013", "B/Washington/02/2019",
    "A/Darwin/9/2021", "A/Wisconsin/588/2019", "B/Phuket/3073/2013", 
    "B/Austria/1359417/2021",
    "A/Darwin/9/2021", "A/Darwin/6/2021", "A/Wisconsin/67/2022", "B/Phuket/3073/2013", 
    "B/Austria/1359417/2021"
  )
) %>%
  mutate(subtype = subtype_labeller(strain))


  # For each year, list prior-current vaccine strain pair in long format
  # (i.e., for each year/prior year combination, indicate the pair of strains (for each subtype) for that pair of years (strain_pair)
  vax_strain_pairs <- tibble(year = 2:4) %>%
    rowwise() %>%
    mutate(previous_year = list(1:(year - 1))) %>%
    unnest(previous_year) %>%
    expand_grid(subtype = subtype_levels) %>%
    arrange(desc(subtype)) %>%
    mutate(strain_pair = c(
      "A/Hong Kong/45/2019-A/Cambodia/e0826360/2020",
      "A/Hong Kong/45/2019-A/Darwin/6/2021",
      "A/Cambodia/e0826360/2020-A/Darwin/6/2021",
      "A/Hong Kong/45/2019-A/Darwin/6/2021",
      "A/Cambodia/e0826360/2020-A/Darwin/6/2021",
      "A/Darwin/6/2021-A/Darwin/6/2021",
      "A/Hawaii/70/2019-A/Wisconsin/588/2019",
      "A/Hawaii/70/2019-A/Wisconsin/588/2019",
      "A/Wisconsin/588/2019-A/Wisconsin/588/2019",
      "A/Hawaii/70/2019-A/Wisconsin/67/2022",
      "A/Wisconsin/588/2019-A/Wisconsin/67/2022",
      "A/Wisconsin/588/2019-A/Wisconsin/67/2022",
      rep("B/Phuket/3073/2013-B/Phuket/3073/2013", 6),
      "B/Washington/02/2019-B/Washington/02/2019",
      "B/Washington/02/2019-B/Austria/1359417/2021",
      "B/Washington/02/2019-B/Austria/1359417/2021",
      "B/Washington/02/2019-B/Austria/1359417/2021",
      "B/Washington/02/2019-B/Austria/1359417/2021",
      "B/Austria/1359417/2021-B/Austria/1359417/2021"
      )) %>%
    # Create a label for each pair
    separate(strain_pair, sep = "-", remove = F, into = c("strain1", "strain2")) %>%
    mutate(across(c(strain1, strain2), ~ as.character(strain_labeller(.)))) %>%
    mutate(across(c(strain1, strain2), ~ str_remove_all(., "\\s*\\([^\\)]+\\)"))) %>%
    mutate(pair_label = paste0(strain1, "-", strain2)) %>%
    select(-strain1, -strain2)


treatment_colors <- tibble(treatment = treatment_levels) %>%
  mutate(n_prior_vaccinations = str_count(treatment, "V") - 1) %>%
  left_join(tibble(
    n_prior_vaccinations = -1:4,
    color = c("gray", brewer.pal(9, "Blues")[c(3, 5, 7, 8, 9)])
  )) %>%
  select(-n_prior_vaccinations)

treatment_colors <- setNames(treatment_colors$color, nm = treatment_colors$treatment)
year_4_nprior_vax_colors <- treatment_colors[c("PPPV", "PPVV", "PVVV", "VVVV")]
names(year_4_nprior_vax_colors) <- 0:3


# =====  Figure specifications =====
subtype_colors = c(H3N2 = "#6A3DA9", H1N1 = "#fe7b3f", B = "#A93D49")

default_figure_font_size = 8
small_font_size = 8 # We've made these the same.

baseline_figure_settings <- theme(
   legend.text = element_text(size = default_figure_font_size),
   legend.title = element_text(size = default_figure_font_size),
   axis.text = element_text(size = default_figure_font_size),
   axis.title = element_text(size = default_figure_font_size),
   strip.text = element_text(size = default_figure_font_size),
   axis.line = element_line(linewidth = 0.3)
)

main_text_point_size <- 0.2
main_text_alpha <- 0.2

boxplot_point_size <- 0.8
boxplot_point_alpha <- 0.4
boxplot_point_stroke <- 2
boxplot_point_shape <- 16
boxplot_line_width <- 0.2
boxplot_median_width <- 0.75

default_col_width <- 0.5

main_text_width_full <- 7.25

supp_fig_witdh_full <- 8.5 - 2 # 8.5 inches with 1 inch margins
supp_fig_height_full <- 11 - 2 # 11 inches with 1 inch margins