library(tidyverse)
library(readxl)
library(writexl)
source("R/base_functions.R")

# ================================== File paths ================================

data_folder_path <- "~/Box/DRIVE_data"

treatment_assignment_file <- file.path(data_folder_path, "participant_data/treatment_assignment.csv")

superparticipant_file_drive1 <- file.path(data_folder_path, "participant_data/DRIVE1_super_participants.csv")

drive_1_sera_file <- file.path(data_folder_path, "sample_usage/DRIVE1_sera.csv")
drive_2_sera_file <- file.path(data_folder_path, "sample_usage/DRIVE2_sera.csv")

FRNT_titer_file_main_format <- file.path(data_folder_path, "FRNT/H3N2_single_table_main_format.csv")
FRNT_titer_file_20250711_single_measurements <- file.path(data_folder_path, "FRNT/DRIVE_FRNT_results_20250711_single_measurements.csv")
FRNT_titer_file_20250711_double_measurements <- file.path(data_folder_path, "FRNT/DRIVE_FRNT_results_20250711_double_measurements.csv")

H1N1_FRNT_titers_file <- file.path(data_folder_path, "FRNT/H1N1_FRNT/DRIVE_H1N1_FRNT_20250801.csv")
H3N2_HAI_file <- file.path(data_folder_path, "HAI_titers/H3N2_HAI/H3N2_Y4_HAI.csv")

HAI_titers_batches_1to4_files <- list.files(path = file.path(data_folder_path, "HAI_titers/batches_1-4/"),
                                            pattern = 'titers_batch_[0-9]\\.csv', full.names = T)

HAI_titers_batches_5to7_file <- file.path(data_folder_path, "HAI_titers/batches_5-7/HAI_titers_batches_5-7.csv")

luminex_file <- file.path(data_folder_path, "Luminex/DRIVE_year3_luminex_results_summary.csv")

# Dates of covid vaccination (DRIVE-I only)
cov2_vdates_file <- file.path(data_folder_path, "participant_data/DRIVE_1_covid_vaxdates.csv")

HK_flu_surveillance_file <- file.path(data_folder_path, "Infections/HK_govt_flu_surveillance.csv")

NAI_infections_H1N1_path <- "results/NAI_inferred_infections/Inferred_N1_infections.csv"
NAI_infections_H3N2_path <- "results/NAI_inferred_infections/Inferred_N2_infections.csv"

timepoint_days_and_covariates_DRIVE1_file <- file.path(data_folder_path, "participant_data/DRIVE_1_timepoint_days_and_other_covariates.csv")
timepoint_days_and_covariates_DRIVE2_file <- file.path(data_folder_path, "participant_data/DRIVE_2_timepoint_days_and_other_covariates.csv")
vaccination_dates_DRIVE1_file <- file.path(data_folder_path, "participant_data/DRIVE_1_vaccination_dates.csv")
vaccination_dates_DRIVE2_file <- file.path(data_folder_path, "participant_data/DRIVE_2_vaccination_dates.csv")

age_and_sex_DRIVE1_file <- file.path(data_folder_path, "participant_data/age_and_sex_DRIVE1.csv")
age_and_sex_DRIVE2_file <- file.path(data_folder_path, "participant_data/age_and_sex_DRIVE2.csv")

flu_and_covid_PCR_DRIVE1_file <- file.path(data_folder_path, "Infections/flu_and_covid_PCR_DRIVE_1.csv")
flu_and_covid_PCR_DRIVE2_file <- file.path(data_folder_path, "Infections/flu_and_covid_PCR_DRIVE_2.csv")


# Scores from Tsang lab
tsang_scores_path <- file.path(data_folder_path, "Yale_scores/DRIVE_donors_module_scores.csv")

# ==================== Reading / initial processing functions ==================

# Function linking ids to treatment assignment
read_treatment_assignment <- function(treatment_assignment_file){
  treatment_assignment <- read_csv(treatment_assignment_file, show_col_types = FALSE) %>%
    dplyr::rename(pID = DRIVE_pid)
  
  treatment_codes <- treatment_assignment %>%
    select(intervention.group, drive, Y1, Y2, Y3, Y4, Y5) %>%
    arrange(intervention.group) %>%
    unique() %>%
    mutate(across(matches('Y'), function(x){str_sub(x, 1, 1)})) %>%
    mutate(treatment_Y1 = toupper(Y1),
           treatment_Y2 = toupper(paste0(Y1, Y2)),
           treatment_Y3 = toupper(paste0(Y1, Y2, Y3)),
           treatment_Y4 = toupper(paste0(Y1, Y2, Y3, Y4)),
           treatment_Y5 = toupper(paste0(Y1, Y2, Y3, Y4, Y5))) %>%
    select(intervention.group, drive, matches('treatment_')) %>%
    mutate(treatment_Y5 = case_when(
      str_detect(treatment_Y5, 'NA') ~ NA_character_,
      T ~ treatment_Y5
    ))
  
  return(
    left_join(treatment_assignment, treatment_codes)
  )
}

read_superparticipant_status <- function(superparticipant_file_drive1){
  superparticipant_status <- read_csv(superparticipant_file_drive1) %>%
    pivot_longer(cols = !matches('DRIVE_pid'), names_to = 'year_timepoint', values_to = 'provided_sample') %>%
    separate(year_timepoint, into = c('year', 'timepoint'), sep = '\\.') %>%
    mutate(year = as.integer(str_remove(year, 'y')),
           timepoint = as.integer(str_remove(timepoint, 'd'))) %>%
    rename(pID = DRIVE_pid)
}

# Reads precise timepoint days, vaccination dates
read_timepoint_days <- function(timepoint_days_and_covariates_file,
                                vaccination_dates_file){
  
  vaccination_dates <- read_csv(vaccination_dates_file) %>%
    rename(pID = DRIVE_pid) %>%
    filter(pID != "#N/A") %>%
    pivot_longer(cols = matches('vdate'),
                 names_transform = ~as.integer(str_remove(.x, "vdate.y")),
                 names_to = "year",
                 values_to = "year_vdate") %>%
    mutate(year_vdate = lubridate::dmy(year_vdate))
  
  tp_days <- read_csv(timepoint_days_and_covariates_file) %>%
    rename(pID = DRIVE_pid) %>%
    # Make it explicit that counting starts on d0 of year 1:
    mutate(y1d0 = 0) %>%
    select(pID, matches('y[0-9]+d')) %>%
    # ndays will designate the true number of days since day 0 of year 1
    # timepoint will designate the "official" associated time point 
    # (e.g. ndays for y1 day 30 may be 28, 31, 30, etc.)
    # (ndays for y2 day 0 may be 365, 357, 380, etc.)
    pivot_longer(cols = !any_of("pID"), values_to = 'ndays_since_y1d0') %>%
    separate(name, into = c('year', 'timepoint'), sep = 'd') %>%
    mutate(year = as.integer(str_remove(year, 'y')),
           timepoint = as.integer(timepoint))

  # Add a 365 timepoint equal to the d0 of the subsequent year
  d365_days <- tp_days %>%
    filter(year != 1, timepoint == 0) %>%
    mutate(timepoint = 365, year = year - 1)
  
  tp_days <- tp_days %>%
    bind_rows(d365_days) %>%
    arrange(pID, year, timepoint) %>%
    group_by(pID, year) %>%
    mutate(ndays_since_year_vax = ndays_since_y1d0 -
             ndays_since_y1d0[timepoint == 0]) %>%
    ungroup()
  
  return(tp_days %>%
           left_join(vaccination_dates %>%
                      group_by(pID) %>%
                      mutate(y1_vdate = year_vdate[year == 1]) %>%
                      ungroup()
                      ) %>%
           mutate(timepoint_date = y1_vdate + ndays_since_y1d0) %>%
           select(pID, year, year_vdate, everything()))
  
}

read_other_covariates <- function(timepoint_days_and_covariates_file){
  other_covariates <- read_csv(timepoint_days_and_covariates_file) %>%
    rename(pID = DRIVE_pid) %>%
    select(pID, asthma, smoking, height, weight)

  # Check that all values are either 1, 0 or 88(unknown)
  stopifnot(all(other_covariates$asthma %in% c(1, 0, 88)))
  stopifnot(all(other_covariates$smoking %in% c(1, 0, 88, 99)))

  other_covariates <- other_covariates %>%
    # Values coded as "unknown" (88) or "refused" (99) will be coded as 0 (absent)
    mutate(asthma = asthma == 1,
           smoking = smoking == 1,
           # Height in in cm, hence the /100
          BMI = weight / ((height / 100)^2))
          
  return(other_covariates)
}

read_HAI_titers_batches <- function(HAI_titers_batches_file){
  HAI_titers_batches <- read_csv(HAI_titers_batches_file) %>%
    # "No." column doesn't mean anything
    select(-any_of("No.")) %>%
    # Fill in DRIVE pids (which were entered only in the first row for each ID
    # in some of the raw data files)
    fill(all_of("DRIVE_pid"), .direction = "down") %>%
    rename(pID = DRIVE_pid) %>%
    # Year number accompanied by an "r" in some cases, remove.
    mutate(year = as.integer(str_remove(year,"r")),
           # Removes (73) from once instance where instance where timepoint was
           # coded "D273 (73)"
           timepoint = str_remove(timepoint, " \\(73\\)"), 
           timepoint = as.integer(str_remove(timepoint, "[DdC]"))) %>%
    pivot_longer(cols = matches("/"), names_to = 'strain', values_to = 'titer') %>%
    # Fixes a value incorrectly coded as 540
    mutate(titer = case_when(
      titer == 540 ~ 640,
      titer != 540 ~ titer
    )) %>%
    # In some files a ". " is used instead of / ("B. Austria", etc.)
    mutate(strain = str_replace(strain, "\\. ","/")) %>%
    mutate(log2_titer = log(titer, base = 2),
           titer_type = 'HAI') %>%
    compute_fold_changes() %>%
    mutate(batch = as.character(batch)) %>%
    select(pID, year, timepoint, titer_type, batch, strain, titer, log2_titer, fold_change)
  
  check_titer_values(HAI_titers_batches$titer)
  
  return(HAI_titers_batches)
    
}


# Some files (e.g., from Hensley lab) have a column combining year and time point
# This function splits that column into separate year and time point columns
split_time_point_id <- function(data){
  data %>%
    mutate(timepoint = str_extract(time_point_id, 'D[0-9]+'),
           timepoint = as.integer(str_remove(timepoint, 'D')),
           year = str_extract(time_point_id, 'y[0-9]+'),
           year = as.integer(str_remove(year, 'y')))
}

read_FRNT_titers_main_format <- function(FRNT_titer_file){
  FRNT_titer <- read_csv(FRNT_titer_file, show_col_types = F)
  
  FRNT_titer <- FRNT_titer %>%
    # Re-compute GMT (not all data transfers came with it pre-computed)
    select(-GMT) %>%
    mutate(GMT = case_when(
      is.na(FRNT_assay_2) ~ FRNT_assay_1,
      !is.na(FRNT_assay_2) ~ sqrt(FRNT_assay_1*FRNT_assay_2)
    )) %>%
    mutate(titer = GMT, titer_type = 'FRNT (GMT of 1+ measurements)') %>%
    split_time_point_id() %>%
    rename(pID = DRIVE_pid) %>%
    mutate(log2_titer = log(titer, base = 2)) %>%
    mutate(batch = as.character(batch)) %>%
    select(pID, year, timepoint, titer_type, batch, strain, titer, log2_titer, matches('_assay_'))
  
  check_titer_values(FRNT_titer$titer)
  
  return(FRNT_titer)
}

read_FRNT_titers_second_format <- function(FRNT_titer_file_other_format){
  read_csv(FRNT_titer_file_other_format, show_col_types = FALSE) %>%
    select(-matches("osition")) %>%
    pivot_longer(matches("FRNT"), values_to = "titer", names_to = "strain") %>%
    mutate(titer = case_when(
      titer == "x" ~ NA_character_,
      titer == "<20" ~ "10",
      titer == "≥1280" ~ "1280",
      T ~ titer
    )) %>%
    mutate(Year = as.character(Year)) %>%
    mutate(Year = case_when(
      Year == "z" ~ "2", # Typo. Verified based on time_point_sample_id that this is year 2
      T ~ Year,
    )) %>%
    mutate(Year = as.integer(Year),
           titer = as.integer(titer),
           log2_titer = log2(titer),
           FRNT_assay_1 = NA,
           FRNT_assay_2 = NA,
           titer_type = "FRNT (GMT of 1+ measurements)",
           strain = str_remove(strain, "_FRNT90"),
           strain = str_replace(strain, "HongKong", "Hong Kong")) %>%
    rename(year = Year, timepoint = Timepoint, batch = Batch, pID = DRIVE_pid) %>%
    mutate(timepoint = as.integer(timepoint), batch = as.character(batch)) %>%
    select(pID, year, timepoint, titer_type, batch, strain, titer, log2_titer, matches('assay')) %>%
    # Many rows with missing titers, corresponding to titers measured previously and shown in the 
    # main format table. Remove from here
    filter(!is.na(titer))
}

read_FRNT_titers_third_format <- function(FRNT_titer_file_third_format){
  read_csv(FRNT_titer_file_third_format, show_col_types = FALSE) %>%
    select(-matches("osition")) %>%
    pivot_longer(matches("FRNT"), values_to = "titer", names_to = "strain") %>%
    mutate(titer = case_when(
      titer == "<20" ~ "10",
      titer == "≥1280" ~ "1280",
      T ~ titer
    )) %>%
    mutate(titer = as.integer(titer), 
           replicate = str_extract(strain, "rep[0-9]+_"),
           strain = str_remove(strain, replicate),
           replicate = str_extract(replicate, "[0-9]+")) %>%
    pivot_wider(names_from = "replicate", values_from = "titer", names_prefix = "FRNT_assay_") %>%
    mutate(titer = sqrt(FRNT_assay_1 * FRNT_assay_2),
           log2_titer = log2(titer),
            titer_type = "FRNT (GMT of 1+ measurements)",
           strain = str_remove(strain, "_FRNT90"),
           strain = str_replace(strain, "HongKong", "Hong Kong")) %>%
    rename(year = Year, timepoint = Timepoint, batch = Batch, pID = DRIVE_pid) %>%
    mutate(timepoint = as.integer(timepoint), batch = as.character(batch)) %>%
    select(pID, year, timepoint, titer_type, batch, strain, titer, log2_titer, matches('assay')) %>%
    # Many rows with missing titers, corresponding to titers measured previously and shown in the 
    # main format table. Remove from here
    filter(!is.na(titer))      
}

# FRNT titers to H1N1 also have a slightly different format
read_FRNT_titers_to_H1N1 <- function(H1N1_FRNT_titers_file){
  read_csv(H1N1_FRNT_titers_file, show_col_types =  F) %>%
    select(-matches("osition")) %>%
    pivot_longer(matches("FRNT"), values_to = "titer", names_to = "strain") %>%
    mutate(titer = case_when(
      titer == "<20" ~ "10",
      titer == "≥5120" ~ "5120",
      T ~ titer
    )) %>%
    mutate(titer = as.integer(titer),
           log2_titer = log2(titer),
           strain = str_remove(strain, "_FRNT"),
           titer_type = "FRNT (GMT of 1+ measurements)",
           FRNT_assay_1 = titer, 
           FRNT_assay_2 = NA) %>%
    rename(year = Year, timepoint = Timepoint, batch = Batch, pID = DRIVE_pid) %>%
    mutate(timepoint = as.integer(timepoint), batch = as.character(batch)) %>%
    select(pID, year, timepoint, titer_type, batch, strain, titer, log2_titer, matches('assay'))
}

# Another special format for HAI titers to H3N2
read_HAI_titers_to_H3N2 <- function(H3N2_HAI_file, serum_samples){
  H3N2_HAI_titers <- read_csv(H3N2_HAI_file) %>%
    # Compute the GMT across replicated measurements
    mutate(titer = sqrt(SG_Titer * GL_Titer),
           log2_titer = log2(titer),
           titer_type = "HAI") %>%
    select(-SG_Titer, -GL_Titer) %>%
    rename(sid = sample_id) %>%
    left_join(serum_samples %>%
               select(sid, pID, year, timepoint), by = join_by(sid)) %>%
    compute_fold_changes() %>%
    mutate(fourfold_rise_or_greater = fold_change >= 4,
           titer_40_or_greater = titer >= 40) %>%
    mutate(strain = "A/Darwin/9/2021") %>%
    select(-sid)


 # Stop if any year != 4   
 stopifnot(all(H3N2_HAI_titers$year == 4))
 
 return(H3N2_HAI_titers)
}

read_luminex <- function(luminex_file){
  read_csv(luminex_file) %>%
    rename(pID = DRIVE_pid) %>%
    split_time_point_id() %>%
    select(-matches('intervention'), -matches("sample_id"),
           -prescr_id, -time_point_id) %>%
    select(pID, year, timepoint, everything()) %>%
    pivot_longer(cols = !any_of(c('pID', 'year', 'timepoint')),
                 names_to = 'strain_and_replicate', 
                 values_to = "luminex") %>%
    separate("strain_and_replicate", into = c("strain", "replicate"),
             sep = '_') %>%
    pivot_wider(names_from = 'replicate', 
                values_from = 'luminex', 
                names_prefix = 'luminex_') %>%
    arrange(pID, year, strain, timepoint) %>%
    rename(luminex = luminex_ave) %>%
    mutate(log10_luminex = log(luminex, base = 10))
  
  
}

read_positive_pcr_tests <- function(flu_and_covid_PCR_file, date_format_function){
  positive_tests <- read_csv(flu_and_covid_PCR_file) %>%
    rename(pID = DRIVE_pid) %>%
    # Some participants provided more than one sample on the same date for the 
    # same reported illness (1 saliva, 1 pooled nasal/throat swab)
    # Summarise those instances in a single result
    group_by(pID, cdate) %>%
    # Codes based on the data dictionary
    summarise(
      flu_A_positive = any(pcr_flu_pos == 1),
      flu_B_positive = any(pcr_flu_pos == 2),
      cov2_positive = any(pcr_cov2_pos == 1),
      H1N1_positive = any(pcr_flua_pos == 1),
      H3N2_positive = any(pcr_flua_pos == 2),
      `B/Victoria_positive` = any(pcr_flub_pos == 1),
      `B/Yamagata_positive` = any(pcr_flub_pos == 2)
      ) %>%
    ungroup() %>%
    pivot_longer(cols = matches('positive'),
                 names_to = 'result') %>%
    filter(value == T) %>%
    select(-value) %>%
    mutate(cdate = date_format_function(cdate))
  
  #NOTE: this object has multiple rows for the same test.
  # For isntance, a test positive for H3N2 will have one row
  # with result "flu_A_positive" and another row with result
  # "H3N2_positive"  
  return(positive_tests)
}

# This function takes a data set and adds a d365 whenever possible, using data
# from day 0 in the subsequent year. read_timepoint_days does a similar operation
# to determine the precise dates/n days since vaccination for the d365 time point.
get_d365_from_next_year_d0 <- function(data){
  
  d365_data <- data %>%
    # Take time point 0 for all years except the first year
    filter(timepoint == 0, year != 1) %>%
    # Re-set the timepoint as 365 and subtract 1 from the year
    mutate(timepoint = 365, year = year - 1)
  
  # Bind rows with original dataset
  return(data %>%
           bind_rows(d365_data) %>%
           arrange(year, timepoint))

}

read_serum_samples <- function(serum_samples_file, date_format_function){
  read_csv(serum_samples_file) %>%
    rename(pID = DRIVE_pid) %>%
      # Process Timepoint column
      mutate(year = str_extract(Timepoint, "Year [0-9]+") %>%
              str_remove("Year ") %>%
              as.integer(),
             timepoint = str_extract(Timepoint, "Day [0-9]+") %>%
              str_remove("Day ") %>%
              as.integer()) %>%
      select(-Timepoint) %>%
      mutate(cdate = date_format_function(cdate))
}

read_cov2_vdates <- function(covid_vdates_file){
  read_csv(covid_vdates_file) %>%
    mutate(across(matches("date"), dmy)) %>%
    rename(pID = DRIVE_pid) %>%
    rename_with(~str_replace(., "^d([0-9]+)_date$", "cov2_vdate_d\\1"), matches("^d[0-9]+_date$"))
}

read_HK_flu_surveillance <- function(HK_flu_surveillance_file){
  read_csv(HK_flu_surveillance_file) %>%
    mutate(From = dmy(From), To = dmy(To)) %>%
    mutate(
      H1_pos_1000visits = ILIrate * H1_proportion,
      H3_pos_1000visits = ILIrate * H3_proportion,
      B_pos_1000visits  = ILIrate * B_proportion
  )
}

# ================================= Load data ==================================

# When this script is first called by the script that infers infections NAI titers,
# these files won't exist yet
NAI_inferred_infections <- NULL
if(file.exists(NAI_infections_H1N1_path) & file.exists(NAI_infections_H3N2_path)){
  NAI_inferred_infections <- bind_rows(
    read_csv(NAI_infections_H1N1_path) %>% mutate(result = "H1N1_positive"),
  read_csv(NAI_infections_H3N2_path) %>% mutate(result = "H3N2_positive")
  ) %>%
  # So that we can use the same functions that work on positive_pcr_tets, we
  # rename infection_date as cdate
  rename(cdate = infection_date)
}

superparticipant_status_drive1 <- read_superparticipant_status(superparticipant_file_drive1)
treatment_assignment <- read_treatment_assignment(treatment_assignment_file)
cov2_vdates <- read_cov2_vdates(cov2_vdates_file)

timepoint_days <-
  bind_rows(
    read_timepoint_days(timepoint_days_and_covariates_DRIVE1_file, vaccination_dates_DRIVE1_file),
    read_timepoint_days(timepoint_days_and_covariates_DRIVE2_file, vaccination_dates_DRIVE2_file)
  )


age_and_sex_DRIVE1 <- read_csv(age_and_sex_DRIVE1_file) %>%
  rename(pID = DRIVE_pid,
         age_in_year_1 = scr_age,
         sex = scr_sex) %>%
  mutate(sex = recode_values(sex,
                          1 ~ 'M',
                          2 ~ 'F'))

age_and_sex_DRIVE2 <- read_csv(age_and_sex_DRIVE2_file) %>%
  rename(pID = DRIVE_pid,
         age_in_year_1 = scr_age,
         sex = scr_male) %>%
  mutate(sex = recode_values(sex,
                          1 ~ 'M',
                          0 ~ 'F'))

age_and_sex <- bind_rows(age_and_sex_DRIVE1, age_and_sex_DRIVE2)                                          

rm(age_and_sex_DRIVE1)
rm(age_and_sex_DRIVE2)

# Read additional covariates besides age/sex, infection-derived ones
other_covariates <-
  bind_rows(
    read_other_covariates(timepoint_days_and_covariates_DRIVE1_file),
    read_other_covariates(timepoint_days_and_covariates_DRIVE2_file)
  )

positive_pcr_tests <-
  bind_rows(
    read_positive_pcr_tests(flu_and_covid_PCR_DRIVE1_file, date_format_function = dmy),
    read_positive_pcr_tests(flu_and_covid_PCR_DRIVE2_file, date_format_function = dmy)
  )

# Table of serum samples collected 
serum_samples <- bind_rows(
  read_serum_samples(drive_1_sera_file, date_format_function = mdy),
  read_serum_samples(drive_2_sera_file, date_format_function = mdy)
  ) %>%
  annotate_with_treatment(treatment_assignment) %>%
  annotate_with_timepoint_days(timepoint_days) %>%
  annotate_with_previous_infections(positive_pcr_tests, NAI_inferred_infections) %>%
  annotate_with_age_and_sex(age_and_sex)

FRNT_titers_main_format <- read_FRNT_titers_main_format(FRNT_titer_file_main_format)
FRNT_titers_20250711_single_measurements <- read_FRNT_titers_second_format(FRNT_titer_file_20250711_single_measurements) %>%
  mutate(FRNT_assay_1 = titer)
FRNT_titers_20250711_double_measurements <- read_FRNT_titers_third_format(FRNT_titer_file_20250711_double_measurements)

# Measurements in FRNT_titers_20250711_single_measurements are shown again in FRNT_titers_20250711_double_measurements
# Filter those out to avoid duplication
FRNT_titers_20250711_single_measurements <- FRNT_titers_20250711_single_measurements %>%
  anti_join(FRNT_titers_20250711_double_measurements, by = c("pID", "strain", "year", "timepoint", "batch"))

FRNT_titers <- bind_rows(
  FRNT_titers_main_format,
  FRNT_titers_20250711_single_measurements,
  FRNT_titers_20250711_double_measurements
)

# Check each pID/strain/year/timepoint/batch present in FRNT_titers only once
stopifnot(
  FRNT_titers %>%
    group_by(pID, strain, year, timepoint, batch) %>%
    summarise(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>%
    nrow() == 0
)

FRNT_titers <- FRNT_titers %>%
  compute_fold_changes()

HAI_titers_batches_1to4 <- bind_rows(lapply(as.list(HAI_titers_batches_1to4_files),
                                  read_HAI_titers_batches))

HAI_titers_batches_5to7 <- read_HAI_titers_batches(HAI_titers_batches_5to7_file)


titers <- bind_rows(list(FRNT_titers,
                         HAI_titers_batches_1to4,
                         HAI_titers_batches_5to7)) %>%
  mutate(fourfold_rise_or_greater = fold_change >= 4,
  titer_40_or_greater = titer >= 40) %>%
  get_d365_from_next_year_d0() %>%
  annotate_with_treatment(treatment_assignment) %>%
  annotate_with_timepoint_days(timepoint_days) %>%
  annotate_with_age_and_sex(age_and_sex) %>%
  left_join(other_covariates)

# Annotate and order strains and subtypes  
strain_levels <- unique(titers$strain)

titers <- titers %>%
  annotate_and_order_strains_and_subtypes(strain_levels) %>%
  annotate_with_previous_infections(positive_pcr_tests, NAI_inferred_infections) %>%
  annotate_with_cov2_vaccination(cov2_vdates = cov2_vdates)

# FRNT titers to H1N1 are handled separately, not added to the main titers object
H1N1_FRNT_titers <- read_FRNT_titers_to_H1N1(H1N1_FRNT_titers_file) %>%
  annotate_with_treatment(treatment_assignment) %>%
  annotate_and_order_strains_and_subtypes(strain_levels = NULL) %>%
  annotate_with_timepoint_days(timepoint_days) %>%
  annotate_with_previous_infections(positive_pcr_tests, NAI_inferred_infections)

H3N2_HAI_titers <- read_HAI_titers_to_H3N2(H3N2_HAI_file, serum_samples) %>%
  annotate_with_treatment(treatment_assignment) %>%
  annotate_and_order_strains_and_subtypes(strain_levels = NULL) %>%
  annotate_with_timepoint_days(timepoint_days) %>%
  annotate_with_previous_infections(positive_pcr_tests, NAI_inferred_infections)

luminex <- read_luminex(luminex_file) %>%
  get_d365_from_next_year_d0() %>%
  annotate_with_timepoint_days(timepoint_days) %>%
  annotate_with_treatment(treatment_assignment) %>%
  annotate_with_age_and_sex(age_and_sex)

luminex_vax_strains <- luminex %>%
  #This keeps only vaccine strains, standardizes some strain names to match titer data
  filter_luminex_to_vax_strains() %>%
  annotate_and_order_strains_and_subtypes(strain_levels = NULL) %>%
  annotate_with_previous_infections(positive_pcr_tests, NAI_inferred_infections)

HK_flu_surveillance <- read_HK_flu_surveillance(HK_flu_surveillance_file)