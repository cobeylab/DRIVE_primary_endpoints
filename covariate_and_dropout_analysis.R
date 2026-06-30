library(rlang)
library(survival)
library(coxphf)
source("R/load_data.R")

# Code by A.S.

# =====================================================
#                 Function definitions
# =====================================================

read_covariate_data <- function(covariate_data_file){
  read_csv(covariate_data_file) %>% 
    dplyr::select(!matches("^y.*\\d")) %>%
    rename(pID = DRIVE_pid) %>% 
    #the following operation makes ethnicity, edu, marital_status, hh.member, hypertension, asthma, depression, digestive, hypothryoidism, derm, smoking, current.smoking columns to be changed to a "character" column; 
    #the other columns are numerical
    mutate(across(-c(weight, height), ~ifelse(. == 88, "Unknown", .))) %>% 
    mutate(across(-c(weight, height), ~ifelse(. == 99, "Refused", .))) %>%
    #these are for the binary variables
    mutate(across(c(hypertension, asthma, depression, digestive, hypothryoidism, derm, smoking, current.smoking),  ~case_when(
      . == "1" ~ "Yes", 
      . == "0" ~ "No"))) %>% 
    mutate(ethnicity = case_when(
      ethnicity == 1 ~ "Chinese",
      ethnicity == 2 ~ "Indonesian",
      ethnicity == 3 ~ "Filipino",
      ethnicity == 4 ~ "White",
      ethnicity == 5 ~ "Indian",
      ethnicity == 6 ~ "Pakistani",
      ethnicity == 7 ~ "Nepalese",
      ethnicity == 8 ~ "Japanese",
      ethnicity == 9 ~ "Thai",
      ethnicity == 10 ~ "Other Asian",
      ethnicity == 66 ~ "Other")) %>% 
    mutate(edu = case_when(
      edu == 1 ~ "No school training",
      edu == 2 ~ "Primary school",
      edu == 3 ~ "Form 1-3",
      edu == 4 ~ "Form 4-5",
      edu == 5 ~ "Advanced level",
      edu == 6 ~ "College/Undergraduate",
      edu == 7 ~ "Postgraduate or above")) %>% 
    mutate(marital = case_when(
      marital == 1 ~ "Married",
      marital == 2 ~ "Separated/divorced",
      marital == 3 ~ "Widowed",
      marital == 4 ~ "Never married")) %>% 
    mutate(across(c(hh.member1.agegp, hh.member2.agegp, hh.member3.agegp, hh.member4.agegp, hh.member5.agegp, hh.member6.agegp, hh.member7.agegp), 
                  ~ case_when(
                    . == 1 ~ "0-4 years",
                    . == 2 ~ "5-9 years",
                    . == 3 ~ "10-19 years",
                    . == 4 ~ "20-29 years",
                    . == 5 ~ "30-39 years",
                    . == 6 ~ "40-49 years",
                    . == 7 ~ "50-59 years",
                    . == 8 ~ "60-69 years",
                    . == 9 ~ "70-74 years",
                    . == 10 ~ "75-79 years",
                    . == 11 ~ "80-84 years",
                    . == 12 ~ "85+ years",
                    . == 13 ~ "85+ years",
                    . == 14 ~ "85+ years"
                  ))) %>% 
    mutate(BMI = weight/(0.01*height)^2)
}

read_DRIVE_income_data <- function(DRIVE_income_file){
  read_excel(DRIVE_income_file, sheet = "Data") %>% 
    rename(pID = DRIVE_pid) %>% 
    mutate(income = case_when(
      income == 1 ~ "<$4,000",
      income == 2 ~ "$4,000-5,999", 
      income == 3 ~ "$6,000-7,999",
      income == 4 ~ "$8,000-9,999", 
      income == 5 ~ "$10,000-14,999", 
      income == 6 ~ "$15,000-19,999", 
      income == 7 ~ "$20,000-24,999", 
      income == 8 ~ "$25,000-29,999", 
      income == 9 ~ "$30,000-34,999", 
      income == 10 ~ "$35,000-39,999", 
      income == 11 ~ "$40,000-44,999", 
      income == 12 ~ "$45,000-49,999", 
      income == 13 ~ "$50,000-59,999", 
      income == 14 ~ "$60,000-79,999", 
      income == 15 ~ "$80,000-99,999",
      income == 16 ~ ">$100,000",
      income == 88 ~ "Unknown"
    ))
}

read_DRIVE_med_data <- function(DRIVE_med_file){
  read_excel(DRIVE_med_file, sheet = "Sheet1") %>% 
    rename(pID = DRIVE_pid)
}

read_DRIVE_withdrawal_data <- function(DRIVE_withdrawal_file){
  read_excel(DRIVE_withdrawal_file, sheet = "Sheet1") %>% 
    rename(pID = DRIVE_pid) %>% 
    mutate(wdf_type = case_when(
      wdf_type == 1 ~ "Participant requested withdrawal",
      wdf_type == 2 ~ "Family member requested withdrawal on participant's behalf",
      wdf_type == 3 ~ "Do not meet the inclusion criteria",
      wdf_type == 66 ~ "other")) %>%
    mutate(wdf_type_spec = case_when(
      str_detect(wdf_type_spec, regex("^lost", ignore_case = T)) ~ "lost to follow up",
      T ~ wdf_type_spec
    )) %>% 
    mutate(withdrawal_year = year(wdf_wddate))
}

##generalize covariate comparison for enrolled and dropped-out participants
#select covariates that are categorical: age-group, sex, income, edu, occupation, asthma, hypothryoidism, derm, ethnicity, smoking, marital 
#select covariates that are continuous: age, BMI, hh.member
#do chi-sq test on categorical variables and KS test on continuous variables
#output as a list of covariate dataframes and a list of p.values with statistical tests
cov_comp_discrete <- function(covariate, cov_data_start, cov_data_drop_out){
  covariate_var <- enquo(covariate)
  covariate_name <- as_name(covariate_var)
  
  cov_dist_start <- cov_data_start %>%
    filter(!is.na(!!covariate_var) & !!covariate_var != "Unknown") %>% 
    group_by(!!covariate_var) %>% 
    summarise(num_indivs_start = n()) %>% 
    mutate(percent_start = round(100*num_indivs_start/sum(num_indivs_start), 1)) %>%
    mutate(num_and_per_start = paste(num_indivs_start, "(", percent_start, "%)", sep = " ")) %>% 
    ungroup()
  
  cov_dist_drop <- cov_data_drop_out %>%
    filter(!is.na(!!covariate_var) & !!covariate_var != "Unknown") %>% 
    group_by(!!covariate_var) %>% 
    summarise(num_indivs_drop = n()) %>% 
    mutate(percent_drop = round(100*num_indivs_drop/sum(num_indivs_drop), 1)) %>%
    mutate(num_and_per_drop = paste(num_indivs_drop, "(", percent_drop, "%)", sep = " ")) %>% 
    ungroup()
  
  cov_dist <- left_join(cov_dist_start, cov_dist_drop) 
  cov_dist$num_indivs_drop[is.na(cov_dist$num_indivs_drop)] <- 0
  
  cov_dist_long <- cov_dist %>% 
    pivot_longer(cols = c(num_indivs_start, num_indivs_drop),
                 names_to = "timepoint",
                 values_to = "number")
  
  formula_str <- paste0("number", " ~ ", covariate_name, " + ", "timepoint")
  formula_obj <- as.formula(formula_str)
  
  cont_table_cov <- xtabs(formula_obj, data = cov_dist_long)
  
  if (min(dim(cont_table_cov)) < 2) {
    chisq_test_results <- NA
    chisq_pvalue <- NA
    fisher_test_results <- NA
    fisher_test_pvalue <- NA
  } else {
    chisq_test_results <-  tryCatch(chisq.test(cont_table_cov), error = function(e) NA)
    chisq_pvalue <- if (is.list(chisq_test_results)) round(chisq_test_results$p.value, 5) else NA
    
    fisher_test_results <- tryCatch(fisher.test(cont_table_cov, simulate.p.value=T), error = function(e) NA)
    fisher_test_pvalue <- if (is.list(fisher_test_results)) round(fisher_test_results$p.value, 5) else NA
  }
  
  cov_dist <- cov_dist %>% 
    dplyr::select(-c(num_indivs_start, percent_start, num_indivs_drop, percent_drop)) %>% 
    mutate(chisq_pval = chisq_pvalue, fisher_pval = fisher_test_pvalue, covariate = covariate_name) 
  colnames(cov_dist)[1] <- "atribute"
  
  return(list(cov_dist, chisq_test_results, fisher_test_results))
  
}

#compare different covariates for the Y1 enrolled population and the remaining people after drop-outs for each year
#to check if the randomization across treatment groups holds TRUE for a covariate within a year
cov_comp_treat <- function(covariate, treatment_yr, covariate_data){ #this is for a given year
  covariate_var <- enquo(covariate)
  covariate_name <- as_name(covariate_var)
  treatment_var <- enquo(treatment_yr)
  treatment_name <- as_name(treatment_var)
  
  unique_covariate_values <- covariate_data %>% 
    filter(!is.na(!!covariate_var) & !!covariate_var != "Unknown" & !!covariate_var != "Refused") %>% 
    dplyr::select(!!covariate_var) %>% 
    unique()
  
  #if (length(unique_covariate_values) > 1){
  cov_data_table <- covariate_data %>% 
    filter(!is.na(!!covariate_var) & !!covariate_var != "Unknown" & !!covariate_var != "Refused") %>% 
    group_by(!!treatment_var, !!covariate_var) %>% 
    summarise(num_indivs = n()) %>% 
    ungroup()
  
  cov_data_table_wide <- cov_data_table %>% 
    pivot_wider(names_from = treatment_name, values_from = num_indivs) 
  cov_data_total_by_treatment <- cov_data_table_wide %>% 
    summarise(across(!covariate_name, \(x) sum(x, na.rm = TRUE)))
  
  # prop.test_results <- pairwise.prop.test(x = unlist(cov_data_table_wide[1,-1]),
  #                                         n = unlist(cov_data_total_by_treatment),
  #                                         p.adjust.method = "none")
  
  formula_str <- paste0("num_indivs", " ~ ", covariate_name, " + ", treatment_name)
  formula_obj <- as.formula(formula_str)
  
  cont_table_cov <- xtabs(formula_obj, data = cov_data_table)
  
  if (min(dim(cont_table_cov)) < 2) {
    chisq_test_results <- NA
    chisq_pvalue <- NA
    fisher_test_results <- NA
    fisher_test_pvalue <- NA
  } else {
    chisq_test_results <- tryCatch(chisq.test(cont_table_cov), error = function(e) NA)
    chisq_pvalue <- if (is.list(chisq_test_results)) round(chisq_test_results$p.value, 2) else NA
    
    fisher_test_results <- tryCatch(fisher.test(cont_table_cov, simulate.p.value=TRUE), error = function(e) NA)
    fisher_test_pvalue <- if (is.list(fisher_test_results)) round(fisher_test_results$p.value, 2) else NA
  }
  
  
  return(list(cov_data_table_wide %>% 
                mutate(across(where(is.numeric) & !covariate_name, ~round(.x/sum(.x, na.rm = T)*100,2),.names = "{.col}_pct")) %>% 
                mutate(chisq_pval = chisq_pvalue, fisher_pval = fisher_test_pvalue) %>% 
                mutate(covariate = covariate_name) %>% 
                rename(attribute := !!sym(covariate_name)), 
              chisq_test_results, fisher_test_results)) #,prop.test_results
  # } else {
  #   return(list())
  # }
}

# ================================================================================
# Load covariate, income, medicine usage, withdrawal reason, PCR test results data
# ================================================================================

covariate_data_drive1_file <- file.path(data_folder_path, "/participant_data/DRIVE_1_timepoint_days_and_other_covariates.csv")
income_data_drive_file <- file.path(data_folder_path, "/participant_data/DRIVE_income.xlsx")
med_data_drive_file <- file.path(data_folder_path, "participant_data/DRIVE_postv_med.xlsx")
withdrawal_file_DRIVEI <- file.path(data_folder_path, "participant_data/DRIVE_withdrawal.xlsx")
pcr_test_file <- file.path(data_folder_path, "Infections/flu_and_covid_PCR_DRIVE_1.csv")

covariate_data_drive1 <- read_covariate_data(covariate_data_drive1_file)
income_data_drive <- read_DRIVE_income_data(income_data_drive_file)
med_data_drive <- read_DRIVE_med_data(med_data_drive_file)
#Y4 and Y5 drop-outs are in the withdrawal file data set
withdrawal_reasons_DRIVEI <- read_DRIVE_withdrawal_data(withdrawal_file_DRIVEI) %>% filter(str_detect(pID, "D1"))

pcr_test_results <- read_csv(pcr_test_file) %>%
  rename(pID = DRIVE_pid)  %>%
  mutate(cdate = dmy(cdate)) %>% 
  dplyr::select(c(pID, cdate, pcr_flu_pos, pcr_cov2_pos)) %>% 
  filter(pcr_flu_pos != 77)

pcr_test_results_summary <- pcr_test_results %>% 
  group_by(pID) %>% 
  summarise(infection_status_flu = ifelse(any(pcr_flu_pos %in% c(1, 2), na.rm = T),1,0),
            infection_status_cov2 = as.integer(max(pcr_cov2_pos, na.rm = T) == 1),
            .groups = "drop")

#cov2_vdates loaded from load_data.R
cov2_vdates$num_cov2vax <- rowSums(!is.na(cov2_vdates[, -1]))
min(cov2_vdates$num_cov2vax)
cov2_vdates$cov2_vax_status <- 1

# ======================================================
# Add relevant information to the loaded covariate data
# ======================================================

# treatment_assignment coming from load_data.R
covariate_data_drive1 <- covariate_data_drive1 %>%
  left_join(treatment_assignment %>% filter(drive == 1)) %>%
  left_join(age_and_sex) %>%
  left_join(income_data_drive %>% filter(pID %in% covariate_data_drive1$pID)) %>%
  left_join(med_data_drive %>% filter(pID %in% med_data_drive$pID)) %>%
  left_join(withdrawal_reasons_DRIVEI %>% dplyr::select(pID, wdf_wddate)) %>%  
  left_join(pcr_test_results_summary) %>%
  left_join(cov2_vdates %>% dplyr::select(pID, cov2_vax_status)) %>%
  mutate(cov2_vax_status = coalesce(cov2_vax_status, 0),
         infection_status_flu = coalesce(infection_status_flu, 0),
         infection_status_cov2 = coalesce(infection_status_cov2, 0)) %>%
  mutate(age_group_y1 = case_when(
    age_in_year_1 >= 18 & age_in_year_1 <= 25 ~ '18-25',
    age_in_year_1 > 25 & age_in_year_1 <= 35 ~ '26-35',
    age_in_year_1 > 35 & age_in_year_1 <= 45 ~ '36-45'
  ), BMI_group_y1 = case_when(
    round(BMI, 1) >= 14 & round(BMI,1) < 18.5 ~ '14-18.5',
    round(BMI, 1) >= 18.5 & round(BMI,1) < 25 ~ '18.5-25',
    round(BMI, 1) >= 25 & round(BMI,1) <= 30 ~ '25-30'
  )) %>% 
  #superparticipant_status_drive1 comes from load_data.R
  mutate(superparticipant_status = ifelse(pID %in% superparticipant_status_drive1$pID, 1, 0))


#vaccination and infection status of enrolled and drop-outs
#p ositive_pcr_tests generated by load_data.R
#timepoint_days generated by load_data.R
timepoint_days_DRIVEI <- timepoint_days %>% filter(pID %in% unique(covariate_data_drive1$pID))

cov2_vax_by_timepoint_summary <- timepoint_days_DRIVEI %>% 
  dplyr::select(pID, year, year_vdate, matches("timepoint")) %>%
  #unique() %>%
  left_join(cov2_vdates) %>%
  group_by(pID, year, timepoint, timepoint_date, year_vdate) %>% 
  summarise(
    cov2_vax_before_fluvax = any(c_across(starts_with("cov2_vdate")) < year_vdate, na.rm = T),
    cov2_vax_before_sample = any(c_across(starts_with("cov2_vdate")) < timepoint_date, na.rm = T),
    .groups = "drop")

#annotating timepoint_days_DRIVEI with cov2_vax status to check for infection and vaccination status by timepoint_date
timepoint_days_DRIVEI <- left_join(timepoint_days_DRIVEI, cov2_vax_by_timepoint_summary) %>% 
  annotate_with_previous_infections(positive_pcr_tests, NAI_inferred_infections = NAI_inferred_infections) %>% 
  dplyr::select(pID, year, year_vdate, matches("timepoint"), 
                cov2_vax_before_sample, infection_before_sample_anyflu_PCR, infection_before_sample_cov2_PCR)

# ======================================================
# Steps to find out drop outs by Y4 D182
# ======================================================

Y4D182_date <- timepoint_days_DRIVEI %>%
  filter(year == 4, timepoint == 182) %>%
  filter(!is.na(timepoint_date)) %>%
  pull(timepoint_date) %>%
  max()

dropout_status <- timepoint_days_DRIVEI %>%
  filter(year == 4, timepoint == 182) %>% 
  left_join(withdrawal_reasons_DRIVEI %>%
              dplyr::select(pID, wdf_wddate)) %>%
  mutate(dropout = case_when(
    !is.na(timepoint_date) ~ F, #excluding people who did not have a visit on Y4 D182
    is.na(wdf_wddate) ~ F,      #excluding people who never filled out the withdrawal form
    wdf_wddate >= Y4D182_date ~ F, #excluding people who filled out the withdrawal form later than the Y4 D182 date
    wdf_wddate < Y4D182_date ~ T,
  )) %>%
  dplyr::select(pID, dropout) %>%
  unique()

#table summarizing drop-out status by Y4 D182
dropout_status %>%
  group_by(dropout) %>%
  count()

covariate_data_drive1 <- covariate_data_drive1 %>% left_join(dropout_status)

covariates_for_drop_outs <- covariate_data_drive1 %>% filter(dropout == T)

# ======================================
# Investigating reasons for dropping out
# ======================================

wdf_type_df <- as.data.frame(table(withdrawal_reasons_DRIVEI$wdf_type[withdrawal_reasons_DRIVEI$pID %in% covariates_for_drop_outs$pID])) %>%
  mutate(
    prop = Freq / sum(Freq),
    perc = round(100 * prop, 1),
    Var1_wrapped = str_wrap(Var1, width = 22),
    label = paste0(Var1_wrapped, "\n", perc, "%"),
    ypos = sum(Freq) - 0.8*cumsum(Freq)   # midpoints for labeling
  )

ggplot(wdf_type_df, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(y = ypos, label = label),
            color = "white", size = 4, lineheight = 0.9) +
  scale_fill_manual(values = c("tomato3", "steelblue3", "gold3")) + 
  theme_void() +
  theme(legend.position = "none")

# ===============================================================
#     Steps to find out year-wise drop out participants IDs
# ===============================================================

# find the last year of attendance for the drop-outs
timepoint_days_DRIVEI_drop_outs <- timepoint_days_DRIVEI %>% 
  filter(pID %in% covariates_for_drop_outs$pID) %>% 
  filter(!is.na(year_vdate), !is.na(timepoint_date)) %>% 
  unique()
last_tp_drop_outs_DRIVEI <- timepoint_days_DRIVEI_drop_outs %>% #year in this case is year of drop-out for each participant
  group_by(pID) %>% 
  slice_tail(n = 1) %>%
  ungroup()
covariates_for_drop_outs <- covariates_for_drop_outs %>% left_join(last_tp_drop_outs_DRIVEI)
#creating a list of drop out IDs where the elements in the list is DRIVE-I year specific drop out IDs
drop_out_IDs_year <- split(covariates_for_drop_outs$pID, covariates_for_drop_outs$year)
#combining drop outs in later years with prior years to facilitate easy subsetting year wise
cumulative_drop_out_IDs <- accumulate(drop_out_IDs_year, c)
#filtering out year wise drop outs from the covariate data set and store subsetted data into a list
covariates_after_drop_outs <- map(cumulative_drop_out_IDs, ~filter(covariate_data_drive1, !pID %in% .x))
#filtering out year wise drop outs from the timepoint data and store subsetted data into a list
timepoint_days_DRIVEI_after_drop_outs <- map(cumulative_drop_out_IDs, ~filter(timepoint_days_DRIVEI %>% filter(timepoint == 30), !pID %in% .x))
# subsetting corresponding years 
timepoint_days_DRIVEI_after_drop_outs <- map2(timepoint_days_DRIVEI_after_drop_outs, seq_along(timepoint_days_DRIVEI_after_drop_outs),
                                              ~ .x[.x$year == .y + 1, , drop = FALSE])
covariates_after_drop_outs <- map2(covariates_after_drop_outs, timepoint_days_DRIVEI_after_drop_outs, ~ left_join(.x, .y))
y1_covariates <- left_join(covariate_data_drive1, timepoint_days_DRIVEI %>% filter(year == 1 & timepoint == 30))
#append the Y1 enrolled individuals covariate dataframe to the list of data frames created earlier
covariates_after_drop_outs <- append(list(y1_covariates), covariates_after_drop_outs)
names(covariates_after_drop_outs) <- c("Y1_enrolled", "Y2_enrolled", "Y3_enrolled", "Y4_enrolled", "Y5_enrolled")#, "final")

# ===============================================================
# Covariate comparison between enrolled and drop out participants
# ===============================================================

covariates_to_compare <- c("age_group_y1", "BMI_group_y1", "sex", "income", "edu", "asthma", #, "derm"
                           "smoking", "treatment_Y4", "superparticipant_status", 
                           "infection_status_flu", "infection_status_cov2", "cov2_vax_status") 
cov_comp_enrolled_drop <- do.call(rbind, map(
  covariates_to_compare, 
  ~cov_comp_discrete(!!sym(.x), covariate_data_drive1, covariates_for_drop_outs)[[1]]
))

#$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$
#$$$$$$    p-values for the infection and vaccination status   #$$$$$$
#$$$$$$     variables are estimated from the following code    #$$$$$$
#$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$. #$$$$$$

# ================================================================
#        Time varying survival analysis to associate risk of 
#        dropping out with infection and vaccination status
# ================================================================

# Create the static baseline data 
data_SA1 <- timepoint_days_DRIVEI %>% 
  left_join(withdrawal_reasons_DRIVEI %>%
              dplyr::select(pID, wdf_wddate)) %>% 
  filter(!is.na(timepoint_date)) %>% # Ignore NA sample dates to find the first sample
  group_by(pID) %>%
  summarise(
    first_sample_date = min(timepoint_date, na.rm = TRUE),
    drop_out_date = first(as.Date(wdf_wddate)) # Same for all rows of this ID
  ) %>%
  mutate(
    # Determine the exit date
    # If dropped out BEFORE end of study, exit is dropout_date
    # Otherwise, they are censored at end_of_study_date
    exit_date = if_else(
      !is.na(drop_out_date) & drop_out_date < Y4D182_date,
      drop_out_date,
      Y4D182_date),
    # Calculate total observation time in days
    time_to_exit = as.numeric(exit_date - first_sample_date, units = "days")
  ) %>%
  filter(time_to_exit > 0) %>% #Remove anomalies
  left_join(dropout_status %>% 
              rename(status = dropout) %>% 
              mutate(status = as.numeric(status)))  

# total number of dropouts
sum(data_SA1$status)

# Prepare longitudinal time-varying data
data_SA2 <- timepoint_days_DRIVEI %>%
  filter(!is.na(timepoint_date)) %>%
  group_by(pID) %>%
  mutate(
    first_sample_date = min(timepoint_date, na.rm = TRUE),
    # Calculate sample time in days from their first sample
    sample_time = as.numeric(timepoint_date - first_sample_date, units = "days")
  ) %>%
  ungroup() %>%
  dplyr::select(pID, sample_time, cov2_vax_before_sample, infection_before_sample_anyflu_PCR, infection_before_sample_cov2_PCR)


tv_data_SA <- tmerge(
  data1 = data_SA1, 
  data2 = data_SA1, 
  id = pID, 
  dropout = event(time_to_exit, status)
)

# Add the time-dependent covariates
# tdc() ensures the covariate value applies from `sample_time` onward until the next sample
tv_data_SA <- tmerge(
  data1 = tv_data_SA, 
  data2 = data_SA2, 
  id = pID, 
  cov2vax = tdc(sample_time, cov2_vax_before_sample),
  anyflu = tdc(sample_time, infection_before_sample_anyflu_PCR),
  cov2inf = tdc(sample_time, infection_before_sample_cov2_PCR)
)

# ==========================================
# Fit Cox Models (Separately for each factor)
# ==========================================

#cov2vax turns out to be a perfect predictor of drop out risk in this case
model_cov2vax <- tryCatch(
  {# Attempt standard Cox model
    coxph(Surv(tstart, tstop, dropout) ~ cov2vax, data = tv_data_SA)
  },
  warning = function(w) {
    cat("Notice: coxph threw a warning (", w$message, ").\nFalling back to coxphf...\n", sep="")
    coxphf(Surv(tstart, tstop, dropout) ~ cov2vax, data = tv_data_SA)
  },
  error = function(e) {
    cat("Notice: coxph failed (", e$message, ").\nFalling back to coxphf...\n", sep="")
    coxphf(Surv(tstart, tstop, dropout) ~ cov2vax, data = tv_data_SA)
  })

print(summary(model_cov2vax))

model_anyflu <- coxph(Surv(tstart, tstop, dropout) ~ anyflu, data = tv_data_SA)
print(summary(model_anyflu))

model_cov2inf <- coxph(Surv(tstart, tstop, dropout) ~ cov2inf, data = tv_data_SA)
print(summary(model_cov2inf))

#this is what is actually compared in a survival model
dropout_rate_flu <- tv_data_SA %>%
  group_by(anyflu) %>%
  summarise(Total_Dropouts = sum(dropout == 1, na.rm = TRUE),
            Total_Time_Observed = sum(tstop - tstart, na.rm = TRUE)) %>%
  mutate(Dropout_Rate = Total_Dropouts/Total_Time_Observed)

dropout_rate_cov2 <- tv_data_SA %>%
  group_by(cov2inf) %>%
  summarise(Total_Dropouts = sum(dropout == 1, na.rm = TRUE),
            Total_Time_Observed = sum(tstop - tstart, na.rm = TRUE)) %>%
  mutate(Dropout_Rate = Total_Dropouts/Total_Time_Observed)

dropout_rate_cov2vax <- tv_data_SA %>%
  group_by(cov2vax) %>%
  summarise(Total_Dropouts = sum(dropout == 1, na.rm = TRUE),
            Total_Time_Observed = sum(tstop - tstart, na.rm = TRUE)) %>%
  mutate(Dropout_Rate = Total_Dropouts/Total_Time_Observed)

# ======================================================
# Covariate comparison across treatment groups year wise
# ======================================================

covariates_to_compare_yrwise <- c("age_group_y1", "BMI_group_y1", "sex", "asthma", "smoking", #"derm", 
                                  "cov2_vax_before_sample", "infection_before_sample_anyflu_PCR", "infection_before_sample_cov2_PCR")

treatment_yrs <- c("treatment_Y1", "treatment_Y2", "treatment_Y3", "treatment_Y4")

covariates_after_drop_outs <- map(covariates_after_drop_outs, ~ .x %>%
                                    mutate(across(any_of(covariates_to_compare_yrwise), as.character)))

res <- list(); name_df <- c()
for(i in 1:length(covariates_to_compare_yrwise)){
  for(j in 1:length(treatment_yrs)){ #length(treatment_yrs)
    k <- (i - 1)*length(treatment_yrs) + j #length(treatment_yrs)
    name_df[k] <- paste(covariates_to_compare_yrwise[i], treatment_yrs[j], sep = " ")
    res[[k]] <- map(covariates_after_drop_outs[-length(covariates_after_drop_outs)], #-length(covariates_after_drop_outs)
                    ~cov_comp_treat(!!sym(covariates_to_compare_yrwise[i]), !!sym(treatment_yrs[j]), .x)[[1]]
    ) %>% 
      #adding enrollment year as a column to each dataframe
      imap(~{.x %>% mutate(enrolled_in = sub("_.*", "", .y))})
  }
}

names(res) <- name_df

#combine all Y1 dataframes
Y1_covs <- res[grepl("Y1", names(res))]
#combine all Y2 dataframes
Y2_covs <- res[grepl("Y2", names(res))]
#combine all Y3 dataframes
Y3_covs <- res[grepl("Y3", names(res))]
#combine all Y4 dataframes
Y4_covs <- res[grepl("Y4", names(res))]

# Combine selected dataframes by row
combined_df_Y1 <- bind_rows(lapply(Y1_covs, bind_rows)) #bind_rows(Y1_covs); 
combined_df_Y2 <- bind_rows(lapply(Y2_covs, bind_rows))
combined_df_Y3 <- bind_rows(lapply(Y3_covs, bind_rows))
combined_df_Y4 <- bind_rows(lapply(Y4_covs, bind_rows))

final_df_Y1 <- combined_df_Y1 %>% filter(enrolled_in == "Y1")
final_df_Y2 <- combined_df_Y2 %>% filter(enrolled_in == "Y2")
final_df_Y3 <- combined_df_Y3 %>% filter(enrolled_in == "Y3")
final_df_Y4 <- combined_df_Y4 %>% filter(enrolled_in == "Y4")

#Tables S3,S4,S5, and S6
final_df_Y1 <- final_df_Y1 %>% 
  mutate(P = paste0(P, " (",P_pct, "%)"),
         V = paste0(V, " (",V_pct, "%)")) %>% 
  dplyr::select(c(covariate, attribute, P, V, fisher_pval))

final_df_Y2 <- final_df_Y2 %>% 
  mutate(PP = paste0(PP, " (",PP_pct, "%)"),
         PV = paste0(PV, " (",PV_pct, "%)"),
         VV = paste0(VV, " (",VV_pct, "%)")) %>% 
  dplyr::select(c(covariate, attribute, PP, PV, VV, fisher_pval))

final_df_Y3 <- final_df_Y3 %>% 
  mutate(PPP = paste0(PPP, " (",PPP_pct, "%)"),
         PPV = paste0(PPV, " (",PPV_pct, "%)"),
         PVV = paste0(PVV, " (",PVV_pct, "%)"),
         VVV = paste0(VVV, " (",VVV_pct, "%)")) %>% 
  dplyr::select(c(covariate, attribute, PPP, PPV, PVV, VVV, fisher_pval))

final_df_Y4 <- final_df_Y4 %>% 
  mutate(PPPP = paste0(PPPP, " (",PPPP_pct, "%)"),
         PPPV = paste0(PPPV, " (",PPPV_pct, "%)"),
         PPVV = paste0(PPVV, " (",PPVV_pct, "%)"),
         PVVV = paste0(PVVV, " (",PVVV_pct, "%)"),
         VVVV = paste0(VVVV, " (",VVVV_pct, "%)")) %>% 
  dplyr::select(c(covariate, attribute, PPPP, PPPV, PPVV, PVVV, VVVV, fisher_pval))

output_dir <- "results/covariate_and_dropout_analysis"
dir.create(output_dir, showWarnings = F)

write.csv(cov_comp_enrolled_drop, paste0(output_dir, "/dropout_characteristics.csv"), row.names = F)
write.csv(final_df_Y1, paste0(output_dir, "/group_characteristics_Y1.csv"), row.names = F)
write.csv(final_df_Y2, paste0(output_dir, "/group_characteristics_Y2.csv"), row.names = F)
write.csv(final_df_Y3, paste0(output_dir, "/group_characteristics_Y3.csv"), row.names = F)
write.csv(final_df_Y4, paste0(output_dir, "/group_characteristics_Y4.csv"), row.names = F)