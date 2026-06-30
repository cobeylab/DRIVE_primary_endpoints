library(drc)
library(tidyverse)
library(lubridate)
library(ggtext)
library(patchwork)

# Bulk of code by A.S.

# Requires `serum_samples`, `HK_flu_surveillance` `positive_pcr_tests` from load data
source("R/load_data.R")

#function to fit 2-parameter Log-logistic model
fit_safe.ll2 <- function(dat) {
  tryCatch(
    drm(formula = MeanNAAct ~ NAAct_dil, data = dat, fct = LL.2(upper = 1)), 
    error = function(e) NULL   
  )
}

#function to process raw NAI samples for consistency in column names and get long-data format to fit Log-logistic models
process_NAI_data <- function(ella_assay_data, serum_samples_data, assay_type){
  
  #rename columns all uniformly across data.frames
  ella_assay_data <- ella_assay_data %>% 
    dplyr::select(c(`Sample ID`, matches("Extrapolated Titer.*log2x"), 
                    matches("Titer.*50%"), 
                    #`y1 = NA @ x1`, `y2 = NA @ x2`
                    starts_with("Mean"))) %>% 
    dplyr::select(!starts_with("Fold Change")) %>% 
    rename(sid = `Sample ID`)
  
  if (assay_type == "limited") {
    ella_assay_data <- ella_assay_data %>% 
      rename_with(.cols = matches("Titer.*<"),
                  .fn = ~ paste0("x2 = Titer @ NA < 50%")) %>% 
      rename_with(.cols = matches("Titer.*>"),
                  .fn = ~ paste0("x1 = Titer @ NA > 50%")) %>% 
      rename_with(.cols = matches("Extrapolated Titer.*log2x"),
                  .fn = ~ paste0("Extrapolated Titer [log2x vs y]"))
  } else {
    ella_assay_data <- ella_assay_data
  }
  
  ella_assay_data <- ella_assay_data %>% 
    filter(sid != "Pooled Sera") %>% 
    mutate(across(
      .cols = starts_with("MeanAdjOD"),
      #~100* (`Mean Adjusted Virus OD` - .x)/`Mean Adjusted Virus OD`,
      .fns  = ~if_else(.x < `Mean Adjusted Virus OD`, .x/`Mean Adjusted Virus OD`, 1), 
      .names = "New_MeanNAAct_{.col}"
    ))
  
  serum_samples_sub <- serum_samples_data %>% filter(sid %in% ella_assay_data$sid)
  
  ella_assay_data <- left_join(ella_assay_data, 
                               serum_samples_sub %>% 
                                 dplyr::select(c(sid, pID, year, 
                                                 treatment, timepoint, timepoint_date,
                                                 most_recent_infection_flu_A_PCR,
                                                 most_recent_infection_H3N2_PCR,
                                                 most_recent_infection_H1N1_PCR))) 
  
  ella_assay_data <- ella_assay_data 
  
  if (assay_type == "limited"){
    cols_to_select <- c("New_MeanNAAct_MeanAdjOD_1", "New_MeanNAAct_MeanAdjOD_2", 
                        "New_MeanNAAct_MeanAdjOD_3", "New_MeanNAAct_MeanAdjOD_4")
  } else {
    cols_to_select <- c("New_MeanNAAct_MeanAdjOD_1", "New_MeanNAAct_MeanAdjOD_2", 
                        "New_MeanNAAct_MeanAdjOD_3", "New_MeanNAAct_MeanAdjOD_4",
                        "New_MeanNAAct_MeanAdjOD_5", "New_MeanNAAct_MeanAdjOD_6", 
                        "New_MeanNAAct_MeanAdjOD_7", "New_MeanNAAct_MeanAdjOD_8",
                        "New_MeanNAAct_MeanAdjOD_9", "New_MeanNAAct_MeanAdjOD_10")
  }
  
  ella_assay_data_long <- ella_assay_data %>% 
    pivot_longer(cols = all_of(cols_to_select), 
                 names_to = "NAAct_dil", values_to = "MeanNAAct")
  
  if (assay_type == "limited"){
    ella_assay_data_long <- ella_assay_data_long %>% 
      mutate(NAAct_dil = case_when(
        NAAct_dil == "New_MeanNAAct_MeanAdjOD_1" ~ 20,
        NAAct_dil == "New_MeanNAAct_MeanAdjOD_2" ~ 80,
        NAAct_dil == "New_MeanNAAct_MeanAdjOD_3" ~ 320,
        NAAct_dil == "New_MeanNAAct_MeanAdjOD_4" ~ 1280,
      ))} else {
        ella_assay_data_long <- ella_assay_data_long %>% 
          mutate(NAAct_dil = case_when(
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_1" ~ 20,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_2" ~ 40,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_3" ~ 80,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_4" ~ 160,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_5" ~ 320,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_6" ~ 640,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_7" ~ 1280,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_8" ~ 2560,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_9" ~ 5120,
            NAAct_dil == "New_MeanNAAct_MeanAdjOD_10" ~ 10240
          ))
      }
  
  ella_assay_data_long <- ella_assay_data_long %>% 
    dplyr::select(-starts_with("^MeanNAAct_"))
  
  return(ella_assay_data_long)
}


#function to fit Log-logistic models on processed NAI data to estimate NAI titer
fit_log_logistic_models <- function(NAI_data_long, half_lowest_dil){
  
  dfs_to_fit_group <- NAI_data_long %>%
    group_by(sid) 
  dfs_to_fit <- dfs_to_fit_group %>% group_split()
  
  names(dfs_to_fit) <- dfs_to_fit_group %>% group_keys() %>% pull(sid)

  # drm()'s internal optim() writes convergence-failure text directly to stderr,
  # bypassing R's condition system (so suppressWarnings/tryCatch can't silence it).
  # Redirect the message connection to a temp file while fitting, then restore it.
  err_con <- file(tempfile(), open = "wt")
  sink(err_con, type = "message")
  fits.ll2 <- map(dfs_to_fit, fit_safe.ll2)
  sink(type = "message")
  close(err_con)
  
  conv_df <- tibble(
    sid  = names(fits.ll2),
    converged.ll2  = !map_lgl(fits.ll2, is.null)
  )
  
  NAI_data_long_up <- NAI_data_long %>%
    left_join(conv_df, by = "sid")
  
  coef_df <- tibble(
    sid = names(fits.ll2),
    estimated_titer = map_dbl(fits.ll2, ~ {
      cc <- try(coef(.x), silent = TRUE)
      if (inherits(cc, "try-error") || length(cc) < 2) return(NA) #for the fits that did not converge
      unname(cc[2])
    }),
    coef1_ll2 = map_dbl(fits.ll2, ~ {
      cc <- try(coef(.x), silent = TRUE) 
      if (inherits(cc, "try-error") || length(cc) < 2) return(NA) #for the fits that did not converge
      unname(cc[1])
    })
  )
  
  #handling cases where the logistic regression estimated slope is negative; 
  #i.e., ensuring fitting a upward trend in NA activity vs serum dilution
  coef_df <- coef_df %>% 
    mutate(estimated_titer = ifelse(coef1_ll2 > 0, half_lowest_dil, estimated_titer))
  
  NAI_data_long_up <- NAI_data_long_up %>%
    left_join(coef_df, by = "sid")
  
  NAI_data_long_up <- NAI_data_long_up %>% 
    mutate(estimated_titer = ifelse(converged.ll2 == F, half_lowest_dil, estimated_titer))
  
  #prediction from ll.2 fit
  pred_grid.ll2 <- NAI_data_long_up %>%
    filter(converged.ll2) %>%
    group_by(sid) %>%
    reframe(
      NAAct_dil = seq(1, 10240, length.out = 1000)
    ) %>%
    mutate(
      pred = map2_dbl(sid, NAAct_dil, ~ {
        fit <- fits.ll2[[as.character(.x)]]
        as.numeric(predict(fit, newdata = data.frame(NAAct_dil = .y)))
      })
    )
  
  return(list(coef_df, pred_grid.ll2, NAI_data_long_up))
  
}

#function to calculate fold-change in NAI titer between consecutive time points
calculate_fold_change <- function(df) {
  if(nrow(df) > 1){
    
    df <- df %>% arrange(year, timepoint, .by_group = T) 
    
    #initializing
    df$fc_from_prev_ll2 <- 0
    df$time_bet_sample <- 0
    df$baseline_titer <- 1
    
    #calculating actual values
    for (i in 2:nrow(df)){
      df$baseline_titer[i] <- df$estimated_titer[i-1]
      df$fc_from_prev_ll2[i] <- df$estimated_titer[i]/df$estimated_titer[i-1]
      df$time_bet_sample[i] <- as.numeric(df$timepoint_date[i] - df$timepoint_date[i-1])
    }
    return(df)
    
  }else {
    return(data.frame(timepoint_date = NA))
  }
}

#function to get the NAI titers and fold-change in titer of PCR confirmed infections before and after the infection dates
NAI_pcr_confirmed <- function(fold_change_df, N1_or_N2){
  
  #this is all of titer_df
  infection_dates_for_pcr_confirmed_infections <- fold_change_df %>% 
    filter(!is.na(vline_x)) %>%  #vline_x represents the date of PCR confirmed infections
    group_by(pID) %>% 
    summarise(is_inf_fluA = any(unique(!is.na(most_recent_infection_flu_A_PCR)) == 1),
              is_inf_h1n1 = any(unique(!is.na(most_recent_infection_H1N1_PCR)) == 1),
              is_inf_h3n2 = any(unique(!is.na(most_recent_infection_H3N2_PCR)) == 1)) %>% 
    ungroup()
  
  #for "D10365", NAI titers were measured only for dates before they were PCR-positive 
  
  if(N1_or_N2 == "N1"){
    pcr_identified_infections <- (infection_dates_for_pcr_confirmed_infections %>% 
                                       filter(is_inf_h1n1 == T))$pID #only H1N1 positives
  } else {
    pcr_identified_infections <- (infection_dates_for_pcr_confirmed_infections %>% 
                                       filter(is_inf_h3n2 == T))$pID #only H3N2 positives
  }
    NAI_of_pcrs <- fold_change_df %>% 
      filter(pID %in% pcr_identified_infections) %>% 
      group_by(pID)  %>% 
      #arrange(timepoint_date, .by_group = T) %>% 
      arrange(year, timepoint, .by_group = T) %>% 
      mutate(pcr_infection_date = unique(vline_x[!is.na(vline_x)]),
             before_infection_date = timepoint_date[findInterval(pcr_infection_date, timepoint_date)],
             after_infection_date = timepoint_date[findInterval(pcr_infection_date, timepoint_date) + 1],
             titer_before_infection = estimated_titer[timepoint_date == before_infection_date],
             titer_after_infection = estimated_titer[timepoint_date == after_infection_date],
             fold_change_during_infection = titer_after_infection/titer_before_infection) 
    
    if(N1_or_N2 == "N1"){ 
      NAI_of_pcrs <- NAI_of_pcrs %>% mutate(antigen = "N1")
    } else {
      NAI_of_pcrs <- NAI_of_pcrs %>% mutate(antigen = "N2")
    }
  
  return(NAI_of_pcrs)
  
}

#function to estimate the threshold fold change to call for a boost in NAI titer based on the antibody ceiling line
adaptive_threshold_to_call_for_infection <- function(df, ceiling_line_intercept, ceiling_line_slope){
  
  #initializing
  if(nrow(df) > 1){
    #initializing
    df <- df %>% arrange(year, timepoint, .by_group = T)
    df$threshold <- 1
    
    #calculating actual values
    for (i in 2:nrow(df)){
      #the minimum threshold to detect an infection is 1.5 (2^0.585) fold boost in titer
      df$threshold[i] <- 2^max(log2(min_threshold),(ceiling_line_intercept + ceiling_line_slope*log(df$baseline_titer[i], base = 2)))
    }
    df$boost_ll2 <- ifelse(df$fc_from_prev_ll2 >= df$threshold, 1,0)
    df$boost_ll2[is.na(df$boost_ll2)] <- 0
    
    return(df)
  }else {return(data.frame())}
}

# Based on the titer fold change df, this function first computes the threshold to call for an infection for each pID
# Then for each possibly infected individual, this function gathers the bracketing dates and titers on those dates corresponding to a boost 
# Finally it calculates the possible infection date for each individual based on the surveillance data and computes an interpolated titer
# This function also outputs the true positive and false negative infections from NAI data of PCR confirmed infections
NAI_sensitivity <- function(fold_change_df, NAI_pcr_confirmed_df, N1_or_N2, threshold_boost_intercept){ 
  #this will get rid off the pIDs for whom NAI titer was measured only once
  titer_fold_change <- fold_change_df %>% 
    group_by(pID) %>% 
    group_modify(~adaptive_threshold_to_call_for_infection(.x, threshold_boost_intercept, ceiling_line_slope)) %>% 
    ungroup()
  
  ids_with_boost <- unique(titer_fold_change$pID[titer_fold_change$boost_ll2 == 1])
  identify_infections_from_titer <- titer_fold_change %>% filter(pID %in% ids_with_boost)
  identify_noinfections_from_titer <- titer_fold_change %>% filter(!pID %in% ids_with_boost)
  
  time_of_infection <- identify_infections_from_titer %>% 
    dplyr::select(pID, year, treatment, timepoint, timepoint_date, 
                  most_recent_infection_flu_A_PCR, most_recent_infection_H3N2_PCR, most_recent_infection_H1N1_PCR,
                  estimated_titer, vline_x, fc_from_prev_ll2, time_bet_sample, baseline_titer, threshold, boost_ll2) %>% 
    group_by(pID) %>% 
    #I have to do this step as there is an NA timepoint_date during a boost in NAI titer
    group_modify(~{.x %>%
        mutate(timepoint_date = as.Date(ifelse(is.na(timepoint_date), lag(timepoint_date) + timepoint, timepoint_date)))
    }) %>%
    mutate(before_boost_index = c(boost_ll2[-1] == 1 & boost_ll2[-length(boost_ll2)] == 0 |
                                    boost_ll2[-1] == 1 & boost_ll2[-length(boost_ll2)] == 1, F),
           before_boost_date = ifelse(before_boost_index, timepoint_date, NA),
           after_boost_date = c(ifelse(boost_ll2 == 1, timepoint_date, NA)[-1], NA),
           titer_before_boost = ifelse(before_boost_index, estimated_titer, NA),
           titer_after_boost = c(ifelse(boost_ll2 == 1, estimated_titer, NA)[-1], NA)
    ) %>% 
    ungroup()
  
  ids_with_post_boost_titer_at_half_lowest_dil <- time_of_infection$pID[time_of_infection$titer_after_boost == half_lowest_dil_value] %>% na.omit()
  
  time_of_infection <- time_of_infection %>% 
    filter(!pID %in% ids_with_post_boost_titer_at_half_lowest_dil) 
  
  #infection date is inferred based on the surveilance data
  if(N1_or_N2 == "N1"){
    infection_date_vec <- as.Date(as.numeric(t(rbind(mapply(
      infer_infection_date,
      interval_start_date = as.Date(time_of_infection$before_boost_date),
      interval_end_date = as.Date(time_of_infection$after_boost_date),
      MoreArgs = list(surveillance_data = HK_surveilance_upd, N1_or_N2 = "N1"))
    ))))} else {
      infection_date_vec <- as.Date(as.numeric(t(rbind(mapply(
        infer_infection_date,
        interval_start_date = as.Date(time_of_infection$before_boost_date),
        interval_end_date = as.Date(time_of_infection$after_boost_date),
        MoreArgs = list(surveillance_data = HK_surveilance_upd, N1_or_N2 = "N2"))
      ))))
    }
  
  time_of_infection <- time_of_infection %>% 
    mutate(infection_date = infection_date_vec) 
  
  time_of_infection <- interpolate_titer_at_infection_date(time_of_infection) %>% 
    mutate(threshold_intercept = threshold_boost_intercept)
  
  #for each possible infected individual, counting the number of times they were infected
  number_infections <- time_of_infection %>% 
    group_by(pID) %>% 
    summarise(number_of_infections = length(which(boost_ll2 == 1))) %>% 
    mutate(threshold_intercept = threshold_boost_intercept) %>% 
    ungroup()
  
  #for "D10365", NAI titers were measured only for dates before they were PCR-positive 
  NAI_pcr_infs <- NAI_pcr_confirmed_df %>% 
    group_by(pID) %>% 
    mutate(threshold_during_infection = 2^max(log2(min_threshold),(threshold_boost_intercept + ceiling_line_slope*log(titer_before_infection, base = 2))),
           boost_during_infection = ifelse(fold_change_during_infection >= threshold_during_infection, 1, 0)) %>% 
    ungroup()
  
  true_positive_infections <- unique(NAI_pcr_infs$pID[NAI_pcr_infs$boost_during_infection == 1])
  false_negative_infections <- unique(NAI_pcr_infs$pID[NAI_pcr_infs$boost_during_infection == 0])
  
  true_positive_numbers <- length(true_positive_infections)
  false_negative_numbers <- length(false_negative_infections)
  
  df_return <- data.frame(
    antigen = N1_or_N2,
    num_true_positive = true_positive_numbers,
    num_false_negative = false_negative_numbers,
    total_NAI_based_infections = sum(number_infections$number_of_infections),
    median_threshold = median(titer_fold_change$threshold, na.rm = T),
    threshold_intercept = threshold_boost_intercept
  )
  
  pIDs_possible_infections <- time_of_infection %>% 
    filter(!is.na(infection_date)) %>% 
    dplyr::select(pID, infection_date)
  
  
  return(list(df_return, #summary table of TP and FN infections to estimate sensitivity of the assay
              time_of_infection, 
              pIDs_possible_infections, #IDs of participants with possible infection based on the threshold
              number_infections, #How many times each participant got infected according to NAI boost
              titer_fold_change, 
              true_positive_infections, #IDs of TP infections
              false_negative_infections #IDs of FN infections
  )) 
  
}


week_start <- function(d) ceiling_date(d, unit = "week", week_start = 7) - 7 #week starts on Sunday and ends on Saturday
week_end <- function(d) ceiling_date(d, unit = "week", week_start = 7) - 1  #week starts on Sunday and ends on Saturday


infer_infection_date <- function(surveillance_data, interval_start_date, interval_end_date, N1_or_N2){
  if (is.na(interval_start_date)) return(NA)
  
  interval_start_date <- week_start(interval_start_date) %>% na.omit()
  interval_end_date <- week_end(interval_end_date) %>% na.omit()
  
  case_col <- if (N1_or_N2 == "N1") "H1_pos_1000visits" else "H3_pos_1000visits"
  
  inferred_infection_date <- surveillance_data %>%
    filter(From >= interval_start_date, To <= interval_end_date) %>%
    rename(cases = all_of(case_col)) %>%
    arrange(From) %>%
    mutate(cumulative_cases = cumsum(cases),
           total_cases = sum(cases, na.rm = TRUE),
           frac = cumulative_cases / total_cases) %>%
    filter(frac >= 0.5) %>%      # crossed 50%
    slice(1) %>% 
    mutate(infection_date = as.Date(0.5*(as.numeric(From) + as.numeric(To)))) %>%
    pull(infection_date)
  
  return(inferred_infection_date)
}

interpolate_titer_at_infection_date <- function(df) {
  df <- df %>%
    mutate(frac = as.numeric(infection_date - as.Date(before_boost_date)) / as.numeric(as.Date(after_boost_date) - as.Date(before_boost_date)),
           titer_est = 2^(log(titer_before_boost, base = 2) + frac*(log(titer_after_boost, base = 2) - log(titer_before_boost, base = 2))))
  #pull(titer_est)
}


get_low_circ_intervals <- function(surveillance_data, subtype_col, cutoff_low_circ){
  
  val_round <- round(surveillance_data[[subtype_col]], 2)
  
  df_low_circ <- surveillance_data[val_round <= cutoff_low_circ,
                                   c("Year", "Week", "From", "To", subtype_col)][-c(1:40), ] #the earliest NAI titer was measured on 9th Oct, 2022 (Week 41)
  
  total_burden_during_low_circ <- sum(df_low_circ[[subtype_col]])
  
  df_high_circ <- surveillance_data[val_round > cutoff_low_circ,
                                    c("Year", "Week", "From", "To", subtype_col)] #the earliest NAI titer was measured on 9th Oct, 2022 (Week 41)
  
  total_burden_during_high_circ <- sum(df_high_circ[[subtype_col]])
  
  if (nrow(df_low_circ) == 0) return(df_low_circ[F, ])
  
  #sort by Year, Week
  df_low_circ <- df_low_circ[order(df_low_circ$Year, df_low_circ$Week), ]
  
  idx <- df_low_circ$Year*52 + df_low_circ$Week
  
  #new group whenever this index is not consecutive
  new_group <- c(T, diff(idx) > 1)
  grp <- cumsum(new_group)
  
  # for each group, take first/last row as interval boundaries
  df_return <- do.call(rbind, lapply(split(df_low_circ, grp), function(d){
    data.frame(Year_start = d$Year[1],
               Year_end = d$Year[nrow(d)],
               Week_start = d$Week[1],
               Week_end = d$Week[nrow(d)],
               Week_start_date = d$From[1],
               Week_end_date = d$To[nrow(d)],
               stringsAsFactors = F) %>% 
      mutate(interval_length = as.numeric(Week_end_date - Week_start_date)) %>% 
      filter(interval_length >= 14) #low circulation period should be at least 2 weeks long
  }
  ))
  
  rownames(df_return) <- NULL
  
  df_return %>%
    mutate(ymin = -Inf,
           ymax = Inf)
  
  return(list(df_return, total_burden_during_low_circ, total_burden_during_high_circ, df_low_circ, df_high_circ))
}

#This function first chops up the titer data based on the identified low influenza circulation levels generated from the earlier function
#Then it identifies how many NAI titer boost is identified within low circulation periods
NAI_specificity <- function(titer_df, low_circ_df, threshold_boost_intercept, inferred_inf_pIDs){
  
  d <- titer_df[["timepoint_date"]]
  
  #dividing titer data into low circulation intervals
  titer_dfs_low_circ <- lapply(seq_len(nrow(low_circ_df)), function(i) {
    from_i <- low_circ_df$Week_start_date[i]
    to_i   <- low_circ_df$Week_end_date[i]
    
    titer_df[d >= from_i & d <= to_i, ,drop = F] %>% filter(!is.na(sid))
  })
  
  #calculating titer fold change and any boost during low circulation intervals
  titer_fold_change_low_circ_adapt <- do.call(rbind, 
                                              map(titer_dfs_low_circ, 
                                                  ~(.x %>% 
                                                      group_by(pID) %>% 
                                                      group_modify(~calculate_fold_change(.x)) %>% 
                                                      group_modify(~adaptive_threshold_to_call_for_infection(.x, threshold_boost_intercept, ceiling_line_slope)) %>% 
                                                      ungroup())))
  
  find_boost_low_circ <- titer_fold_change_low_circ_adapt %>% 
    filter(time_bet_sample != 0) %>% 
    group_by(pID) %>% 
    mutate(any_boost = any(boost_ll2 == 1))
  
  #False positives are only among the identified infected ones coming from the sensitivity analysis
  df_return <- data.frame(num_false_positive = sum(find_boost_low_circ$any_boost[find_boost_low_circ$pID %in% inferred_inf_pIDs$pID]),
                          num_true_negative = sum(find_boost_low_circ$any_boost == F),
                          threshold_intercept = threshold_boost_intercept)
  
  return(df_return)
  
}

# Abbreviated-month x-axis labels: show the month (Jan, Mar, ...) on every break,
# adding the year on a second line only for January so it appears once per year.
month_year_labels <- function(x) {
  m <- format(x, "%b")
  ifelse(month(x) == 1, paste0(m, "\n", year(x)), m)
}

# Quarterly breaks anchored to January, so each year's January is guaranteed to be
# a break (and thus carry a year label). `limits` is the axis date range.
quarterly_jan_breaks <- function(limits) {
  start <- floor_date(min(limits), "year")
  end   <- ceiling_date(max(limits), "year")
  seq(start, end, by = "3 months")
}

# Build the inferred-infection histogram and infection-date/surveillance overlay
# for a single subtype, then combine them into the stacked figure. H1 and H3
# differ only in their data, color key, surveillance column and y-axis label.
make_inferred_infections_fig <- function(low_circ, time_of_infection_sub,
                                         color_key, surv_subtype, y_axis_name) {

  n_infections_by_month <- time_of_infection_sub %>%
    select(infection_date) %>%
    filter(!is.na(infection_date)) %>%
    # Recode all infection dates as the first day of the corresponding year/month
    mutate(infection_date = paste(year(infection_date), month(infection_date), "1", sep = "-") %>% ymd()) %>%
    group_by(infection_date) %>%
    count() %>%
    ungroup()

  inferred_inf_dates <- ggplot() +
    geom_rect(data = low_circ,
              aes(xmin = Week_start_date, xmax = Week_end_date, ymin = -Inf, ymax = Inf),
              fill = "grey40",
              alpha = 0.2) +
    geom_col(aes(x = infection_date, y = n), data = n_infections_by_month, fill = subtype_colors[color_key], color = "black", linewidth = 0.2) +
    scale_x_date(breaks = quarterly_jan_breaks,
                 labels = month_year_labels,
                 minor_breaks = NULL) +
    labs(x = element_blank(), y = "Number of inferred infections") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_line(linewidth = 0.2))

  date_range <- ggplot_build(inferred_inf_dates)$layout$panel_scales_x[[1]]$range$range

  infection_date_surv <- ggplot() +
    geom_col(data = HK_surveilance_upd_long %>%
               filter(Year %in% c(2022, 2023, 2024, 2025) & subtype == surv_subtype),
             aes(x = monday_date, y = num_pos_cases_1000visitis), fill = subtype_colors[color_key], color = "black", linewidth = 0.2) +
    geom_rect(data = low_circ,
              aes(xmin = Week_start_date, xmax = Week_end_date, ymin = -Inf, ymax = Inf),
              fill = "grey40",
              inherit.aes = FALSE,
              alpha = 0.2) +
    geom_point(data = time_of_infection_sub,
               aes(timepoint_date, log(estimated_titer, base = 2)*scale_factor, group = pID),
               color = subtype_colors[color_key], alpha = 0.6, size = point_size) +
    geom_line(data = time_of_infection_sub,
              aes(timepoint_date, log(estimated_titer, base = 2)*scale_factor, group = pID),
              color = subtype_colors[color_key], alpha = 0.6, size = 0.3*line_size) +
    geom_point(data = time_of_infection_sub,
               aes(x = as.Date(as.numeric(infection_date)), y = log(titer_est, base = 2)*scale_factor),
               size = point_size, alpha = 0.6) +
    theme(axis.text.x = element_text(size = axis_text)) +
    scale_x_date(
      breaks = quarterly_jan_breaks,
      labels = month_year_labels,
      minor_breaks = NULL,
      limits = as.Date(date_range)
    ) +
    scale_y_continuous(
      name = y_axis_name,
      sec.axis = sec_axis(~ . / scale_factor, name = expression(log[2]("NAI titer")))) +
    theme_bw() +
    xlab("Year") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank())

  (infection_date_surv / inferred_inf_dates) +
    plot_layout(axis_titles = "collect") +
    labs(x = "Year") +  # shared titles
    plot_annotation(tag_levels = "A") &
    theme(
      axis.title = element_text(size = axis_label),
      axis.text = element_text(size = axis_text),
      axis.title.y = element_text(size = axis_label),
      plot.tag = element_text(margin = margin(t = 10, r = 5, b = 18, l = 2), face = "bold", size = 10),
      plot.tag.position = c(0, 1)
    )
}

#parameters for this analysis
low_circ_threshold <- 0.5
min_threshold <- 1.5
half_lowest_dil_value <- 10 #lowest dilution tested in both full and limited dilution assay is 20

#NAI data file path
NAI_titers_Batch1_file <- file.path(data_folder_path, "NAI_titers/DRIVE NAI ELLA Result Sheet_20251208 (Batch 1).xlsx")
NAI_titers_Batch2_file <- file.path(data_folder_path, "NAI_titers/DRIVE NAI ELLA Result Sheet_20251208 (Batch 2).xlsx")
NAI_titers_Batch3_file <- file.path(data_folder_path, "NAI_titers/DRIVE NAI ELLA Result Sheet_20251208 (Batch 3).xlsx")
#NAI data full and limited dilution file paths
NAI_titers_full_dil_file <- file.path(data_folder_path, "NAI_titers/DRIVE NAI ELLA Result Sheet_20251210 (PCR-Positives) (Full Dilution).xlsx")
NAI_titers_lim_dil_file <- file.path(data_folder_path, "NAI_titers/DRIVE NAI ELLA Result Sheet_20251210 (PCR-Positives) (Limited Dilution).xlsx")

# HK_surveillance loaded by load_data.R
HK_surveilance_upd <- HK_flu_surveillance %>%
  filter(Year %in% c(2022, 2023, 2024, 2025))

HK_surveilance_upd_long <- HK_surveilance_upd %>% 
  dplyr::select(c(Year, Week, H1_pos_1000visits, H3_pos_1000visits)) %>% 
  pivot_longer(cols = c(H1_pos_1000visits, H3_pos_1000visits),
               names_to = "subtype",
               values_to = "num_pos_cases_1000visitis") %>% 
  mutate(monday_date = ISOweek::ISOweek2date(paste0(Year, "-W", sprintf("%02d", Week), "-1")))


#gives overall intervals of low circulation periods by compressing consecutive weeks
H1_low_circ <- get_low_circ_intervals(HK_surveilance_upd, "H1_pos_1000visits", low_circ_threshold)[[1]]
H3_low_circ <- get_low_circ_intervals(HK_surveilance_upd, "H3_pos_1000visits", low_circ_threshold)[[1]]

# ======================================
#     Load N1 data from the file
# ======================================

N1_titer_Batch1 <- read_excel(NAI_titers_Batch1_file, sheet = "A Victoria 2570 2019 (N1)")
N1_titer_Batch2 <- read_excel(NAI_titers_Batch2_file, sheet = "A Victoria 2570 2019 (N1)")
N1_titer_Batch3 <- read_excel(NAI_titers_Batch3_file, sheet = "A Victoria 2570 2019 (N1)")
N1_full_dil     <- read_excel(NAI_titers_full_dil_file, sheet = "A Victoria 2570 2019 (N1)") 
N1_lim_dil      <- read_excel(NAI_titers_lim_dil_file, sheet = "A Victoria 2570 2019 (N1)") 

# ======================================
#     Load N2 data from the file
# ======================================

N2_titer_Batch1 <- read_excel(NAI_titers_Batch1_file, sheet = "A Cambodia e0826360 2020 (N2)")
N2_titer_Batch2 <- read_excel(NAI_titers_Batch2_file, sheet = "A Cambodia e0826360 2020 (N2)")
N2_titer_Batch3 <- read_excel(NAI_titers_Batch3_file, sheet = "A Cambodia e0826360 2020 (N2)")
N2_full_dil     <- read_excel(NAI_titers_full_dil_file, sheet = "A Cambodia e0826360 2020 (N2)") 
N2_lim_dil      <- read_excel(NAI_titers_lim_dil_file, sheet = "A Cambodia e0826360 2020 (N2)") 

# ======================================
#       Process loaded N1 data
# ======================================

N1_titer_Batch1p <- process_NAI_data(N1_titer_Batch1, serum_samples, "limited")
N1_titer_Batch2p <- process_NAI_data(N1_titer_Batch2, serum_samples, "limited")
N1_titer_Batch3p <- process_NAI_data(N1_titer_Batch3, serum_samples, "limited")
N1_titer_full_dilp <- process_NAI_data(N1_full_dil, serum_samples, "full")
N1_titer_lim_dilp <- process_NAI_data(N1_lim_dil, serum_samples, "limited")

# ======================================
#       Process loaded N2 data
# ======================================
#for batch 1 N2 data, we need to calculate the MeanAdjODs
N2_titer_Batch1$MeanAdjOD_1 <- rowMeans(N2_titer_Batch1[, c("AdjOD_A1", "AdjOD_B1")], na.rm = T)
N2_titer_Batch1$MeanAdjOD_2 <- rowMeans(N2_titer_Batch1[, c("AdjOD_A2", "AdjOD_B2")], na.rm = T)
N2_titer_Batch1$MeanAdjOD_3 <- rowMeans(N2_titer_Batch1[, c("AdjOD_A3", "AdjOD_B3")], na.rm = T)
N2_titer_Batch1$MeanAdjOD_4 <- rowMeans(N2_titer_Batch1[, c("AdjOD_A4", "AdjOD_B4")], na.rm = T)

N2_titer_Batch1 <- N2_titer_Batch1 %>% 
  rename(`Mean Adjusted Virus OD` = MeanAdjVirusOD)

N2_titer_Batch1p <- process_NAI_data(N2_titer_Batch1, serum_samples, "limited")
N2_titer_Batch2p <- process_NAI_data(N2_titer_Batch2, serum_samples, "limited")
N2_titer_Batch3p <- process_NAI_data(N2_titer_Batch3, serum_samples, "limited")
N2_titer_full_dilp <- process_NAI_data(N2_full_dil, serum_samples, "full")
N2_titer_lim_dilp <- process_NAI_data(N2_lim_dil, serum_samples, "limited")

# ==================================================
#         Estimate NAI titer by fitting 
#  log.logistic regression on processed N1 data
# ==================================================

#Fitting log-logistic model will print error message for the non-convergence cases
#We are bypassing this error by using try-catch and assigning NAI titer = 10 (half of lowest dilution tested)

processed_N1_b1 <- fit_log_logistic_models(N1_titer_Batch1p, half_lowest_dil_value)
processed_N1_b2 <- fit_log_logistic_models(N1_titer_Batch2p, half_lowest_dil_value)
processed_N1_b3 <- fit_log_logistic_models(N1_titer_Batch3p, half_lowest_dil_value)
processed_N1_full_dil <- fit_log_logistic_models(N1_titer_full_dilp, half_lowest_dil_value)
processed_N1_lim_dil <- fit_log_logistic_models(N1_titer_lim_dilp, half_lowest_dil_value)

# ==================================================
#         Estimate NAI titer by fitting 
#  log.logistic regression on processed N2 data
# ==================================================

processed_N2_b1 <- fit_log_logistic_models(N2_titer_Batch1p, half_lowest_dil_value)
processed_N2_b2 <- fit_log_logistic_models(N2_titer_Batch2p, half_lowest_dil_value)
processed_N2_b3 <- fit_log_logistic_models(N2_titer_Batch3p, half_lowest_dil_value)
processed_N2_full_dil <- fit_log_logistic_models(N2_titer_full_dilp, half_lowest_dil_value)
processed_N2_lim_dil <- fit_log_logistic_models(N2_titer_lim_dilp, half_lowest_dil_value)


N1_all_batches <- list(processed_N1_b1, processed_N1_b2, processed_N1_b3, processed_N1_lim_dil)
N1_titer_updated_dfs <- lapply(N1_all_batches, `[[`, 3)   
N1_titer_updated_combined_df <- do.call(rbind, N1_titer_updated_dfs) %>% 
  mutate(year_timepoint = paste("year", year, "day", timepoint)) #this has all serum_dilutions for each `sid`

N1_titer_updated_combined_df_clean <- N1_titer_updated_combined_df %>% 
  distinct(pID, year_timepoint, .keep_all = T) #no of rows = 1501, 312 unique pIDs

N2_all_batches <- list(processed_N2_b1, processed_N2_b2, processed_N2_b3, processed_N2_lim_dil) 
N2_titer_updated_dfs <- lapply(N2_all_batches, `[[`, 3)   
N2_titer_updated_combined_df <- rbind(rbind(N2_titer_updated_dfs[[1]], N2_titer_updated_dfs[[2]][,-c(9:12)]), 
                                            N2_titer_updated_dfs[[3]][,-c(9:12)], N2_titer_updated_dfs[[4]][,-c(9:12)]) %>% 
  mutate(year_timepoint = paste("year", year, "day", timepoint))

N2_titer_updated_combined_df_clean <- N2_titer_updated_combined_df %>% 
  distinct(pID, year_timepoint, .keep_all = T) #no of rows = 1508, 314 unique pIDs

N1_titer_updated_combined_df_clean <- N1_titer_updated_combined_df_clean %>%
  mutate(vline_x = coalesce(most_recent_infection_flu_A_PCR, most_recent_infection_H1N1_PCR))

N2_titer_updated_combined_df_clean <- N2_titer_updated_combined_df_clean %>%
  mutate(vline_x = coalesce(most_recent_infection_flu_A_PCR, most_recent_infection_H3N2_PCR))

# ==================================================
#        Details of PCR confirmed infections
# ==================================================

fluA_positive <- subset(positive_pcr_tests, result %in% c("flu_A_positive", "H3N2_positive", "H1N1_positive")) %>% 
  distinct(pID, .keep_all = T)

subtyped_unsubtyped_pcr <- subset(positive_pcr_tests, result %in% c("flu_A_positive", "H3N2_positive", "H1N1_positive")) %>% 
  group_by(pID) %>% 
  summarise(n_occur = n(),
            distinct_vals = paste(sort(unique(result)), collapse = ", ")) %>% 
  mutate(subtype = str_split(distinct_vals, pattern = ", ", simplify = T)[,2])


pIDs_with_NAI_titer_measurements <- unique(c(unique(N2_titer_updated_combined_df_clean$pID), 
                                             unique(N1_titer_updated_combined_df_clean$pID))) #total 320 pIDs have NAI titer measured for either N1 or N2 or both

#dataframe for these 314 unique pIDs, measured N1, measured N2 (2 indivs for which N2 was measured but no N1 probably)

num_pcrs_with_NAI_meas <- length(intersect(pIDs_with_NAI_titer_measurements, fluA_positive$pID)) #27 total cases
pcrs_with_NAI_meas <- intersect(pIDs_with_NAI_titer_measurements, fluA_positive$pID)


N1_fold_change <- N1_titer_updated_combined_df_clean %>% 
  group_by(pID) %>% 
  group_modify(~calculate_fold_change(.x))

N2_fold_change <- N2_titer_updated_combined_df_clean %>% 
  group_by(pID) %>% 
  group_modify(~calculate_fold_change(.x))

# ======================================================
#     NAI titers for PCR confirmed infections
# ======================================================

NAI_N1_infections <- NAI_pcr_confirmed(fold_change_df = N1_fold_change, N1_or_N2 = "N1")
NAI_N2_infections <- NAI_pcr_confirmed(fold_change_df = N2_fold_change, N1_or_N2 = "N2")

#for some reasons, the following column were "characters", not "numeric", but this doesn't cause any problem in analysis
NAI_N2_infections$`x2 = Titer @ NA < 50%` <- as.numeric(NAI_N2_infections$`x2 = Titer @ NA < 50%`)
NAI_pcr_confirmed_infs <- rbind(NAI_N1_infections[,-c(10:13)], NAI_N2_infections)

NAI_pcr_confirmed_upd <- NAI_pcr_confirmed_infs %>% 
  filter(!pID %in% subtyped_unsubtyped_pcr$pID[subtyped_unsubtyped_pcr$subtype == ""]) %>% 
  distinct(pID, antigen, .keep_all = T)

NAI_pcr_confirmed_outside_infection_date <- rbind(NAI_N1_infections[,-c(10:13)] %>% 
                                                    group_by(pID) %>%
                                                    slice_tail(n = -1)%>% 
                                                    filter(!timepoint_date %in% after_infection_date) %>% 
                                                    ungroup(),
                                                  NAI_N2_infections %>% 
                                                    group_by(pID) %>%
                                                    slice_tail(n = -1) %>% 
                                                    filter(!timepoint_date %in% after_infection_date) %>% 
                                                    ungroup())

# ==========================================================
# Comparing NAI titers with or w/o PCR confirmed infections
# ==========================================================

t_test_between_fc_during_inf_or_no_inf <- t.test(log(NAI_pcr_confirmed_outside_infection_date$fc_from_prev_ll2, base = 2),
                                                 log(NAI_pcr_confirmed_upd$fold_change_during_infection, base = 2), var.equal = F)

p_value_fc <- t_test_between_fc_during_inf_or_no_inf$p.value
print(paste("p-value:", round(p_value_fc, 4)))

df_fc <- data.frame(
  value = c(log(NAI_pcr_confirmed_outside_infection_date$fc_from_prev_ll2, base = 2), 
            log(NAI_pcr_confirmed_upd$fold_change_during_infection, base = 2)),
  group = rep(c("outside of detected infection \nby PCR", "during infection detected \nby PCR"), 
              c(length(NAI_pcr_confirmed_outside_infection_date$fc_from_prev_ll2), 
                length(NAI_pcr_confirmed_upd$fold_change_during_infection)))
)

df_fc_summary <- df_fc %>%
  group_by(group) %>%
  summarise(mean_val = mean(value),
            sd_val = sd(value),
            se_val = sd(value)/sqrt(n()),
            n = n(),
            .groups = "drop")

df_fc <- df_fc %>% mutate(group = case_when(
  group == "during infection detected \nby PCR" ~ "Yes",
  group == "outside of detected infection \nby PCR" ~ "No"))

# ======================================================
#     Fitting the antibody ceiling line to NAI titer
# ======================================================

fit_line_fb_titer <- lm(log(fold_change_during_infection, base = 2) ~ log(titer_before_infection, base = 2), 
                        data = NAI_pcr_confirmed_upd)

coefs <- coef(fit_line_fb_titer)
a <- coefs[1]   # intercept
b <- coefs[2]   # slope

# ======================================================
#      NAI analysis sensitivity and specificity
# ======================================================

threshold_intercept_series <- c(3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5,7, 7.5, 8, 8.5)
ceiling_line_slope <- b 

sensitivity_res_N1 <- map(threshold_intercept_series, 
                          ~ NAI_sensitivity(fold_change_df = N1_fold_change, 
                                            NAI_pcr_confirmed = NAI_N1_infections,
                                            N1_or_N2 = "N1", 
                                            threshold_boost_intercept = .x))

sensitivity_res_N2 <- map(threshold_intercept_series, 
                          ~ NAI_sensitivity(fold_change_df = N2_fold_change, 
                                            NAI_pcr_confirmed = NAI_N2_infections, 
                                            N1_or_N2 = "N2", 
                                            threshold_boost_intercept = .x))

specificity_res_N1 <- map(threshold_intercept_series,
                          ~ NAI_specificity(titer_df = N1_titer_updated_combined_df_clean, 
                                            low_circ_df = H1_low_circ, 
                                            threshold_boost_intercept = .x,
                                            inferred_inf_pIDs = sensitivity_res_N1[[.x]][[3]]))

specificity_res_N2 <- map(threshold_intercept_series,
                          ~ NAI_specificity(titer_df = N2_titer_updated_combined_df_clean, 
                                            low_circ_df = H3_low_circ, 
                                            threshold_boost_intercept = .x,
                                            inferred_inf_pIDs = sensitivity_res_N2[[.x]][[3]]))

df_sensitivity_N1 <- lapply(sensitivity_res_N1, `[[`, 1) 
df_sensitivity_N1 <- do.call(rbind, df_sensitivity_N1)

df_sensitivity_N2 <- lapply(sensitivity_res_N2, `[[`, 1) 
df_sensitivity_N2 <- do.call(rbind, df_sensitivity_N2)

df_sensitivity <- rbind(df_sensitivity_N1, df_sensitivity_N2)

df_specificity_N1 <- do.call(rbind, specificity_res_N1) %>% mutate(antigen = "N1")
df_specificity_N2 <- do.call(rbind, specificity_res_N2) %>% mutate(antigen = "N2")

df_specificity <- rbind(df_specificity_N1, df_specificity_N2)

data_final <- left_join(df_sensitivity, df_specificity)

data_final <- data_final %>% 
  group_by(antigen, threshold_intercept) %>% 
  mutate(sensitivity = num_true_positive/(num_true_positive + num_false_negative),
         specificity = num_true_negative/(num_true_negative + num_false_positive)) %>% ungroup()

chosen_threshold_N1 <- data_final %>% 
  filter(antigen == "N1" & round(sensitivity,2) >= 0.9 & round(specificity,2) >= 0.9) %>% 
  dplyr::select(threshold_intercept)

chosen_threshold_N2 <- data_final %>% 
  filter(antigen == "N2" & round(sensitivity,2) >= 0.9 & round(specificity,2) >= 0.9) %>% 
  dplyr::select(threshold_intercept)

#this has to be hand-coded as there is a slight increase in specificity
#for 3.5 than 3, keeping sensitivity the same
chosen_threshold_N2 <- chosen_threshold_N2[2,1] 

# ======================================================
#      Possible infections from the NAI analysis
# ======================================================

N1_index <- which(df_sensitivity_N1$threshold_intercept == chosen_threshold_N1$threshold_intercept)
N2_index <- which(df_sensitivity_N2$threshold_intercept == chosen_threshold_N2$threshold_intercept)

#3rd element in the return list contains pIDs and infection dates of possible infected
N1_possible_infections <- sensitivity_res_N1[[N1_index]][[3]] 
N2_possible_infections <- sensitivity_res_N2[[N2_index]][[3]]

#2nd element in the return list contains titer of infected before and after infection date and the corresponding dates
N1_time_of_infection <- sensitivity_res_N1[[N1_index]][[2]]
N2_time_of_infection <- sensitivity_res_N2[[N2_index]][[2]]

#this is done for plotting only the section of NAI rise for each possible infected
N2_time_of_infection_sub <- N2_time_of_infection %>% 
  group_by(pID) %>% 
  slice(which(timepoint_date %in% as.Date(before_boost_date) | timepoint_date %in% as.Date(after_boost_date))) %>% 
  ungroup()

N1_time_of_infection_sub <- N1_time_of_infection %>% 
  group_by(pID) %>% 
  slice(which(timepoint_date %in% as.Date(before_boost_date) | timepoint_date %in% as.Date(after_boost_date))) %>% 
  ungroup()

# ======================================================
#                    Making figures 
# ======================================================

#figure specifics
axis_label <- default_figure_font_size
axis_text <- small_font_size
legend_text <- small_font_size
legend_label <- default_figure_font_size
anno_text <- 3

point_size <- 1
line_size <- 1

fc_comp <- ggplot()  +
  geom_boxplot(data = df_fc, aes(x = group, y = value), outliers = F) +
  geom_jitter(data = df_fc, aes(x = group, y = value),
              width = 0.1, alpha = 0.8, size = point_size) +
  labs(x = "PCR-confirmed infection", 
       y = "NAI fold change (log<sub>2</sub>)") +
  annotate("text", x = 1.5, y = max(df_fc$value)*0.95, 
           label = paste("p =", signif(p_value_fc, digits = 2)), size = anno_text) +
  theme_bw() +
  guides(color = "none") +
  theme(
    axis.text.x = element_text(size = axis_text),
    axis.text.y = element_text(size = axis_text),
    axis.title.y = ggtext::element_markdown(size = axis_label),
    axis.title.x = element_text(size = axis_label),
    panel.grid = element_blank())

ct <- cor.test(log(NAI_pcr_confirmed_upd$titer_before_infection, base = 2), 
               log(NAI_pcr_confirmed_upd$fold_change_during_infection, base = 2), method = "pearson")

r  <- unname(ct$estimate)
p  <- ct$p.value

# format p-value nicely
p_lab <- if (p < 0.001) "p < 0.001" else paste0("p = ", signif(p, 2))

lab_corr <- paste0("r = ", round(r, 2), "\n", p_lab)

# illustrating ceilling effect
antibody_ceiling <- ggplot() +
  geom_point(data = NAI_pcr_confirmed_upd,
             aes(log(titer_before_infection, base = 2), log(fold_change_during_infection, base = 2), 
                 shape = factor(antigen)),
             #color = pID), 
             size = point_size, stroke = 2) +
  scale_shape_manual(values = c(1,5)) +
  theme_bw() +
  geom_smooth(data = NAI_pcr_confirmed_upd,
              aes(log(titer_before_infection, base = 2), log(fold_change_during_infection, base = 2)),
              method = "lm", se = T, colour = "black", linewidth = 1) + 
  annotate("text",
           x = 5.2, y = 5.6,
           label = lab_corr,
           hjust = 0, vjust = 0,
           size = anno_text) +
  labs(x = "NAI titer before infection (log<sub>2</sub>)",
       y = "", shape = "antigen") +
  theme(
    #text = element_text(family = "Arial"),
    axis.title.x = ggtext::element_markdown(size = axis_label),  # x-axis label
    axis.title.y = element_text(size = axis_label),  # y-axis label
    axis.text.x  = element_text(size = axis_text),  # x tick labels
    axis.text.y  = element_text(size = axis_text),   # y tick labels
    legend.text = element_text(size = legend_text),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank()
  ) 

#fig 1 - NAI analysis
fig1 <- fc_comp + plot_spacer() + antibody_ceiling +
  plot_layout(widths = c(2, 0.03, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 10),
    plot.tag.position = c(0, 1)
  ) &
  coord_cartesian(ylim = c(-3, 8)) #common y-axis limits for two panels

#visualizing sensitivity and specificity of ELLA Assay
N1_Spec_sens <- ggplot() +
  geom_point(data = data_final %>% filter(antigen == "N1"), aes(threshold_intercept, sensitivity), color = "orangered3", size = 3, shape = 1) +
  geom_point(data = data_final %>% filter(antigen == "N1"), aes(threshold_intercept, specificity), color = "slateblue3", size = 3, shape = 1) +
  geom_line(data = data_final %>% filter(antigen == "N1"), aes(threshold_intercept, sensitivity), color = "orangered3") +
  geom_line(data = data_final %>% filter(antigen == "N1"), aes(threshold_intercept, specificity), color = "slateblue3") +
  theme_bw() +
  geom_vline(xintercept = chosen_threshold_N1$threshold_intercept, linetype = "dashed", color = "grey70") +
  ylab("Sensitivity and specificity") + #ELLA assay sensitivity and specificity
  xlab("Fold change threshold at baseline log<sub>2</sub> NAI titer = 0") + #Threshold intercept for the antibody ceiling line
  ggtitle("N1 antigen") +
  annotate("text",
           x = 7, y = 0.2,          # coordinates in data space
           label = "Sensitivity",
           hjust = 0, vjust = 1,
           size = anno_text, color = "orangered3") +
  annotate("text",
           x = 7, y = 0.95,          # coordinates in data space
           label = "Specificity",
           hjust = 0, vjust = 1,
           size = anno_text, color = "slateblue3") +
  theme(axis.title.x = ggtext::element_markdown(size = axis_label),  # x-axis label
        axis.title.y = element_text(size = axis_label),  # y-axis label
        axis.text.x  = element_text(size = axis_text),  # x tick labels
        axis.text.y  = element_text(size = axis_text),
        plot.title = element_text(size = axis_label),
        panel.grid = element_blank())   # y tick labels

N2_Spec_sens <- ggplot() +
  geom_point(data = data_final %>% filter(antigen == "N2"), aes(threshold_intercept, sensitivity), color = "orangered3", size = 3, shape = 1) +
  geom_point(data = data_final %>% filter(antigen == "N2"), aes(threshold_intercept, specificity), color = "slateblue3", size = 3, shape = 1) +
  geom_line(data = data_final %>% filter(antigen == "N2"), aes(threshold_intercept, sensitivity), color = "orangered3") +
  geom_line(data = data_final %>% filter(antigen == "N2"), aes(threshold_intercept, specificity), color = "slateblue3") +
  theme_bw() +
  geom_vline(xintercept = chosen_threshold_N2$threshold_intercept, linetype = "dashed", color = "grey70") +
  ylab("Sensitivity and specificity") +
  xlab("Fold change threshold at baseline log<sub>2</sub> NAI titer = 0") + #Threshold intercept for the antibody ceiling line
  ggtitle("N2 antigen") +
  annotate("text",
           x = 7, y = 0.2,          # coordinates in data space
           label = "Sensitivity",
           hjust = 0, vjust = 1,
           size = anno_text, color = "orangered3") +
  annotate("text",
           x = 7, y = 0.95,          # coordinates in data space
           label = "Specificity",
           hjust = 0, vjust = 1,
           size = anno_text, color = "slateblue3") +
  theme(axis.title.x = ggtext::element_markdown(size = axis_label),  # x-axis label
        axis.title.y = element_text(size = axis_label),  # y-axis label
        axis.text.x  = element_text(size = axis_text),  # x tick labels
        axis.text.y  = element_text(size = axis_text),
        plot.title = element_text(size = axis_label),
        panel.grid = element_blank())

#fig 2 - NAI analysis
fig2 <- (N1_Spec_sens | N2_Spec_sens) + 
  plot_layout(axis_titles = "collect")

#scale factor to plot NAI titer dynamics on top of surveilance data
scale_factor <- 2 * max(HK_surveilance_upd_long$num_pos_cases_1000visitis, na.rm = T) / 
  max(log(N2_titer_updated_combined_df_clean$estimated_titer, base = 2), na.rm = T)


fig3 <- make_inferred_infections_fig(
  low_circ = H1_low_circ,
  time_of_infection_sub = N1_time_of_infection_sub,
  color_key = "H1N1",
  surv_subtype = "H1_pos_1000visits",
  y_axis_name = "H1 positive cases per 1000 visits"
)

fig4 <- make_inferred_infections_fig(
  low_circ = H3_low_circ,
  time_of_infection_sub = N2_time_of_infection_sub,
  color_key = "H3N2",
  surv_subtype = "H3_pos_1000visits",
  y_axis_name = "H3 positive cases per 1000 visits"
)

# Note: Figure numbering here is specific to this script, does not match numbering in the manuscript
results_dir <- "results/NAI_inferred_infections/"
dir.create(results_dir, recursive = T, showWarnings = F)

ggsave(
  filename = "NAI_fold_change_PCR_confirmed.pdf",
  plot = fig1,
  path = results_dir,
  width = supp_fig_witdh_full,
  height = 3.25/1.4,
  units = "in"
)

ggsave(
  filename = "NAI_sensitivity_specificity.pdf",
  plot = fig2,
  path = results_dir,
  width = supp_fig_witdh_full,
  height = 3.25/1.4,
  units = "in"
)

ggsave(
  filename = "NAI_inferred_infections_H1.pdf",
  plot = fig3,
  path = results_dir,
  width = supp_fig_witdh_full,
  height = 5.5,
  units = "in"
)

ggsave(
  filename = "NAI_inferred_infections_H3.pdf",
  plot = fig4,
  path = results_dir,
  width = supp_fig_witdh_full,
  height = 5.5,
  units = "in"
)

# Export list of inferred infections 
write.csv(N1_possible_infections, "results/NAI_inferred_infections/Inferred_N1_infections.csv", row.names = F)
write.csv(N2_possible_infections, "results/NAI_inferred_infections/Inferred_N2_infections.csv", row.names = F)
