library(tidyverse)
library(cowplot)
library(ggh4x)
theme_set(theme_cowplot())
library(GGally)
library(imprinting)
library(reporttools) 
library(ggpubr) 


# ===== Load data ======

source("R/base_functions.R")
source("R/load_data.R")

dir.create("results/endpoint_analyses", recursive = T, showWarnings = F)

titers <- titers %>%
  # Keep a single measurement per strain/timepoint/titer type
  # prioritizing batches that measured the same person in multiple years
  select_titer_batches() %>%
  filter(year >= 2) %>%
  select(-matches('FRNT_assay')) %>%
  arrange(pID, strain, timepoint)

# If using this code later for drive 2, review.
stopifnot(all(unique(titers$drive) == 1))

# =====  Distribution of time since vaccination on day 30  =====
#                     (exploratory/interactive)

titers %>%
  filter(timepoint == 30) %>%
  select(pID, year, ndays_since_year_vax) %>%
  unique() %>%
  ggplot(aes(x = ndays_since_year_vax)) +
    geom_histogram(binwidth = 1) +
    facet_wrap("year") +
    geom_vline(xintercept = c(28, 28 + 7, 28 + 14))

titers %>%
  filter(timepoint == 30) %>%
  select(pID, year, ndays_since_year_vax) %>%
  unique() %>%
  ggplot(aes(x = ndays_since_year_vax, color = factor(year))) +
    stat_ecdf() +
    geom_vline(xintercept = c(28, 28 + 7, 28 + 14), linetype = "dashed", alpha = 0.5) +
    xlab("Days since year vaccination") +
    ylab("Cumulative proportion") +
    scale_color_discrete(name = "Year")

titers %>%
  filter(timepoint == 30) %>%
  select(pID, year, ndays_since_year_vax) %>%
  unique() %>%
  group_by(year) %>%
  summarise(
    `< 28`    = mean(ndays_since_year_vax < 28),
    `[28, 42]` = mean(ndays_since_year_vax >= 28 & ndays_since_year_vax <= 42),
    `> 42`    = mean(ndays_since_year_vax > 42)
  ) %>%
  print()


# ====== Visualizing influenza circulation =======
vaccination_times <- timepoint_days %>%
  # filter to DRIVE I
  filter(str_detect(pID, "D1"),
         !is.na(year_vdate)) %>%
  select(pID, year, year_vdate) %>%
  unique()

midpoint_vaccination_dates <- vaccination_times %>%
  group_by(year) %>%
  summarise(midpoint_vdate = min(year_vdate) + (max(year_vdate) - min(year_vdate)) / 2)

# For plotting only: each round spans from the first day of the ISO week of its
# earliest vaccination to the first day of the ISO week of its latest vaccination.
round_vaccination_spans <- vaccination_times %>%
  group_by(year) %>%
  summarise(xmin = snap_to_iso_week_start(min(year_vdate)),
            xmax = snap_to_iso_week_start(max(year_vdate)))

drive_1_flu_infection_dates <- positive_pcr_tests %>%
  filter(str_detect(pID, "D1"), result != "cov2_positive") %>%
  arrange(cdate)

# Check that the minimum cdate in drive_1_flu_infection_dates is 
# after the maximum year_vdate in year 3
min(drive_1_flu_infection_dates$cdate) > max(vaccination_times$year_vdate[vaccination_times$year == 3])

# What's more important: Check that all participants gave their year 3 day 30 sample before they were infected
titers %>%
  filter(timepoint == 30, year == 3) %>%
  filter(infection_before_sample_anyflu_PCR | infection_before_sample_H3N2_NAI | infection_before_sample_H1N1_NAI) %>%
  nrow() == 0

last_y3d30_date <- titers %>%
  filter(timepoint == 30, year == 3) %>%
  filter(timepoint_date == max(timepoint_date)) %>%
  pull(timepoint_date)

# Year labels centered within each calendar year, to be used in place of
# scale_x_date's default labels (which sit directly under the Jan 1 tick).
# Centered on the full calendar year (not the span of available data), so a
# partial final year is labeled at the same relative position as full years.
HK_year_label_positions <- HK_flu_surveillance %>%
    filter(Year >= 2019) %>%
    distinct(Year) %>%
    mutate(x = ymd(paste0(Year, "-01-01")) + (ymd(paste0(Year + 1, "-01-01")) - ymd(paste0(Year, "-01-01"))) / 2)

HK_flu_timeseries <- HK_flu_surveillance %>%
    select(Year, Week, From, To, matches("proportion")) %>%
    select(-AandB_proportion) %>%
    pivot_longer(cols = matches("proportion"),
                 names_to = "subtype", values_to = "proportion") %>%
    filter(Year >= 2019) %>%
    mutate(subtype = case_when(
        subtype == "H1_proportion" ~ "H1N1",
        subtype == "H3_proportion" ~ "H3N2",
        subtype == "B_proportion" ~ "B"
    )) %>%
    mutate(subtype = factor(subtype, levels = c("H3N2", "H1N1", "B"))) %>%
    ggplot(aes(x = From, y = proportion, color = subtype)) +
    geom_rect(
        data = round_vaccination_spans,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "gray90"
    ) +
    geom_line(aes(group = subtype), linewidth = 0.7) +
    scale_x_date(date_breaks = "1 year", labels = NULL) +
    xlab(NULL) +
    ylab("Respiratory specimens positive\nfor influenza in Hong Kong") +
    scale_color_manual(name = "", values = subtype_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.3, 0.1)) +
    baseline_figure_settings +
    theme(
        legend.position = "top",
        legend.justification = "left",
        legend.box = "horizontal",
        legend.box.spacing = margin(0.3)
    ) +
    geom_text(
        data = midpoint_vaccination_dates %>% mutate(y = if_else(year %in% c(3, 5), 0.265, 0.29)),
        aes(x = midpoint_vdate, y = y, label = paste("Round", year)),
        color = "gray50",
        size = default_figure_font_size,
        size.unit = "pt",
        vjust = 0
    ) +
    geom_point(
        # Snap to ISO-week start for plotting only
        data = drive_1_flu_infection_dates %>% mutate(cdate = snap_to_iso_week_start(cdate)),
        aes(x = cdate, y = 0.323),
        inherit.aes = FALSE,
        shape = 16,
        alpha = 0.5,
        size = 1.5
    ) +
    geom_text(
        data = drive_1_flu_infection_dates %>%
            summarise(x = min(cdate) + (max(cdate) - min(cdate)) / 2),
        aes(x = x, y = 0.335, label = "PCR-confirmed infections"),
        inherit.aes = FALSE,
        hjust = 0.5, vjust = 0,
        size = default_figure_font_size,
        size.unit = "pt",
        lineheight = 0.8
    ) +
    geom_text(
        data = HK_year_label_positions,
        aes(x = x, y = -Inf, label = Year),
        inherit.aes = FALSE,
        color = "black",
        size = default_figure_font_size,
        size.unit = "pt",
        vjust = 1.6
    ) +
    coord_cartesian(ylim = c(NA, 0.3), clip = "off")

save_plot("results/HK_flu_surveillance_timeseries.pdf",
          HK_flu_timeseries,
          base_height = 1.9, base_width = 4.25)
save_plot("results/HK_flu_surveillance_timeseries.svg",
          HK_flu_timeseries,
          base_height = 1.9, base_width = 4.25)          


#====== Analysis of FRNT/HAI titers ======

# Table of GMTs
GMT_table <- titers %>%
  flag_vaccine_strains() %>%
  filter(is_vaccine_strain) %>%
  annotate_years_with_vaccine_updates() %>% # Irrelevant here but required for next function
  annotate_with_n_prior_vax() %>%
  mutate(n_previous_vax = case_when(
    str_detect(treatment, "V") ~ as.character(n_previous_vax),
    !str_detect(treatment, "V") ~ "placebo"
  )) %>%
  mutate(n_previous_vax = factor(n_previous_vax, levels = c("placebo", as.character(0:3)))) %>%
  group_by(timepoint, n_previous_vax, year, subtype) %>%
  summarise(GMT = round(2^mean(log2_titer), 2),
            IQR_lower = round(quantile(2^log2_titer, 0.25), 2),
            IQR_upper = round(quantile(2^log2_titer, 0.75), 2)) %>%
  ungroup() %>%
  mutate(value = paste0(GMT, " (", IQR_lower, "-", IQR_upper, ")")) %>%
  mutate(year_subtype = paste(year, subtype, sep = "-")) %>%
  select(-year, -subtype, -GMT, -matches("IQR")) %>%
  pivot_wider(names_from = "year_subtype", values_from = "value")

write_csv(GMT_table, "results/endpoint_analyses/GMT_table.csv")

# Table of pairwise comparisons
full_comparison_table <- run_pairwise_wilcoxon(
  titers, response_var = "log2_titer", p.adjust.method = "holm") %>%
  flag_vaccine_strains() %>%
  filter(is_vaccine_strain) %>%
  annotate_and_order_strains_and_subtypes(strain_levels = strain_levels)

write_csv(full_comparison_table, "results/endpoint_analyses/full_comparison_table.csv")

# In year 4, ratio of VVVV against all other non-placebo groups
VVVV_comparisons_Y4 <- full_comparison_table %>%
  filter(year == 4, timepoint == 30,
         (treatment_1 == "VVVV" | treatment_2 == "VVVV")) %>%
  filter(treatment_2 != "PPPP") %>%
  mutate(GMT_ratio = GMT_treatment_1 / GMT_treatment_2) %>%
  mutate(percent_lower = (1 - GMT_ratio) * 100)

print(VVVV_comparisons_Y4)

VVVV_comparisons_Y4 %>%
  summarise(
    mean = mean(percent_lower),
    min = min(percent_lower),
    max = max(percent_lower)
  )

write_csv(VVVV_comparisons_Y4, "results/endpoint_analyses/VVVV_comparisons_Y4.csv")

# Plots of comparisons between groups
# First, FRNT/HAI titers
titers_plot_y2 <- plot_postvax_response(response_data = titers,
                                        response_var = 'log2_titer',
                                        year = 2) +
  ylab("Antibody titer") +
  theme(legend.position = "none")

save_plot("results/endpoint_analyses/Y2_postvax_titers.pdf",
          titers_plot_y2,
          base_height = supp_fig_height_full / 1.7,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y2_postvax_titers.png",
          titers_plot_y2,
          base_height = supp_fig_height_full / 1.7,
          base_width = supp_fig_witdh_full, dpi = 300)
                                      
titers_plot_y3 <- plot_postvax_response(response_data = titers,
                                        response_var = 'log2_titer',
                                        year = 3)  +
  ylab("Antibody titer") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("results/endpoint_analyses/Y3_postvax_titers.pdf",
          titers_plot_y3,
          base_height = supp_fig_height_full,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y3_postvax_titers.png",
          titers_plot_y3,
          base_height = supp_fig_height_full,
          base_width = supp_fig_witdh_full, dpi = 300)

titers_plot_y4 <- plot_postvax_response(response_data = titers,
                                        response_var = 'log2_titer',
                                        year = 4)  +
  ylab("Antibody titer") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("results/endpoint_analyses/Y4_postvax_titers.pdf",
          titers_plot_y4,
          base_height = supp_fig_height_full,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y4_postvax_titers.png",
          titers_plot_y4,
          base_height = supp_fig_height_full,
          base_width = supp_fig_witdh_full, dpi = 300)

titers_plot_y4_days0_30_182 <- plot_postvax_response(
  response_data = titers %>% filter(timepoint %in% c(0, 30, 182)),
  response_var = "log2_titer",
  year = 4)

# Main text figure combining years 3 and 4
years_3_4_postvax <- plot_years_3_4_postvax(titers, response_var = "log2_titer") +
  ylab("Antibody titers")

save_plot("results/endpoint_analyses/years_3_4_postvax_titers.pdf",
          years_3_4_postvax,
          base_height = 4.3, base_width = main_text_width_full)

save_plot("results/endpoint_analyses/years_3_4_postvax_titers.svg",
          years_3_4_postvax,
          base_height = 4.3, base_width = main_text_width_full)


# Fraction of participants with >=4 fold rises
y2_fraction_fourfold_or_greater <- plot_fraction_by_group(
  titers, response_var = "fourfold_rise_or_greater", year = 2)

y3_fraction_fourfold_or_greater <- plot_fraction_by_group(
  titers, response_var = "fourfold_rise_or_greater", year = 3)

y4_fraction_fourfold_or_greater <- plot_fraction_by_group(
  titers, response_var = "fourfold_rise_or_greater", year = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("results/endpoint_analyses/Y2_fraction_fourfold_or_greater.pdf",
          y2_fraction_fourfold_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y2_fraction_fourfold_or_greater.png",
          y2_fraction_fourfold_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)

save_plot("results/endpoint_analyses/Y3_fraction_fourfold_or_greater.pdf",
          y3_fraction_fourfold_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y3_fraction_fourfold_or_greater.png",
          y3_fraction_fourfold_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)

save_plot("results/endpoint_analyses/Y4_fraction_fourfold_or_greater.pdf",
          y4_fraction_fourfold_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y4_fraction_fourfold_or_greater.png",
          y4_fraction_fourfold_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)

# Fraction of participants titers at or above 40
y2_fraction_titer_40_or_greater <- plot_fraction_by_group(
  titers, response_var = "titer_40_or_greater", year = 2)

y3_fraction_titer_40_or_greater <- plot_fraction_by_group(
  titers, response_var = "titer_40_or_greater", year = 3)

y4_fraction_titer_40_or_greater <- plot_fraction_by_group(
  titers, response_var = "titer_40_or_greater", year = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot("results/endpoint_analyses/Y2_fraction_titer_40_or_greater.pdf",
          y2_fraction_titer_40_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y2_fraction_titer_40_or_greater.png",
          y2_fraction_titer_40_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)

save_plot("results/endpoint_analyses/Y3_fraction_titer_40_or_greater.pdf",
          y3_fraction_titer_40_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y3_fraction_titer_40_or_greater.png",
          y3_fraction_titer_40_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)

save_plot("results/endpoint_analyses/Y4_fraction_titer_40_or_greater.pdf",
          y4_fraction_titer_40_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y4_fraction_titer_40_or_greater.png",
          y4_fraction_titer_40_or_greater,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)


# Plotting FRNT titers to H1N1 in year 4 separately
Y4_postvax_H1N1_FRNT <- plot_postvax_response(response_data = H1N1_FRNT_titers, response_var = "log2_titer",
                      vaccine_strains_only = F, year = 4) +
  ylab('FRNT titer')

save_plot("results/endpoint_analyses/Y4_postvax_H1N1_FRNT.pdf", Y4_postvax_H1N1_FRNT, base_height = 6, base_width = 6)
save_plot("results/endpoint_analyses/Y4_postvax_H1N1_FRNT.png", Y4_postvax_H1N1_FRNT, base_height = 6, base_width = 6, dpi = 300)


# Same for HAI titers to H3N2 in year 4
Y4_postvax_H3N2_HAI <- plot_postvax_response(response_data = H3N2_HAI_titers, response_var = "log2_titer",
                      vaccine_strains_only = F, year = 4) +
  ylab('HAI titer')

save_plot("results/endpoint_analyses/Y4_postvax_H3N2_HAI.pdf", Y4_postvax_H3N2_HAI, base_height = 8, base_width = 5)
save_plot("results/endpoint_analyses/Y4_postvax_H3N2_HAI.png", Y4_postvax_H3N2_HAI, base_height = 8, base_width = 5, dpi = 300)

Y4_fraction_titer_40_or_greater_H3N2_HAI <- plot_fraction_by_group(H3N2_HAI_titers, response_var = "titer_40_or_greater", year = 4) +
  ylab(bquote("Fraction of participants with HAI titers" ~ symbol("\263") * "40"))

save_plot(
  "results/endpoint_analyses/Y4_fraction_titer_40_or_greater_H3N2_HAI.pdf",
  Y4_fraction_titer_40_or_greater_H3N2_HAI,
  base_height = supp_fig_height_full / 3,
  base_width = supp_fig_witdh_full
)

ggsave(
  filename = "results/endpoint_analyses/Y4_fraction_titer_40_or_greater_H3N2_HAI.png",
  plot = Y4_fraction_titer_40_or_greater_H3N2_HAI,
  height = supp_fig_height_full / 3,
  width = supp_fig_witdh_full,
  units = "in",
  dpi = 300
)

#========== Analysis of Luminex titers =========
luminex_vax_plot <- plot_postvax_response(luminex_vax_strains,
                                          response_var = 'log10_luminex',
                                          year = 3)

save_plot("results/endpoint_analyses/Y3_postvax_luminex_vax_strains.pdf",
          luminex_vax_plot,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full)
save_plot("results/endpoint_analyses/Y3_postvax_luminex_vax_strains.png",
          luminex_vax_plot,
          base_height = supp_fig_height_full - 2,
          base_width = supp_fig_witdh_full, dpi = 300)
                
# total luminex vs titers
luminex_vs_titers_plot <- plot_luminex_vs_titer_correlations(titers, luminex_vax_strains)
save_plot("results/endpoint_analyses/luminex_vs_titers.pdf", luminex_vs_titers_plot, base_height = 8, base_width = 8)
save_plot("results/endpoint_analyses/luminex_vs_titers.png", luminex_vs_titers_plot, base_height = 8, base_width = 8, dpi = 300)

# Exporting de-identified data
timepoints_and_prior_infections <- titers %>%
  select(pID, year, treatment, timepoint, ndays_since_year_vax, matches('infection_before_sample')) %>%
  select(-matches("recent"), -matches("subtype_matched")) %>%
  unique() %>%
  mutate(binned_days_since_intervention = case_when(
    timepoint == 30 & ndays_since_year_vax < 28  ~ "<28 days",
    timepoint == 30 & ndays_since_year_vax >= 28 & ndays_since_year_vax <= 42 ~ "28-42 days",
    timepoint == 30 & ndays_since_year_vax > 42  ~ ">42 days"
  )) %>%
  select(pID, year, treatment, timepoint, binned_days_since_intervention, 
         matches("PCR"), matches("NAI"))

write_csv(timepoints_and_prior_infections, file = "results/data_sharing/timepoints_and_prior_infections.csv")

