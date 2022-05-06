# Queries for BSAI Greenland turbot Relative Population Numbers (RPNs) and
# RPN-weighted length frequencies from the AFSC longline survey.
# Contacts jane.sullivan@noaa.gov, kevin.siwicke@noaa.gov
# Updated May 2022
# Generated using  R version 4.1.2

# set up ----
libs <- c("tidyverse", "RODBC", "zoo")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

theme_set(theme_minimal(base_size = 10))

username_akfin <- '' # fill in your akfin creds
password_akfin <- ''

channel_akfin <- odbcConnect("akfin", uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

# RPN queries ----

# BSAI RPNs (BS = odd years, AI = even years). Full time series starts in 1997.
# The raw 1996 AI in area 15 (NE AI Slope) was lost, though we have area estimates. 
rpn <- sqlQuery(channel_akfin, query = paste0("
                select    year, council_management_area as fmp, rpn, rpn_var
                from      afsc.lls_fmp_subarea_all_strata
                where     country = 'United States' and
                          species_code = '10115' and
                          year >= 1996 and 
                          ((council_management_area in 'Bering Sea' and mod(year, 2) <> 0) or
                          (council_management_area in ('Aleutians') and mod(year, 2) = 0))
                order by  year asc
                ")) %>% 
  rename_all(tolower) 

write_csv(rpn, 'data/bsai_turbot_rpns.csv')

# Geographic area RPNs Note that I included 1996 here. The LLS team doesn't have a
# strong recommendation on whether or not these data should be used.
arearpns <- sqlQuery(channel_akfin, query = ("
                select    *
                from      afsc.lls_area_rpn_all_strata
                where     species_code = '10115' and country = 'United States' and 
                          council_sablefish_management_area in ('Aleutians', 'Bering Sea') and 
                          exploitable = 1
                order by  year asc
                ")) %>% 
  rename_all(tolower) 

write_csv(arearpns, 'data/bsai_turbot_arearpns.csv')

arearpns <- arearpns %>% 
  group_by(year, geographic_area_name) %>%
  mutate(se = sqrt(sum(rpn_var, na.rm = TRUE)))

# geographic area RPNs ----

ggplot(arearpns, aes(x = year, y = rpn / 1e3)) +
  geom_ribbon(aes(x = year,
                  ymin = (rpn - 1.96 * se) / 1e3,
                  ymax = (rpn + 1.96 * se) / 1e3), 
              col = NA,
              alpha = 0.25) + 
  geom_point() +
  geom_line() +
  facet_wrap(~geographic_area_name) +
  theme_minimal(base_size = 13) +
  labs(x = 'Year', y = 'RPN')

ggsave("results/lls_bsai_turbot_arearpn.png", height = 8, width = 7, unit = 'in')

# interp functions ----

# different methods for calculating off-year surveys

# use long-term mean to get average proportion by area to approximate BS or AI in missing survey years
f_ltmean <- function(df, # cols = fmp ('Aleutians', 'Bering Sea'), year, rpn
                     method_lab = "ltmean_olddat", # label for method
                     mean_syr = 1996, # starting year for long-term mean calcs
                     index_syr = 1996 # starting year for index
                     ) {
  
  df <- df %>% 
    right_join(expand_grid(fmp = unique(df$fmp),
                           year = min(df$year):max(df$year)))
  
  df <- df %>% 
    pivot_wider(id_cols = year, values_from = rpn, names_from = fmp) %>% 
    left_join(df %>% 
                group_by(fmp) %>% 
                mutate(mean_rpn = mean(rpn[year >= mean_syr], na.rm = TRUE)) %>% 
                ungroup() %>% 
                mutate(prop = mean_rpn / sum(unique(mean_rpn), na.rm = TRUE)) %>% 
                filter(!is.na(rpn)) %>% 
                select(year, prop)) %>% 
    mutate(bsairpn = ifelse(!is.na(Aleutians), Aleutians / prop, `Bering Sea` / prop),
           Aleutians = ifelse(is.na(Aleutians), bsairpn - `Bering Sea`, Aleutians),
           `Bering Sea` = ifelse(is.na(`Bering Sea`), bsairpn - Aleutians, `Bering Sea`)) %>% 
    pivot_longer(cols = c(Aleutians, `Bering Sea`), names_to = "fmp", values_to = "rpn") %>% 
    select(year, fmp, rpn) %>% 
    mutate(method = method_lab) %>% 
    arrange(fmp, year) %>% 
    filter(year >= index_syr)
  
  return(df)
}

# use linear interpolation to get rpns/cvs in off-survey years; end years set
# equal to nearest year
f_linapprox <- function(df, # cols = fmp ('Aleutians', 'Bering Sea'), year, rpn, rpn_var
                     method_lab = "lininterp_newdat", # label for method
                     index_syr = 1996 # starting year for index
) {
  
  df <- df %>% 
    mutate(rpn_cv = sqrt(rpn_var) / rpn) %>% 
    right_join(expand_grid(fmp = unique(df$fmp),
                           year = unique(df$year))) %>% 
    filter(year >= index_syr) %>% 
    arrange(fmp, year) %>% 
    group_by(fmp) %>% 
    mutate(rpn = zoo::na.approx(rpn, maxgap = 20, rule = 2),
           rpn_cv = zoo::na.approx(rpn_cv, maxgap = 20, rule = 2),
           rpn_var = ifelse(is.na(rpn_var), (rpn_cv * rpn)^2, rpn_var)) %>%
    mutate(method = method_lab) %>% 
    ungroup()
  
  return(df)
}

# old data from akfin -----

#index values taken from 2020 safe. assumed cv=0.2
bs <- c(59328, 63144, 50975, 46616, 23107, 18074, 27850, 16184, 21166, 21001, 20792, 10403)
ai <- c(39262, 37784, 22037, 20170, 15115, 6331, 5374, 3347, 7639, 6315, 3367, 5672, 5697)
yrs <- 1996:2020

rpn_old <- data.frame(fmp = 'Aleutians',
           year = yrs[yrs %% 2 == 0],
           rpn = ai) %>% 
  bind_rows(data.frame(fmp = 'Bering Sea',
                       year = yrs[yrs %% 2 != 0],
                       rpn = bs)) 

write_csv(rpn_old, 'data/bsai_turbot_oldrpns_2020safe.csv')

# compare BS and AI ----

compare_rpns <- f_ltmean(rpn_old, method_lab = 'statusquo_meanratio_olddata', 
                         mean_syr = 1996, index_syr = 1996) %>% 
  bind_rows(f_ltmean(rpn, method_lab = 'statusquo_meanratio_newdat', 
                     mean_syr = 1996, index_syr = 1996)) %>% 
  bind_rows(f_linapprox(rpn, method_lab = 'new_linearapprox_newdat', index_syr = 1996))

ggplot(compare_rpns, aes(x = year, y = rpn / 1e3, 
                         col = method, lty = method, shape = method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~fmp) +
  labs(x = "Year", y = "RPN",
       title ="Greenland turbot Relative Population Numbers",
       subtitle = "Comparing methods for interpolating missing survey years (odd = BS, even = AI)")

ggsave("results/turbot_rpnmethods_fmp.png", units = "in", 
       width = 8, height = 5)  

# compare full BSAI index ----

index <- compare_rpns %>% 
  group_by(year, method) %>% 
  dplyr::summarise(rpn = sum(rpn, na.rm = TRUE),
                   rpn_var = sum(rpn_var, na.rm = TRUE)) %>% 
  # assume fixed cv = 0.2 for old 'longterm mean' method
  mutate(rpn_var = ifelse(is.na(rpn_var) | rpn_var == 0, (rpn * 0.2)^2, rpn_var),
         lci = rpn - 1.96 * sqrt(rpn_var),
         uci = rpn + 1.96 * sqrt(rpn_var)) #%>% View()

ggplot(index, aes(x = year, y = rpn / 1e3,
                  ymin = lci / 1e3, ymax = uci / 1e3, 
                  fill = method, group = method,
                  col = method, lty = method, shape = method)) +
  geom_ribbon(col = 'white', alpha = 0.2) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "RPN",
       title ="Greenland turbot Relative Population Numbers with 95% CI",
       subtitle = "Comparing methods for interpolating missing survey years (odd = BS, even = AI)")

ggsave("results/turbot_rpnmethods.png", units = "in", 
       width = 8, height = 5)  

# Japanese survey ----

# CAUTION! Not a lot of great information available on these data, survey
# methods, etc. Work needs to be done here to make sure consistent stations were
# sampled each year (see table at the end of the
# Database_Background_Instructions on AKFIN Answers for some info on this).

jpn <- sqlQuery(channel_akfin, query = paste0("
                select    country as survey, year, council_management_area as fmp, rpn, rpn_var
                from      afsc.lls_fmp_subarea_all_strata
                where     country = 'Japan' and
                          species_code = '10115' and
                          council_management_area in ('Bering Sea', 'Aleutians') 
                order by  year asc
                ")) %>% 
  rename_all(tolower)

jpn <- jpn %>% 
  bind_rows(f_linapprox(rpn, method_lab = 'linearapprox_newdat') %>% 
              select(-rpn_cv, - method) %>% 
              mutate(survey = "U.S.")) %>% 
  mutate(lci = rpn - 1.96 * sqrt(rpn_var),
         uci = rpn + 1.96 * sqrt(rpn_var)) %>% 
  # loos like 1980 lci drops below 0
  mutate(lci = ifelse(lci < 0, 0, lci))

ggplot(jpn, aes(x = year, y = rpn / 1e3, 
                ymin = lci / 1e3, ymax = uci / 1e3, 
                group = interaction(survey, fmp),
                fill = fmp, col = fmp, 
                lty = fmp, shape = fmp)) +
  geom_ribbon(col = 'white', alpha = 0.2) +
  geom_point() +
  geom_line() +
  labs(x = NULL, y = "RPN",
       fill = NULL, col = NULL, 
       lty = NULL, shape = NULL,
       title ="Greenland turbot Relative Population Numbers with 95% CI",
       subtitle = "Japanese survey 1979-1994, U.S. survey 1996-2021")

ggsave("results/turbot_historical_rpns.png", units = "in", 
       width = 8, height = 5)  

# Lengths ----

# Methods to obtain RPN-weighted lengths were updated for 2021. Records with
# length = 999 now exist to account for instances when there was catch in a
# stratum but no lengths collected. The 999 lengths ensure RPNs sum properly to
# the area-level but should not be included in the length compositions.
# Currently, the user must join on `area_view_lls` in order to filter out areas
# not used in RPN calculations (i.e. `exploitable = 1`).

lens <- sqlQuery(channel_akfin, query = ("
                select    *
                from      afsc.lls_length_rpn_by_area_all_strata
                where     species_code = '10115' and 
                          country = 'United States' and 
                          year >= 1997 and
                          council_sablefish_management_area in ('Bering Sea', 'Aleutians') and
                          length < 999
                order by  year asc
                ")) %>% 
  rename_all(tolower)

areaview <- sqlQuery(channel_akfin, query = ("
                select distinct   council_sablefish_management_area, council_management_area, 
                                  fmp_management_area, geographic_area_name, 
                                  exploitable, area_code
                from              afsc.lls_area_view
                ")) %>% 
  rename_all(tolower)

lens <- lens %>% 
  left_join(areaview, by = c("area_code", "geographic_area_name", "council_sablefish_management_area")) %>% 
  rename(fmp = council_management_area)

# Filter out samples from areas that are not used in RPN calculations
lens <- lens %>% filter(exploitable == 1)  

# RPN-weighted lengths in the BSAI
lensum <- lens %>% 
  group_by(fmp, year, length) %>% 
  summarize(rpn = sum(rpn))

ggplot(lensum, aes(x = length, y = rpn / 1e3, fill = fmp)) +
  geom_area(position = "identity", alpha = 0.25, size = 0.1) +
  facet_wrap(~ year, scales = 'free_y', ncol = 2) +
  labs(x = 'Length', y = 'RPN', fill = NULL)

ggsave("results/turbot_lengths.png", units = "in", 
       width = 14, height = 17)

# RPN-weighted lengths in the BSAI by sex - THE SEX-SPECIFIC TURBOT DATA PRIOR
# TO 2021 ARE ERRORS AND ARE ON THE LIST TO CORRECT IN 2022. I've kept the code
# in here to test that it does gets fixed.

lensum <- lens %>%
  mutate(sex = case_when(sex == 1 ~ 'male',
                         sex == 2 ~ 'female',
                         sex == 3 ~ 'unknown')) %>%
  group_by(fmp, year, sex, length) %>%
  summarize(rpn = sum(rpn))

lensum <- lensum %>% 
  mutate(sex = ifelse(year < 2021, 'unknown', sex))

ggplot(lensum, aes(x = length, y = rpn / 1e3, fill = sex)) +
  geom_area(position = "identity", alpha = 0.25, size = 0.1) +
  facet_wrap(~ year, scales = 'free_y', ncol = 2) +
  labs(x = 'Length', y = 'RPN', fill = NULL)

ggsave("results/turbot_lengthsex.png", units = "in", 
       width = 14, height = 17)
