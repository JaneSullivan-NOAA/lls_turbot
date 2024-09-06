# Killer whale depredation
# Contact: jane.sullivan@noaa.gov
# Updated May 2022
# Generated using  R version 4.1.2

# set up ----
libs <- c("tidyverse", "RODBC")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

theme_set(theme_minimal(base_size = 10))

username_akfin <- '' # fill in your akfin creds
password_akfin <- ''

channel_akfin <- odbcConnect("akfin", uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

# Raw catch query ----

# BSAI LLS catch by skate for turbot including zero catches (BS = odd years, AI
# = even years). Full time series starts in 1997. 

catch <- sqlQuery(channel_akfin, query = paste0("
                select    *
                from      afsc.lls_catch_summary_with_nulls_mv
                where     country = 'United States' and
                          year > 1996 and 
                          species_code = '10115' and
                          ((council_sablefish_management_area in 'Bering Sea' and mod(year, 2) <> 0) or
                          (council_sablefish_management_area in ('Aleutians') and mod(year, 2) = 0))
                order by  year asc
                ")) %>% 
  rename_all(tolower) 

write_csv(catch, 'data/raw_turbot_catch_with_zeros.csv')
catch <- read_csv('data/raw_turbot_catch_with_zeros.csv')
glimpse(catch)

# column descriptions start on p 6:
# https://akfinbi.psmfc.org/analyticsRes/Documentation/Database_Background_Instructions_AKFIN_20210915.pdf

# rpn_filter: Used to filter out skates for relative population number (RPN),
# relative population weight (RPW), and CPUE calculations. A “k” or “K” is used
# to code for killer whale depredation by each hachi. A “g” is used when the
# skate of gear was not fishing effectively

nrow(catch)

catch <- catch %>% 
  mutate(fishing_event_id = paste0(cruise_number, '_', station_number, '_', hachi)) 

# Sets depredated by KW ----
depred <- catch %>% 
  group_by(year, council_sablefish_management_area) %>% 
  dplyr::summarize(nsets = length(unique(fishing_event_id))) %>% 
  left_join(catch %>% 
              filter(rpn_filter %in% c('k', 'K')) %>% 
              group_by(year, council_sablefish_management_area) %>% 
              dplyr::summarize(nsets_depred = length(unique(fishing_event_id)))) %>% 
  mutate(nsets_depred = replace_na(nsets_depred, 0),
         propn_sets_depred = nsets_depred / nsets)

depred %>% 
  ggplot(aes(x = year, y = propn_sets_depred)) +
  geom_line() +
  geom_point() +
  facet_wrap(~council_sablefish_management_area ) +
  theme_bw(base_size = 13) +
  labs(y = "Proportion sets depredated")

ggsave("results/turbot_kwdepred_sets.png", units = "in",
       width = 7, height = 5, bg = 'white')

sum <- depred %>% 
  rename(fmp = council_sablefish_management_area) %>% 
  left_join(read_csv('data/bsai_turbot_rpns.csv')) %>% 
  select(year, fmp, propn_sets_depred, rpn) %>% 
  group_by(fmp) %>% 
  mutate(std_rpn = (rpn - mean(rpn)) / sd(rpn),
         std_depred = ifelse(is.na(propn_sets_depred), NA,
                             (propn_sets_depred - mean(propn_sets_depred, na.rm = TRUE)) / 
                               sd(propn_sets_depred, na.rm = TRUE)))

label <- sum %>% 
  group_by(fmp) %>% 
  summarize(cor = cor(std_rpn, std_depred, method = 'pearson'))  %>% 
  ungroup() %>% 
  mutate(year = 2009, 
         value = 2,
         corlab = paste0("Pearson's = ", round(cor, 2)))

sum %>% 
  rename(`Standardized KW depred` = std_depred,
         `Standardized turbot RPN` = std_rpn) %>% 
  pivot_longer(cols = c(`Standardized KW depred`, `Standardized turbot RPN`)) %>%
  ggplot(aes(x = year, y = value)) +
  geom_line(aes(col = name, lty = name)) +
  geom_point(aes(col = name, shape = name)) +
  geom_text(data = label, 
            aes(x = year, 
                y = value,
                label = corlab),
            size = 2.5) +
  facet_wrap(~fmp, ncol = 1) +
  labs(x = NULL, y = NULL,
       col = NULL, shape = NULL, lty = NULL) +
  ggthemes::scale_color_colorblind() +
  theme_bw(base_size = 13) 
  
ggsave("results/turbot_kwdepred_correlation.png", units = "in",
       width = 7, height = 6, bg = 'white')

