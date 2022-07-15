# Find rivers in the Powell center database that are autotrophic at the annual 
# timescale, complie all daily data 

# A Carter
# 2021-09

setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
library(tidyverse)
library(lubridate)

# read in compiled Powell center data:
dat <- 
  read_csv('../loticlentic_synthesis/data/powell_data_import/compiled_daily_model_results.csv')

# Remove negative GPP, posative ER, and Rhat > 1.05
filtered <- dat %>%
  mutate(GPP = case_when(GPP.Rhat >1.05 ~ NA_real_,
                         GPP >= 0 ~ GPP,
                         TRUE ~ NA_real_),
         ER = case_when(ER.Rhat >1.05 ~ NA_real_,
                        ER <= 0 ~ ER,
                        TRUE ~ NA_real_))

# only keep years that have at least 25% of the days
exc <- filtered %>%
  mutate(year = year(date)) %>%
  group_by(site_name, year) %>%
  summarize(GPP_days = sum(!is.na(GPP))/365,
            ER_days = sum(!is.na(ER))/365) %>%
  ungroup() %>%
  filter(GPP_days < 0.25 | ER_days < 0.25) %>%
  mutate(siteyear = paste0(site_name, year))
  
filtered <- filtered %>%
  mutate(year = year(date),
         siteyear = paste0(site_name, year)) %>%
  filter(!(siteyear %in% unique(exc$siteyear)) ) %>%
  select(-siteyear)

write_csv(filtered, 'data/filtered_powell_data_25.csv') 

# This is the most updated filled dataset taht was used in the SP paper, 
# I got it from Mike V 11/2021:

phil_data <- read_csv('../loticlentic_synthesis/data/phil_powell_data/lotic_gap_filled_unlisted.csv')
tmp <- phil_data %>%
  group_by(sitecode) %>%
  mutate(GPP = zoo::na.approx(GPP, x = Date, na.rm = F ))
# linear interpolation does not match Phil's gap filling. This needs to be 
# addressed before deciding how to get annual NEP if I am going to use that.

plot(filtered$ER, filtered$GPP)
points(phil_data$ER, phil_data$GPP, col = alpha(2, 0.3), pch = 20)

# readRDS('../oxygen_timeseries/data/gapPhilled_data.rds')
# is there more metadata somewhere?
site_dat <- 
  readRDS('../loticlentic_synthesis/data/phil_powell_data/lotic_site_info_full.rds') %>%
  rename(sitecode = Site_ID)

#### start here ####
# Summarize by site years
annual <- dat %>% 
  mutate(year = year(date)) %>%
  group_by(site_name, year) %>%
  summarize(across(starts_with(c('GPP', 'ER')), sum, na.rm = T),
            n = n()) %>%
  ungroup() %>%
  select(-ends_with(c('Rhat', 'n_eff')))

f_ann <- filled_dat %>% 
  group_by(sitecode, year) %>%
  summarize(across(c('GPP_C_filled', 'ER_C_filled'), sum, na.rm = T),
            n = n()) %>%
  ungroup() %>%
  mutate(NEP_C_filled = GPP_C_filled + ER_C_filled)

jpeg('figures/distribution_annual_NEP.jpeg')

  plot(density(f_ann$NEP_C_filled, na.rm = T), xlim = c(-2000, 500),
       main = "Annual NEP of powell center site years",
       xlab = "g O2/m2/y")
  qs = quantile(f_ann$NEP_C_filled, probs=c(0.05, 0.95))
  dd = density(f_ann$NEP_C_filled, na.rm = T)
  ddo = order(dd$x)
  xdens = dd$x[ddo]
  ydens = dd$y[ddo]
  xdens_lt = xdens[xdens <= qs[1]]
  ydens_lt = ydens[xdens <= qs[1]]
  polygon(c(xdens_lt, rev(xdens_lt)), c(ydens_lt, rep(0, length(ydens_lt))),
          col='gray', border='gray')
  xdens_ut = xdens[xdens >= qs[2]]
  ydens_ut = ydens[xdens >= qs[2]]
  polygon(c(xdens_ut, rev(xdens_ut)), c(ydens_ut, rep(0, length(ydens_ut))),
          col='lightgreen', border='lightgreen')
  lines(dd)
  legend('topleft', legend=c('heterotrophic group', 'autotrophic group'), 
         bty='n', fill = c('gray', 'lightgreen'), text.font=2, cex=1)
  
dev.off()

dd <- left_join(f_ann, site_dat, by = c('sitecode')) %>%
  mutate(trophic = case_when(NEP_C_filled > 0 ~ 'autotrophic', 
                             NEP_C_filled < 0 ~ 'heterotrophic', 
                             TRUE ~ NA_character_)) 

write_csv(dd, 'data/filled_annual_metabolism_with_trophic_state.csv')

auto <- dd %>% select(-n, -Source, -ends_with('filled')) %>%
  right_join(filled_dat, by = c('sitecode', 'year')) %>%
  filter(trophic == 'autotrophic')

unique(auto$sitecode)
write_csv(auto, 'data/autotrophic_siteyears_daily_filled_data.csv')


# make dataframe with the top 5% autotrophic and heterotrophic steams

ahrivs <- dd %>% 
  filter(NEP_C_filled >= qs[2] | NEP_C_filled <= qs[1]) %>%
  select(-n, -Source, -ends_with('filled')) %>%
  left_join(filled_dat, by = c('sitecode', 'year'))

ahrivs %>% group_by(sitecode, year, trophic) %>%
  summarize(NEP = sum(GPP_C_filled) + sum(ER_C_filled)) %>%
  distinct() %>% arrange(sitecode, year)

write_csv(ahrivs, 'data/autohetero_quantiles_siteyears_daily_filled_data.csv')
