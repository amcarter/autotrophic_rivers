# Evaluate seasonality of autotrophy
# A Carter
# 9/26/2021

# setup ####
library(tidyverse)
library(lubridate)
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
dat <- read_csv('data/autotrophic_siteyears_daily_filled_2.csv')
sp_full <- read_csv('C:/Users/Alice Carter/git/loticlentic_synthesis/data/sp_full.csv')
# group data by month, determine:
#   1. which month has peak autotrophy and heterotrophy
#   2. what the cumulative NEP is each month
#   3. How the seasonality of flow/lenticity conditions compares

# Compare the timing of peak metabolic fluxes across the heterotrophic and autotrophic site years

sitemonths = dat %>%
  rename(GPP = GPP_C_filled,
         ER = ER_C_filled,
         NEP = NEP_C) %>%
  mutate(month = month(Date)) %>%
  group_by(sitecode, year, month) %>%
  summarize(across(c(GPP, ER, NEP), mean, na.rm = TRUE),
            .groups = 'drop') %>%
  group_by(sitecode, month) %>%
  summarize(across(c(GPP, ER, NEP), #mean, na.rm = TRUE),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))),
            n = n(),
            .groups = 'drop')

sitemonths_all = sp_full %>%
  filter(!(sitecode %in% unique(dat$sitecode))) %>%
  rename(GPP = GPP_C_filled,
         ER = ER_C_filled) %>%
  mutate(month = month(Date), 
         NEP = GPP + ER) %>%
  group_by(sitecode, Year, month) %>%
  summarize(across(c(GPP, ER, NEP), mean, na.rm = TRUE),
            .groups = 'drop') %>%
  group_by(sitecode, month) %>%
  summarize(across(c(GPP, ER, NEP), #mean, na.rm = TRUE),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE))),
            n = n(),
            .groups = 'drop')


maxmonths_etc_gpp = sitemonths %>%
  group_by(sitecode) %>%
  filter(GPP_mean == max(GPP_mean),
         n == max(n)) %>%
  summarize(GPP_mean = first(GPP_mean),
            n = first(n),
            month = first(month),
            .groups = 'drop')
maxmonths_etc_gpp_all = sitemonths_all %>%
  group_by(sitecode) %>%
  filter(GPP_mean == max(GPP_mean),
         n == max(n)) %>%
  summarize(GPP_mean = first(GPP_mean),
            n = first(n),
            month = first(month),
            .groups = 'drop')
maxmonths_etc_nep = sitemonths %>%
  group_by(sitecode) %>%
  filter(NEP_mean == max(NEP_mean),
         n == max(n)) %>%
  summarize(NEP_mean = first(NEP_mean),
            n = first(n),
            month = first(month),
            .groups = 'drop')
maxmonths_etc_nep_all = sitemonths_all %>%
  group_by(sitecode) %>%
  filter(NEP_mean == max(NEP_mean),
         n == max(n)) %>%
  summarize(NEP_mean = first(NEP_mean),
            n = first(n),
            month = first(month),
            .groups = 'drop')

tt <- table(maxmonths_etc_gpp$month)
maxmonths_gpp = c(tt[1],'2'=0, tt[2:7], '9'=0, tt[8], '11'=0, tt[9])
maxmonths_gpp_all = c('1' = 0, table(maxmonths_etc_gpp_all$month), '12' = 0)

maxmonths_etc_er = sitemonths %>%
  group_by(sitecode) %>%
  filter(ER_mean == max(ER_mean),
         n == max(n)) %>%
  summarize(ER_mean = first(ER_mean),
            n = first(n),
            month = first(month),
            .groups = 'drop') 

tt <- table(maxmonths_etc_er$month)
maxmonths_er = c(tt[1:3],'4'=0, tt[4], '6'=0, tt[5:6], '9'=0, tt[7:9])

maxmonths = matrix(c(maxmonths_gpp, maxmonths_gpp_all),
                       nrow = 2, byrow = TRUE,
                       dimnames = list(c('GPP', 'ER'),
                                       month.abb))
jpeg(width=8, height=8, units='in', res=300, quality=100, type='cairo',
     filename='figures/peak_months.jpg')

plot(density(maxmonths_etc_gpp_all$month), main = 'Timing of peak GPP',
     xlab = 'month', xaxt = 'n', ylim = c(0.01, 0.26), lty = 2, bty = 'n')
axis(1, at = seq(1:12), labels = month.abb)
lines(density(maxmonths_etc_gpp$month))
legend(x = 10, y = 0.23, legend=c('autotrophic', 'heterotrophic'),
       lty = c(1,2), bty = 'n', border = NA, 
       xpd = T)

dev.off()

jpeg(width=8, height=8, units='in', res=300, quality=100, type='cairo',
     filename='figures/peak_months_nep.jpg')

plot(density(maxmonths_etc_nep_all$month), main = 'Timing of peak NEP',
     xlab = 'month', xaxt = 'n', ylim = c(0, 0.19), xlim = c(0.5,12.5), lty = 2)
axis(1, at = seq(1:12), labels = month.abb)
lines(density(maxmonths_etc_nep$month))
legend(x = 10, y = 0.17, legend=c('autotrophic', 'heterotrophic'),
       lty = c(1,2), bty = 'n', border = NA, 
       xpd = T)

dev.off()
