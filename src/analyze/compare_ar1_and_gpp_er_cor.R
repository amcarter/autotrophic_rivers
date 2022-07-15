library(tidyverse)
library(lubridate)
setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
dat <- read_csv('data/filled_annual_metabolism_with_trophic_state.csv')
ah <- read_csv('data/autohetero_quantiles_siteyears_daily_filled_data.csv')
sp_full <- read_csv('C:/Users/Alice Carter/git/loticlentic_synthesis/data/sp_full.csv') %>%
  mutate(NEP_C_filled = GPP_C_filled + ER_C_filled)

list <- dat %>% select(sitecode, year, trophic, GPP_C_filled, ER_C_filled) %>%
  distinct() %>%
  mutate(AR1 = NA_real_)

for( i in 750:nrow(list)){
  d <- sp_full %>% filter(sitecode == list$sitecode[i], Year == list$year[i])
  a <- ar(d$NEP_C_filled, order.max = 1)
  list$AR1[i] <- a$ar
}
list$cor <- NA_real_

for( i in 1:nrow(list)){
  d <- sp_full %>% filter(sitecode == list$sitecode[i], Year == list$year[i])
  list$cor[i] <- cor(d$GPP_C_filled, -d$ER_C_filled)
  }

list$NEP <- list$GPP_C_filled + list$ER_C_filled
plot(list$NEP, list$cor)
plot(list$GPP_C_filled, list$cor)

ggplot(list, aes(cor, AR1, col = log(GPP_C_filled))) +
  geom_point(size = 2) +
  xlab('GPP:ER correlation')+
  ylab('ar1 for NEP') +
  scale_color_gradientn(colors = rev(terrain.colors(3)))

sp_full %>% filter(sitecode == 'nwis_10133800') %>%
  ggplot(aes(DOY, GPP_C_filled)) +
  geom_point() +
  geom_point(aes(y = ER_C_filled))+
  facet_wrap(~Year, scales = 'free')

