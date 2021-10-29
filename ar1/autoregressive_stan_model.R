# Write a linear autoregressive state space model of GPP driven by 
#   light for an autotrophic river using stan.
# 10/8/2021
# A Carter

# setup ####
library(tidyverse)
library(lme4)
library(rstan)
library(gdata)
library(bayesplot)

setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')
# read in autotrophic stream data 
# select the five most autotrophic (highest GPP) streams
auto <- read_csv('data/autotrophic_siteyears_daily_filled_2.csv')

annual <- auto %>%
  group_by(sitecode, year) %>%
  summarize(NEP = sum(NEP_C, na.rm = T),
            GPP = sum(GPP_C_filled, na.rm = T),
            par_av = sum(!is.na(light_PAR))/365,
            prod_quo = NEP/GPP) %>%
  arrange(desc(prod_quo))

hist(annual$prod_quo)
ggplot(annual, aes(GPP, NEP, col = NEP/GPP))+
  geom_point(size = 2)

# there are two groups here that I think are interesting in terms of autotrophy:
# 1. High GPP, over 1500 gC/m2/y
# 2. High fraction of GPP converted to NEP, over 30%

# High GPP rivers:
high <- annual %>%
  filter(GPP > 1500 |prod_quo > 0.3, par_av > 0.5) %>%
  arrange(desc(GPP))

auto <- auto %>% filter(sitecode %in% unique(high$sitecode))
ann_auto <- auto %>% 
  group_by(sitecode, year) %>%
  summarize(GPP = sum(GPP_C_filled),
            NEP = sum(NEP_C),
            NEP_GPP = NEP/GPP) %>%
  ungroup() 

plot(ann_auto$GPP, ann_auto$NEP_GPP)
identify(ann_auto$GPP, ann_auto$NEP_GPP, labels = ann_auto$sitecode)

# five of the six site years are from the same river,
#   the Pecos River. It looks like they moved the sensor after 3 years from
#   Girvin to Sheffield, TX. Not sure what to do with that yet because the 
#   years at Girvin (2012-14) all have a huge jump in Q in late fall but
#   not the years at Sheffield (2015-16)
pecos <- auto %>% 
  filter(sitecode == high$sitecode[2])

snake <- auto %>%  
  filter(sitecode == high$sitecode[1])

# look at the data
# the snake river desn't even have a complete year. Maybe it does at some
# other sites but for now I will work on Pecos
plot(pecos$Date, pecos$GPP_C_filled, type = 'l')
par(new = T)
plot(pecos$Date, log(pecos$discharge_m3s), col = 2, type = 'l')
ggplot(pecos, aes(DOY, log(discharge_m3s), col = year))+
  geom_point()
ggplot(pecos, aes(DOY, GPP_C_filled, col = year))+
  geom_point()
ggplot(pecos, aes(DOY, light_PAR, col = year)) +
  geom_point()
ggplot(pecos, aes(light_PAR, GPP_C_filled, col = log(discharge_m3s))) +
  geom_point()

# We need the unfilled timeseries for this analysis
pecos_filled <- pecos
dat <- read_csv('../loticlentic_synthesis/data/powell_data_import/compiled_daily_model_results.csv')

pecos <- dat %>%
  select(sitecode = site_name, Date = date, 
         starts_with(c('GPP', 'ER'))) %>%
  right_join(pecos_filled, by = c('sitecode', 'Date')) %>%
  arrange(Date)

plot(pecos$Date, pecos$GPP, pch = 20)
lines(pecos$Date, pecos$GPP_C_filled, col = 2)

plot(density(pecos$GPP, na.rm = T))
plot(density(pecos$light_PAR, na.rm = T))
# remove values for GPP less than zero
pecos <- pecos %>% 
  mutate(across(starts_with('GPP'), 
         ~ case_when(GPP < 0 ~ NA_real_,
                     TRUE ~ .))) 

# get unfilled data for all eight of the autotrophic rivers selected above:
aa <- dat %>%
  select(sitecode = site_name, Date = date,
         starts_with(c('GPP', 'ER'))) %>%
  right_join(auto, by = c('sitecode', 'Date')) %>%
  arrange(sitecode, Date)

aa <- aa %>%
  mutate(across(starts_with('GPP'),
                ~case_when(GPP < 0 ~ NA_real_,
                           TRUE ~ .)),
         across(starts_with('ER'),
                ~case_when(ER > 0 ~ NA_real_,
                           TRUE ~ .))) %>%
  select(-NEP_C, -ends_with('filled'))

write_csv(aa, 'data/selected_autotrophic_rivers_daily.csv')
# try a linear model:

pecos_lm <- pecos %>%
  mutate(GPP_pre = c(pecos$GPP[1], pecos$GPP[1:(nrow(pecos)-1)]))

mm <- lm(GPP ~ GPP_pre + light_PAR, data = pecos_lm)
summary(mm)
  
# Stan model ####
# prepare data for a stan model 
# stan models cannot accommodate NAs, but filling data is bad. so I'm just 
#   going to remove rows with missing data

# Model equation:
# GPP_t = a + b0 * GPP_t-1 + b1 * PAR_t + sigma_proc 
# GPP_obs_t = GPP_t + sigma_obs

pecos_nona <- pecos_lm %>%
  select(GPP, GPP_pre, light_PAR) %>%
  filter(!is.na(GPP),
         !is.na(GPP_pre),
         !is.na(light_PAR)) %>%
  mutate(across(everything(), scale))
GPP <- pecos_nona$GPP[,1]
GPP_pre <- pecos_nona$GPP_pre[,1] 
light_par <- pecos_nona$light_PAR[,1]
N <- length(GPP)
stan_dat <- list(N = N, GPP = GPP, GPP_pre = GPP_pre, light_par = light_par)

write("// Stan model for simple autoregressive linear regression

data {
 int < lower = 1 > N; // Sample size
 vector[N] GPP_pre; // Predictor1
 vector[N] light_par; // Predictor2
 vector[N] GPP; // Outcome
}

parameters {
 real alpha; // Intercept
 real < lower = 0, upper = 1> beta1; 
 real beta2; 
 real < lower = 0 > sigma_obs;
 real < lower = 0 > sigma_proc;
}

model {
 GPP ~ normal(alpha + beta1 * GPP_pre + beta2 * light_par + sigma_proc, sigma_obs);
 alpha ~ normal(0,1);
 beta1 ~ uniform(0,1);
 beta2 ~ normal(0,1);
 sigma_obs ~ normal(0,1);
 sigma_proc ~ normal(0,1);
}

generated quantities {
} // The posterior predictive distribution",

"src/ar1/ar1_model.stan")
 

stanc('src/ar1/ar1_model.stan')
ar1_mod <- 'src/ar1/ar1_model.stan'

fit <- stan(file = ar1_mod, data = stan_dat, warmup = 500, iter = 10000, 
            chains = 4, cores = 2, thin = 1)
fit

posterior <- extract(fit)
str(posterior)


traceplot(fit)
stan_dens(fit)
plot(fit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon")


acf(GPP)
pacf(GPP)

