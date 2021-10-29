# code from Bob Hall

setwd('C:/Users/Alice Carter/git/autotrophic_rivers/')

library(rstan)
library(tidyverse)

# SAM 5 lag
sink("SAM.stan")

cat("
    
    data {
    int <lower = 0> N;
    vector [N] P;
    vector [N] R;
    int <lower = 0> nweight;
     vector [nweight] alpha;
    
    }
    
    parameters {
    real <lower =0> a;
    real b;
    simplex [nweight] w; 
    real <lower=0> sigma;
    }
    
    transformed parameters{
    
    vector  [N] Pant;
    Pant[1:5]= P[1:5];


for (i in 6:N){
  vector  [nweight] Pvec;
    for(j in 1:nweight){
    Pvec[j]=w[j]*P[i-(j-1)];
    }
    Pant[i]=sum(Pvec);
    
    }

    }
    
    model {
    for ( i in 6:N){
    R[i] ~ normal(b+a*Pant[i], sigma); // likelihood
    }


    b~normal(0,5); //priors
    a~normal(0,1);
    w~dirichlet(alpha);
    }
    generated quantities{
    
    }"
    ,fill=TRUE)
sink()

# SAM 5 lag
sink("SAM_GPP.stan")

cat("
    
    data {
    int <lower = 0> N;
    vector [N] P;
    int <lower = 0> nweight;
     vector [nweight] alpha;
    
    }
    
    parameters {
    real <lower =0> a;
    real b;
    simplex [nweight] w; 
    real <lower=0> sigma;
    }
    
    transformed parameters{
    
    vector  [N] Pant;
    Pant[2:5]= P[2:5];


for (i in 6:N){
  vector  [nweight] Pvec;
    for(j in 1:nweight){
    Pvec[j]=w[j]*P[i-(j-1)];
    }
    Pant[i]=sum(Pvec);
    
    }

    }
    
    model {
    for ( i in 6:N){
    P[i] ~ normal(b+a*Pant[i], sigma); // likelihood
    }


    b~normal(0,5); //priors
    a~normal(0,1);
    w~dirichlet(alpha);
    }
    generated quantities{
    
    }"
    ,fill=TRUE)
sink()

klamath<-read.table("data/KlamMetab.txt", header=T, sep="")

klamathsv<-klamath[klamath$site=="SV",]
klamathwe<-klamath[klamath$site=="WE",]

klamathsv6<-klamathsv[klamathsv$year_numb==6,]

klamathsv6$GPP[111:112]<-c(13.2,13.5)
klamathsv6$ER[111:112]<-c(-9.5,-9.0)
klamathsv6<-klamathsv6[-180,]

sam_SV_all_dat<- list(R=-klamathsv6$ER, P=klamathsv6$GPP, N=length(klamathsv6$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_SV_all_dat, warmup = 500, iter = 1000, 
     chains = 4, cores = 4)

k <- plot(fit, pars = "w")+
  labs(title = "Klamath") +
  xlim(0,1)

# read in autotrophic data ####
dat <- read_csv('data/selected_autotrophic_rivers_daily.csv')

east <- filter(dat, sitecode == 'nwis_10133650', year >= 2012)
east$GPP <- zoo::na.approx(east$GPP, na.rm = F)
east$ER <- zoo::na.approx(east$ER, na.rm = F)
east <- filter(east, !is.na(GPP))

sam_east_all_dat<- list(R=-east$ER, P=east$GPP, N=length(east$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_east_all_dat, warmup = 500, iter = 1000, 
     chains = 4, cores = 4)
# sam_east_gpp_dat<- list( P=east$GPP, N=length(east$GPP), nweight=5, alpha=rep(2,5))
# fit <- stan(file = 'SAM_GPP.stan', data = sam_east_gpp_dat, warmup = 500, iter = 1000, 
#      chains = 4, cores = 4)



plot(east$GPP, east$ER, ylab=expression(paste("ER", "g O"[2], "m"^{-2}, "d"^{-1})), xlab=expression(paste("GPP", "g O"[2], "m"^{-2}, "d"^{-1})), pch=16, col="blue")
lines(c(2,12), c(-2,-12))


print(fit, pars=c("a", "b", "w"))

e <- plot(fit, pars = "w")+
  ggtitle("East Canyon Creek") +
  xlim(0,1)
east %>% filter(year == '2012') %>%
ggplot( aes(Date, scale(GPP) + 3)) +
  geom_line(col = 'forestgreen', size = 1.2) +
  theme_classic() +
  geom_line(aes(y = scale(discharge_m3s) - 6)) +
  geom_line(aes(y = scale(ER)), size = 1.2, col = 'sienna')+
  geom_point(aes(y = light_PAR*10/max(east$light_PAR, na.rm = T)),
             col = 'goldenrod') +
  ylab('relative values') 
pecos <- filter(dat, sitecode == 'nwis_08446500') %>%
  mutate(GPP <- zoo::na.approx(GPP, na.rm = F),
         ER <- zoo::na.approx(ER, na.rm = F))
pecos <- filter(pecos, !is.na(GPP))

sam_pecos_all_dat<- list(R=-pecos$ER, P=pecos$GPP, N=length(pecos$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_pecos_all_dat, warmup = 500, iter = 1000, 
     chains = 4, cores = 4)
p <- plot(fit, pars = "w")+
  ggtitle("Pecos River") +
  xlim(0,1)

# grand ####
grand <- filter(dat, sitecode == 'nwis_04119400') %>%
  mutate(GPP <- zoo::na.approx(GPP, na.rm = F),
         ER <- zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP), !is.na(ER))

sam_grand_all_dat<- list(R=-grand$ER, P=grand$GPP, N=length(grand$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_grand_all_dat, warmup = 500, iter = 1000, 
     chains = 4, cores = 4)
g <- plot(fit, pars = "w")+
  ggtitle("Grand River") +
  xlim(0,1)

N <- nrow(grand)
plot(grand$GPP, grand$ER)
points(grand$GPP[1:(N-1)], grand$ER[2:N], col = 2)
points(grand$GPP[1:(N-2)], grand$ER[3:N], col = 3)
points(grand$GPP[1:(N-3)], grand$ER[4:N], col = 4)
points(grand$GPP[1:(N-4)], grand$ER[5:N], col = 5)


# snake ####
snake <- filter(dat, sitecode == 'nwis_13173600') %>%
  mutate(GPP <- zoo::na.approx(GPP, na.rm = F),
         ER <- zoo::na.approx(ER, na.rm = F)) %>%
  filter(!is.na(GPP), !is.na(ER))

sam_snake_all_dat<- list(R=-snake$ER, P=snake$GPP, N=length(snake$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_snake_all_dat, warmup = 500, iter = 1000, 
     chains = 4, cores = 4)
s <- plot(fit, pars = "w")+
  ggtitle("Snake River") +
  xlim(0,1)

ggplot(snake, aes(Date, scale(GPP) + 3)) +
  geom_line(col = 'forestgreen', size = 1.2) +
  theme_classic() +
  geom_line(aes(y = scale(discharge_m3s) - 2)) +
  geom_line(aes(y = scale(ER)), size = 1.2, col = 'sienna')+
  geom_point(aes(y = light_PAR*10/max(snake$light_PAR, na.rm = T)),
             col = 'goldenrod') +
  ylab('relative values') 
 

plot(1,1, type = 'n')  
         
legend('topright',
       legend = c('productivity', 'respiration', 'discharge', 'light'),
      col = c('forestgreen', 'sienna', 'black', 'goldenrod'),
         pch = c(NA, NA, NA, 20),
         lty = c(1,1,1, NA), lwd = c(2, 2, 1, NA))

ggpubr::ggarrange(e, p, g, s)
Pant<-summary(fit, pars = c("Pant"), probs = c(0.5))$summary

par(mfcol = c(1,2), mai=c(0.8,0.7,0.1,0.1), mgp=c(2,1,0))
plot(east$GPP, east$ER, ylab=expression(paste("ER", "g O"[2], "m"^{-2}, "d"^{-1})), 
     xlab=expression(paste("GPP", " g O"[2], "m"^{-2}, "d"^{-1})), pch=16, col="blue")

plot(Pant[,1], east$ER, ylab="", xlab=expression(paste("Pant", " g O"[2], "m"^{-2}, "d"^{-1})), pch=16, col="green")


# Gleon data ####
gleon <- read_csv('data/GLEON_metabolism_export.csv')
mendota <- gleon %>%
  filter(lakeName == 'Mendota') %>%
  select(solarDay, year, GPP = GPP_mgO2m2, ER = ER_mgO2m2)%>%
  filter(!is.na(GPP))

sam_men_all_dat<- list(R=mendota$ER, P=mendota$GPP, N=length(mendota$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_men_all_dat, warmup = 500, iter = 1000, 
            chains = 4, cores = 4)
plot(fit, pars = "w")
traceplot(fit)

acton <- gleon %>%
  filter(lakeName == 'Acton') %>%
  select(solarDay, year, GPP = GPP_mgO2m2, ER = ER_mgO2m2)%>%
  filter(!is.na(GPP))

sam_act_all_dat<- list(R=acton$ER, P=acton$GPP, N=length(acton$GPP), nweight=5, alpha=rep(1,5))
fit <- stan(file = 'SAM.stan', data = sam_act_all_dat, warmup = 500, iter = 1000, 
            chains = 4, cores = 4)
plot(fit, pars = "w")

ggplot(snake, aes(Date, GPP)) + 
  # geom_point() +
  geom_line(aes(Date, (discharge_m3s)))
