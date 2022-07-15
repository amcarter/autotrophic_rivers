# Wrappers and helpers to prep data and run SAM models 
library(tidyverse)
library(zoo)

# prep data ####
# for now, I am just filling in missing data and cutting the ts to remove
# leading and ending NAs. This will have to be addressed later XXXX
prep_sam_dat <- function(dat, nweight = 5){
  d <- dat %>%
    mutate(GPP = zoo::na.approx(GPP, na.rm = F),             
           ER = zoo::na.approx(ER, na.rm = F)) %>%
    filter(!is.na(GPP))
  
  sam_dat <- list(R = -d$ER, P = d$GPP, N = length(d$GPP), 
                  nweight = nweight, alpha = rep(1, nweight))
  
  return(sam_dat)
  
}

# run models ####
ast %>%