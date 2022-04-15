# Functions for metabolism models
# A Carter
# 4/2022

# Data Prep Functions ####
calc_litter <- function(annual_litter, LAI){
  LAI <- zoo::na.approx(LAI, na.rm = F)
  diffLAI <- c(0, diff(LAI))
  diffLAI <- case_when(diffLAI < 0 ~ diffLAI,
                       TRUE ~ 0)
  litter <- annual_litter * diffLAI/sum(diffLAI, na.rm = T)
  
  return(litter)
}
