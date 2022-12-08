#Requirements
requirements <- c("survival", "survminer")

#CRAN repository
repos <- "http://cran.us.r-project.org"

#install and load missing requirements
for (pack in requirements){
  if( !is.element(pack, .packages(all.available = TRUE)) ) {
    #install.packages(pack, repos = repos)
    install.packages(pack)
  }
  library(pack, character.only = TRUE)
}

#Random number seed
rm(list=ls())
set.seed(10079)

#get data
rawdata <- read.csv("DOD-AA-3H-BL-SUR-ANA-2021.csv")

#define variables
covariates <- names(rawdata[3:23])


# treat missing below detection
rawdata[is.na(rawdata)] <- 0.

#treat zeros with half of the minimun
rawdata[covariates] <- lapply(rawdata[covariates], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/2))


