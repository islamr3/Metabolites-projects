rm(list=ls())
set.seed( 2011 )
library("adapr")
library('survival')
library('survminer')
library('dplyr')
library('CBPS')
library('purrr')
library('gtsummary')
library('gt')
library('cobalt')
library("officer")
library("rvg")
library('lme4')
library('lmerTest')
library('effects')
library('emmeans')
library("rstudioapi") 
library(survival)
library(rms)
library(survivalROC)
library(Hmisc)
library(epicalc)
library(epiDisplay)
library("readxl")

#Load pbc dataset
data(pbc)

#Create event sumarry
pbc1 <- within(pbc, {
  ## transplant (1) and death (2) are considered events, and marked 1
  event <- as.numeric(status %in% c(1,2))
  
  ## Create a survival vector
  Surv <- Surv(time, event)
})


#pbc2 <- within(pbc1, {
#  outcome2yr <- NA
#  outcome2yr[(event == 1) & (time <= 2 * 365)] <- 1 # event+ within two years
#  outcome2yr[(event == 0) | (time  > 2 * 365)] <- 0 # otherwise
#})

## Null model
coxph.null<- coxph(formula = Surv ~ 1, data = pbc1)
## Model with age and sex
coxph.age.sex<- coxph(Surv ~ age + sex, data = pbc1)
## Model with age, sex, and albumin
coxph.age.sex.albumin  <- coxph(Surv ~ age + sex + albumin, data = pbc1)

## Put linear predictors ("lp") into pbc dataset
pbc1$lp.null            <- predict(coxph.null, type = "lp")
pbc1$lp.age.sex         <- predict(coxph.age.sex, type = "lp")
pbc1$lp.age.sex.albumin <- predict(coxph.age.sex.albumin, type = "lp")


## Define a function
fun.survivalROC <- function(lp, t) {
  res <- with(pbc1,
              survivalROC(Stime        = time,
                          status       = event,
                          marker       = get(lp),
                          predict.time = t,
                          method       = "KM"))       # KM method without smoothing
  
  ## Plot ROCs
  with(res, plot(TP ~ FP, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC)))
  abline(a = 0, b = 1, lty = 2)
  res}



## 2 x 5 layout
layout(matrix(1:10, byrow = T, ncol = 5))

## Model with age and sex
res.survivalROC.age.sex <- lapply(1:10 * 365.25, function(t) {
  fun.survivalROC(lp = "lp.age.sex", t)})
