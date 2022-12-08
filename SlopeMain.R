rm(list=ls())
set.seed( 2011 )
library("adapr")
#install.packages("adapr")
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
#library(epicalc)
library(epiDisplay)
library("readxl")
library('ggpubr')
lirary('ggplot2')


setwd("/Users/rabiul/fwdcodeusedinanalysis")
#redcap_data_adnin
redcap <-  read.csv("RedCapUp2.csv",stringsAsFactors = T) 

redcapegfr<-redcap%>%
  dplyr::select(patientid,egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,
                sex_v3y0,albuminuria_v3y0,adnine_v3y0)



###################################LowerTertile for Male eGFR selection for Normal
EGFR<-redcapegfr%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)%>%
  mutate(quartile = ntile(adnine_v3y0, 2))%>%
  filter(quartile==1)

dim(EGFR)

#write.csv(x=EGFR, file="LowTertile_Male_Normal.csv")



###################################UpTertile for Male eGFR selection for Normal
EGFR<-redcapegfr%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)%>%
  mutate(quartile = ntile(adnine_v3y0, 2))%>%
  filter(quartile==2)

dim(EGFR)

#write.csv(x=EGFR, file="UPTertile_Male_Normal.csv")




###################################LowerTertile for All Patients eGFR selection for Normal
EGFR<-redcapegfr%>%
  #filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)%>%
  mutate(quartile = ntile(adnine_v3y0, 2))%>%
  filter(quartile==1)

dim(EGFR)

#write.csv(x=EGFR, file="LowTertile_Normal.csv")


###################################UperTertile for All Patients eGFR selection for Normal
EGFR<-redcapegfr%>%
  #filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)%>%
  mutate(quartile = ntile(adnine_v3y0, 2))%>%
  filter(quartile==2)

dim(EGFR)

#write.csv(x=EGFR, file="UpTertile_Normal.csv")


###################################All Patients eGFR selection for Normal
EGFR<-redcapegfr%>%
  #filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
#write.csv(x=EGFR, file="AllEGFR_Normal.csv")


###################################All Patients eGFR selection for Micro
EGFR<-redcapegfr%>%
  #filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
#write.csv(x=EGFR, file="AllEGFR_Micro.csv")


###################################All Patients eGFR selection for Macro
EGFR<-redcapegfr%>%
  #filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0!=1)%>%
  filter(albuminuria_v3y0!=2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
write.csv(x=EGFR, file="AllEGFR_Macro.csv")

###################################Male eGFR selection for Normal
EGFR<-redcapegfr%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
#write.csv(x=EGFR, file="MaleEGFR_Normal.csv")

###################################Male eGFR selection for Micro
EGFR<-redcapegfr%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)

#write.csv(x=EGFR, file="MaleEGFR_Micro.csv")


###################################Male eGFR selection for Macro
EGFR<-redcapegfr%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0!=1)%>%
  filter(albuminuria_v3y0!=2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)

#write.csv(x=EGFR, file="MaleEGFR_Macro.csv")


#############################################FeMale eGFR selection for Normal

EGFR<-redcapegfr%>%
  filter(sex_v3y0==2)%>%
  filter(albuminuria_v3y0==1)%>%
  #filter(albuminuria_v3y0!=2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
write.csv(x=EGFR, file="FemaleEGFR_Normal.csv")


#############################################Female eGFR selection for Micro

EGFR<-redcapegfr%>%
  filter(sex_v3y0==2)%>%
  filter(albuminuria_v3y0==2)%>%
  #filter(albuminuria_v3y0!=2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
#write.csv(x=EGFR, file="FemaleEGFR_Micro.csv")




#############################################Female eGFR selection for Macro

EGFR<-redcapegfr%>%
  filter(sex_v3y0==2)%>%
  filter(albuminuria_v3y0!=1)%>%
  filter(albuminuria_v3y0!=2)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)
write.csv(x=EGFR, file="FemaleEGFR_Macro.csv")




################For matlab data slope analysis
redcap <-  read.csv("RedCapUp2.csv",stringsAsFactors = T) %>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>% #ArrDelay/max(ArrDelay)
  mutate(adnine_n=abs(scale(adnine_v3y0)))%>% 
  mutate(adnine_w=adnine_v3y0-mean(adnine_v3y0,na.rm = TRUE)/sd(adnine_v3y0,na.rm = TRUE))%>%#-mean(adnine_v3y0)/sd(adnine_v3y0,na.rm = FALSE))%>%
  dplyr::select(patientid,egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,egfr_cric_v19y8,sa_allc_ckdepi_renal1_v5y1, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_v3y0,adnine_n,adnine_w)%>%
  #rename_all(function(x) gsub('ientid', '', x))%>%
  na.omit()


redcap_data_adnin<-redcap%>%
  dplyr::select(patientid,sa_allc_ckdepi_renal1_v5y1, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_n,adnine_w)%>%
  na.omit()

#Create event sumarry
redcap_data3<-within(redcap_data_adnin,{
  #Create event variable
  event <- as.numeric(sa_allc_ckdepi_renal1_v5y1 %in% c(1,2)) #############
  #Create a survival vectro
  surv<-Surv(sa_allc_time_esrd_v3y0,event)
})

#Adenine+Clinical parameters
logity1 <- glm(event ~ age_integer_v3y0+hemoglobin_a1c_v3y0+map_v3y0+bmi_v3y0+egfr_cric_v3y0+log_acr+adnine_v3y0+adnine+adnine_n, data = redcap_data3, family = binomial)
lroc(logity1, graph = T)$auc

#Clinical parameters
logity2 <- glm(event ~ age_integer_v3y0+hemoglobin_a1c_v3y0+map_v3y0+bmi_v3y0+egfr_cric_v3y0+log_acr, data = redcap_data3, family = binomial)
lroc(logity2, graph = T)$auc

anova(logity1, logity2, test = "LRT")
