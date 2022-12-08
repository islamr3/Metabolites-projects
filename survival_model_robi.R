rm(list=ls())
set.seed( 2011 )
library(survival)
library(survminer)
library(lubridate)
library(dplyr)
library("car")
library(cmprsk)
library(tidyverse)

library(forestplot)
library(magrittr)
library(sjPlot)


redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)

redData<-redcap_data_adnin%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0)%>%
  na.omit()

redData1<-redData%>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))


#Create event sumarry
data1<-within(redData1,{
  #Create event variable
  event <- as.numeric(sa_allc_ckdepi_renal1_v3y0 %in% c(1,2)) #############
  #Create a survival vectro
  surv<-Surv(sa_allc_time_esrd_v3y0,event)
})

#survival model
Surdata2<-data1%>%
  mutate(quartile=ntile(adnine,2))

Ndata2<-Surdata2%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  na.omit()


sfit <- survfit(surv~quartile, data=Ndata2)
ggsurvplot(sfit, fun = "cumhaz",pval = TRUE)
ggsurvplot(sfit, fun = "cumhaz",pval = TRUE,legend.title = "Quartile",
           legend.labs = c("Lower", "Upper"),palette = c("blue", "black"))



cox <- coxph(surv ~ age + ha1c + map + bmi + age + egfr+adnine , data = Ndata2)
plot(cox, xlab="Days")

ggsurvplot(cox)


summary(cox)
cox_fit <- survfit(cox)
##https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)
redData<-redcap_data_adnin%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  #dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0)%>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_ckdepi_time_renal3_v3y0,sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0)%>%
  filter(sa_allc_ckdepi_renal1_v3y0 !=9)%>% 
  na.omit()

redData1<-redData%>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>%
  mutate(quartile=ntile(adnine,2))



#Create event sumarry
data1<-within(redData1,{
  #Create event variable
  event <- as.numeric(sa_allc_ckdepi_renal1_v3y0 %in% c(1,2)) #############
  #Create a survival vectro
  #surv<-Surv(sa_allc_time_esrd_v3y0,event)
})



#Scale all Clinical parameters
redData2<-data1%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  mutate(log_acr = scale(log_acr))%>%
  mutate(adnine = scale(adnine))%>%
  dplyr::select(patientid,sa_allc_ckdepi_renal1_v3y0, sa_allc_ckdepi_time_renal3_v3y0,sex_v3y0,age,ha1c,map,bmi,log_acr,egfr,adnine,quartile)






#data2<-data1%>%
#  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0,map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,quartile,event)


#Plot survival Curve
sfit <- survfit(Surv(sa_allc_ckdepi_time_renal3_v3y0,  event)~quartile, data=redData2)
sfit
summary(sfit)
#plot(sfit)
library(survminer)
ggsurvplot(sfit)



#Plot cumulative incidence according to quartile

#plot(sfit, fun = "cumhaz")
ggsurvplot(sfit, fun = "cumhaz",pval = TRUE)



+ylab = c("Survival probability")




#Kaplan-Meier Curve for Red Cap data Survival
ggsurvplot(sfit,
           conf.int=TRUE, # add confidence intervals
           pval=TRUE, # show the p-value for the log-rank test
           risk.table=TRUE, # show a risk table below the plot
           legend.labs=c("Quartile 1", "Quartile 2"), # change group labels
           legend.title="Quartile",  # add legend title
           palette=c("dodgerblue4", "orchid2"), # change colors of the groups
           title="Kaplan-Meier Curve for Kidney Dialysis Survival", # add title to plot
           risk.table.height=.2)



#KM_fit <- survfit(Surv(time, status == 2) ~ sex, data = lung)


#Plot cumulative incidence according to quartile

#library(cmprsk)
ci_fit <- 
  cuminc(
    ftime = redData2$sa_allc_ckdepi_time_renal3_v3y0, 
    fstatus = redData2$event,
    group = redData2$quartile,
    cencode = 2
  )
ci_fit[["Tests"]]



