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

redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)

redData<-redcap_data_adnin%>%
  #filter(sex_v3y0==1)%>%
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

#Create tertile
data2<-data1%>%
  mutate(quartile=ntile(adnine,3))


#Selecting lower quartile
lowqdata3<-data2%>%
  filter(quartile==1)%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  na.omit()

#Selecting lower quartile
middledata3<-data2%>%
  filter(quartile==2)%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  na.omit()


#Selecting lower quartile
upperdata3<-data2%>%
  filter(quartile==3)%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  na.omit()

model1 <- coxph( Surv(sa_allc_time_esrd_v3y0, event) ~ age+ha1c+map+bmi+egfr+log_acr+adnine,
                 data = lowqdata3 )
summary(model1)


model2 <- coxph( Surv(sa_allc_time_esrd_v3y0, event) ~ age+ha1c+map+bmi+egfr+log_acr+adnine,
                 data = middledata3 )
summary(model2)


model3 <- coxph( Surv(sa_allc_time_esrd_v3y0, event) ~ age+ha1c+map+bmi+egfr+log_acr+adnine,
                 data = upperdata3 )

summary(model3)

#ggforest(model2)

#Forest plot based on ggplot2
forestData=function(model1){
  m1<-summary(model1)
  m2<-m1[["conf.int"]]
  m3<-melt(rownames(m2))
  names(m3)<-"covariate"
  
  rownames(m2)<-NULL
  Fdata1<-cbind(m3,m2)
  
  return(Fdata1)
}

Q=rep(NA,dim(mdata1)[1])

ter1<-c("Q1","Q1","Q1")
ter2<-c("Q2","Q2","Q3")
ter3<-c("Q3","Q3","Q4")

map_to_tertile = function(x){
  
  quartiles = quantile(x, na.rm = T, probs = seq(0, 1, by = 1/3))
  
  vect_return = rep(NA, length(x)) %>%
    (function(vect) ifelse(x >= quartiles[1] & x < quartiles[2], 'Q1', vect)) %>%
    (function(vect) ifelse(x >= quartiles[2] & x < quartiles[3], 'Q2', vect)) %>%
    (function(vect) ifelse(x >= quartiles[3], 'Q3', vect)) %>% 
    as.factor()
  
  return(vect_return)
}
  
mdata1=forestData(model1)
mdata2=forestData(model2)
mdata3=forestData(model3)

Forestmodels<-rbind(mdata1,mdata2,mdata3)

m1<-summary(model1)
m2<-m1[["conf.int"]]
m3<-melt(rownames(m2))
names(m3)<-"covariate"

rownames(m2)<-NULL
Fdata1<-cbind(m3,m2)

library(forestplot)
library(magrittr)
library(sjPlot)

plot_models(model1, model2,model3, std.est = "std2", transform=NULL,show.values = TRUE, dot.size=1, vline.color = "black", m.labels = c("Lower", "Middle","Upper"), ci.lvl = 0.95)+ylab("Odds Ratio")+scale_colour_manual(
  values = c("black","blue","#009999"))



gg <- plot_models(m19, type="std2", show.p=TRUE, ci.lvl=.95, dot.size=5, 
                  line.size=2, vline.color = "black", axis.lim=c(.001, 2.25)) + 
  theme_bw(base_size=24) +
  theme(panel.grid.major.y=element_blank()) +
  ggtitle("") + 
  ylab("Odds Ratio")

