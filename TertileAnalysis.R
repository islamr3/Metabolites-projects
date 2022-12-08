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



model3 <- glm( event ~ age+ha1c+map+bmi+egfr+log_acr+adnine,
               data = lowqdata3 )
summary(model1)


model4 <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine,
              data = upperdata3 )
summary(model2)

plot_models(model3, model4,vline.color = "blue", std.est = "std2", m.labels = c("Lower", "Upper"))





Lower<-lowqdata3%>%
  dplyr::select(age,ha1c,map,bmi,egfr,log_acr,adnine,)




forestplot(
  df = Lower,
  estimate = beta,
  logodds = FALSE,
  colour = trait,
  xlab = "1-SD increment in cardiometabolic trait
  per 1-SD increment in biomarker concentration"
)

# Plot KM by Tertile
fit2 = survfit(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1)~lactate_tertile, data = model_data)

# fit survival model with weights
coxmod2 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile,
                weights=w2,
                data=model_data)






qdf<-df%>%
  mutate(quartile=ntile(Adnine,2))  #creating quartile catagorical variable


EGFR<-redcap_data_adnin%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12)

redcap_data_adnin1 = read.csv("/Users/rabiul/fwdcodeusedinanalysis/RedCapUp2.csv") %>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0) %>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>% #ArrDelay/max(ArrDelay)
  mutate(adnine_n=abs(scale(adnine_v3y0)))%>% 
  mutate(adnine_w=adnine_v3y0-mean(adnine_v3y0,na.rm = TRUE)/sd(adnine_v3y0,na.rm = TRUE))%>%#-mean(adnine_v3y0)/sd(adnine_v3y0,na.rm = FALSE))%>%
  dplyr::select(patientid,sa_allc_ckdepi_renal1_v5y1, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_v3y0,adnine_n,adnine_w)%>%
  #rename(ESKD='sa_allc_ckdepi_renal1_v9y3')%>%
  #rename_all(function(x) gsub('_v9y3', '', x))%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  na.omit()


#CRIC_data_ccombine = read_excel("/Users/rabiul/fwdcodeusedinanalysis/CRIC_data_ccombine1.xlsx")
#dplyr::select(patientid, 'Adenine (nM/uM)') #%>%
#rename(adnine='Adenine (nM/uM)')

#CRIC_data_ccombine = read_excel("/Users/rabiul/fwdcodeusedinanalysis/CRIC_AA_Combined_final_Creatinine_Normalized_adenine at nM vs mM.xlsx")%>%
#  dplyr::select(patientid, 'Adenine (nM/uM)') #%>%
#rename(adnine='Adenine (nM/uM)')

redcap_data_adnin1 = read.csv("/Users/rabiul/fwdcodeusedinanalysis/RedCapUp2.csv") %>%
  #dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0) %>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>% #ArrDelay/max(ArrDelay)
  mutate(adnine_n=abs(scale(adnine_v3y0)))%>% 
  mutate(adnine_w=adnine_v3y0-mean(adnine_v3y0,na.rm = TRUE)/sd(adnine_v3y0,na.rm = TRUE))%>%#-mean(adnine_v3y0)/sd(adnine_v3y0,na.rm = FALSE))%>%
  dplyr::select(patientid,sa_allc_ckdepi_renal1_v5y1, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_v3y0,adnine_n,adnine_w)%>%
  #rename(ESKD='sa_allc_ckdepi_renal1_v9y3')%>%
  rename_all(function(x) gsub('_v9y3', '', x))%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==4)%>%
  #colnames(redcap_data_adnin1)[2] <- "ESKD" %>%
  #colnames(redcap_data_adnin1)[colnames(redcap_data_adnin1) == 'sa_allc_ckdepi_renal1_v9y3'] <- 'ESKD'
  #filter(albuminuria_v3y0!=3)%>%
  #filter(albuminuria_v3y0!=4)%>%
  #filter(sa_allc_ckdepi_renal1 !=9)%>%   ##############
na.omit()

redcap_data_adnin2 = read.csv("/Users/rabiul/fwdcodeusedinanalysis/RedCapUp2.csv") %>%
  #dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0) %>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>% #ArrDelay/max(ArrDelay)
  mutate(adnine_n=abs(scale(adnine_v3y0)))%>% 
  mutate(adnine_w=adnine_v3y0-mean(adnine_v3y0,na.rm = TRUE)/sd(adnine_v3y0,na.rm = TRUE))%>%#-mean(adnine_v3y0)/sd(adnine_v3y0,na.rm = FALSE))%>%
  dplyr::select(patientid,sa_allc_ckdepi_renal1_v5y1, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_v3y0,adnine_n,adnine_w)%>%
  #rename(sa_allc_ckdepi_renal1_v23y10=ESKD)%>%
  rename_all(function(x) gsub('_v13y5', '', x))%>%
  filter(albuminuria_v3y0!=3)%>%
  filter(albuminuria_v3y0!=4)%>%
  filter(sa_allc_ckdepi_renal1 !=9)%>%   ##############
na.omit()


redcap_data_adnin3 = read.csv("/Users/rabiul/fwdcodeusedinanalysis/RedCapUp2.csv") %>%
  #dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0) %>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>% #ArrDelay/max(ArrDelay)
  mutate(adnine_n=abs(scale(adnine_v3y0)))%>% 
  mutate(adnine_w=adnine_v3y0-mean(adnine_v3y0,na.rm = TRUE)/sd(adnine_v3y0,na.rm = TRUE))%>%#-mean(adnine_v3y0)/sd(adnine_v3y0,na.rm = FALSE))%>%
  dplyr::select(patientid,sa_allc_ckdepi_renal1_v23y10, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_v3y0,adnine_n,adnine_w)%>%
  #rename(sa_allc_ckdepi_renal1_v23y10=ESKD)%>%
  rename_all(function(x) gsub('_v23y10', '', x))%>%
  filter(albuminuria_v3y0!=3)%>%
  filter(albuminuria_v3y0!=4)%>%
  filter(sa_allc_ckdepi_renal1 !=9)%>%   ##############
na.omit()



redcap_data_adnin<-rbind(redcap_data_adnin1,redcap_data_adnin2,redcap_data_adnin3,all = TRUE)

redcap_data_adnin<-distinct(redcap_data_adnin)


#Create event sumarry
redcap_data3<-within(redcap_data_adnin,{
  #Create event variable
  event <- as.numeric(sa_allc_ckdepi_renal1 %in% c(1,2)) #############
  #Create a survival vectro
  surv<-Surv(sa_allc_time_esrd_v3y0,event)
})

redcap_data3 <- data.frame(redcap_data3)

#Adenine+Clinical parameters
logity1 <- glm(event ~ age_integer_v3y0+hemoglobin_a1c_v3y0+map_v3y0+bmi_v3y0+egfr_cric_v3y0+log_acr+adnine+adnine_v3y0+adnine_n+adnine_w, data = redcap_data3, family = binomial)
logity1 <- glm(event ~ age_integer_v3y0+egfr_cric_v3y0, data = redcap_data3, family = binomial)
lroc(logity1, graph = T)$auc

#Clinical parameters
logity2 <- glm(event ~ age_integer_v3y0+hemoglobin_a1c_v3y0+map_v3y0+bmi_v3y0+egfr_cric_v3y0+log_acr, data = redcap_data3, family = binomial)
lroc(logity2, graph = T)$auc


########################################3

lr1=lroc(logity1, graph = T)
AUCdata1=lr1$diagnostic.table

lr2=lroc(logity2, graph = T)
AUCdata2=lr2$diagnostic.table


AUCdata1 <- data.frame(AUCdata1)
AUCdata2<- data.frame(AUCdata2)


rocy <- cbind(AUCdata1$Sensitivity, AUCdata1$X1.Specificity, AUCdata2$Sensitivity, AUCdata2$X1.Specificity)
rocy <- as.data.frame(rocy)

sojib<-c("Ad and CP"="blue", "CP" = "black")
p2 <- ggplot(rocy) + 
  geom_line(aes(y = V1, x = V2, color="Ad and CP")) +
  geom_line(aes(y = V3, x = V4,color="CP"))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="1-Specificity",
       y="Sensitivity",
       color="Legend")+
  scale_color_manual(values = sojib)
ggplotly(p2)
ggsave("9.png", width = 4, height = 3, units = "in", dpi = 300)


redcap_data3%>% count(event==1)
anova(logity1, logity2, test = "LRT")

log_adenine_percentiles = test %>%
  quantile(na.rm = T) %>%
  data.frame() %>%
  rename(Value = 1) %>%
  tibble::rownames_to_column('Percentile') %>%
  gt()


#Grouping dataset based on quartile method on adenine
temp <- redcap_data_adnin1 %>% mutate(quartile = ntile(adnine_v3y0, 2))
