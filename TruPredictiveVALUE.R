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

#Were you able to generate the data for positive and negative 
#predictive values for urine adenine to predict ESKD outcomes 
#in the NA group at 3 years.

redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)

#Select the normal Albumine for male, Female and male+Female
redData<-redcap_data_adnin%>%
  filter(sex_v3y0==1)%>%
  #filter(albuminuria_v3y0==1)%>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v9y3, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0)%>%
  na.omit()

#Estimate log ACR and adnine ACR in units of mg/mmol
redDataR<-redData%>%
  mutate(acr=((urine_albumin_v3y0/urine_creatinine_v3y0)*1000)/88.42) %>%
  #mutate(acr=urine_albumin_v3y0/acr)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))

#Selection of normal albuminuria based on acr group
redData1<-redDataR%>%
  filter(acr<=3)%>%
  na.omit()


#Create event sumarry
data1<-within(redData1,{
  #Create event variable
  event <- as.numeric(sa_allc_ckdepi_renal1_v3y0 %in% c(1,2)) #############
  #Create a survival vectro
  surv<-Surv(sa_allc_time_esrd_v9y3,event)
})


#Normalize data
data2<-data1%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  na.omit()

#Check data with lebel PV=0.096
ggplot(data2, aes(x = event, fill = event))+geom_bar()
length(which(data2$event==1))
length(which(data2$event==0))
length(which(data2$sa_allc_ckdepi_renal1_v3y0==9))


#fit logistic regression model
#model1 <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr, family="binomial", data=data2)

#http://www.sthda.com/english/wiki/cox-proportional-hazards-model
#cox <- coxph(surv ~ age + ha1c + map + bmi + age + egfr+adnine, data = data2)
#ggsurvplot(survfit(cox), data=data,ggtheme = theme_minimal())
#library(regclass)
#confusion_matrix(model)
#data2$group <- data2$adnine>0
#cox <- coxph(surv ~ age + ha1c + map + bmi + age + egfr+group, data = data2)
ntile(data1$adenine,3)

data2<-data1%>%
  ntile(data1$adenine,3)

#Cutoff value
quantile(data1$adnine, c(1/3, 2/3))
#adnine>10.10 if more than 


LowerQ<-data2%>%
  filter(quartile==1)


# Survival curves
#fit <- survfit(cox, data = data2)
#ggsurvplot(fit, conf.int = TRUE, ggtheme = theme_minimal())

model1 <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine, family="binomial", data=data2)  
library(regclass)
library(caret)
library(InformationValue)
library(ISLR)

#use model to predict probability of default
pred<- predict(model1, data2, type="response")

optimal <- optimalCutoff(data2$event, pred)[1]

#create confusion matrix
r=confusionMatrix(data2$event, pred,threshold = round(optimal, 2))

TP=r[1,1]
FP=r[1,2]
FN=r[2,1]
TN=r[2,2]

PPV1=TP/(TP+FP)
PPV1
NPV1=TN/(TN+FN)
NPV1


#Traing and testing model results
library(caret)
#PPV1<- numeric()
PPV1 <- c()
NPV1<-c()
index <- createDataPartition(data2$event, p = 0.7, list = FALSE,times = 5)
for(i in 1:5){
  testData <- data2[index[,i], ]
  trainData  <- data2[-index[,i], ]
  #Use test and train data partitions however you desire...
  model <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr, family="binomial", data=testData)
  pred=predicted <- predict(model, testData, type="response")
  optimal <- optimalCutoff(testData$event, pred)[1]
  r=confusionMatrix(testData$event,pred,threshold = optimal)
  #print(r)
  
  TP=r[1,1]
  FP=r[1,2]
  FN=r[2,1]
  TN=r[2,2]
  PPV1[i]=TP/(TP+FP)
  #print(PPV1)
  
  NPV1[i]=TN/(TN+FN)
  #print(NPV1)
}

print(mean(PPV1))
print(sd(PPV1))
print(sd(PPV1)/mean(PPV1))

print(mean(NPV1))
print(sd(NPV1))
print(sd(NPV1)/mean(NPV1))












#Randomly shuffle the data
yourData<-data2[sample(nrow(data2)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==2,arr.ind=TRUE)
  testData <- yourData[testIndexes, ]
  trainData <- yourData[-testIndexes, ]
  #Use test and train data partitions however you desire...
  model <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine, family="binomial", data=trainData)
  pred=predicted <- predict(model, testData, type="response")
  optimal <- optimalCutoff(data2$event, predicted)[1]
  r=confusionMatrix(testData$event, pred,threshold = optimal)
  r
}


library(caret)

# define training control
train_control <- trainControl(method = "cv", number = 2)

# train the model on training set
model <- train(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine,
               data = data2,
               trControl = train_control,
               method = "glm",
               family=binomial())

# print cv scores
summary(model)


library(caret)
index <- createDataPartition(data2$event, p = 0.7, list = FALSE,times = 5)
for(i in 1:5){
  testData <- data2[index[,i], ]
  trainData  <- data2[-index[,i], ]
  #Use test and train data partitions however you desire...
  model <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine, family="binomial", data=testData)
  pred=predicted <- predict(model, testData, type="response")
  optimal <- optimalCutoff(testData$event, pred)[1]
  r=confusionMatrix(testData$event,pred,threshold = optimal)
  print(r)
}




for(i in 1:2){
  set.seed(42)
  #Segement your data by fold using the which() function 
  index <- createDataPartition(data2$event, p = 0.7, list = FALSE,times = 3)
  testData <- data2[index, ]
  trainData  <- data2[-index, ]
  #Use test and train data partitions however you desire...
  model <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine, family="binomial", data=trainData)
  pred=predicted <- predict(model, testData, type="response")
  optimal <- optimalCutoff(data2$event, predicted)[1]
  r=confusionMatrix(testData$event, pred,threshold = optimal)
  r
}


index <- createDataPartition(data2$event, p = 0.7, list = FALSE)
train_data <- data2[index, ]
test_data  <- data2[-index, ]



#use model to predict probability of default
pred=predicted <- predict(model, data2, type="response")

optimal <- optimalCutoff(data2$event, predicted)[1]

#create confusion matrix
r=confusionMatrix(data2$event, pred,threshold = round(optimal, 2))

TP=r[1,1]
FP=r[1,2]
FN=r[2,1]
TN=r[2,2]

PPV1=TP/(TP+FP)
PPV1
NPV1=TN/(TN+FN)
NPV1

#Check train and test label with data
ggplot(train_data, aes(x = event, fill = event))+geom_bar()
ggplot(test_data, aes(x = event, fill = event))+geom_bar()


#fit logistic regression model
model <- glm(event ~ age+ha1c+map+bmi+egfr+log_acr+adnine, family="binomial", data=train_data)
#use model to predict probability of default
pred=predicted <- predict(model, test_data, type="response")


pred[pred >.5] = "Up"

table(pred,event)

confusionMatrix(test_data$event, pred)

confusionMatrix() 



#Create tertile
data2<-data1%>%
  mutate(quartile=ntile(adnine,2))


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

#Selecting upper quartile
upperdata3<-data2%>%
  filter(quartile==2)%>%
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
                 data = upperdata3 )
summary(model2)

#ggforest(model2)


library(forestplot)
library(magrittr)
library(sjPlot)


plot_models(model1, model2, std.est = "std2", transform=NULL,show.values = TRUE, vline.color = "black", m.labels = c("Lower", "Upper"), ci.lvl = 0.95)+ylab("Odds Ratio")+scale_colour_manual(
  values = c("black","blue"))


p3<-plot_model(model1,std.est = "std2", transform=NULL,show.values = TRUE, vline.color = "black", m.labels = c("Lower", "Upper"), ci.lvl = 0.95)+ylab("Odds Ratio")


p4<-plot_model(model2, std.est = "std2", transform=NULL,show.values = TRUE, vline.color = "black", m.labels = c("Lower", "Upper"), ci.lvl = 0.95)+ylab("Odds Ratio")

v<-ggarrange(p3, p4, ncol=1, nrow=1) 

#Use ggplot in here

library(broom)
library(forestmangr)
model_output <- tidy(model1)
out_conf <- tidy(model1, conf.int = TRUE)
lm_model_out <- round_df(out_conf, digits=2)
lm_model_out1 <- lm_model_out[-1,]

p1<-ggplot(lm_model_out, aes(x=reorder(term, estimate), y=estimate)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), 
                width = 0.2,size  = 1,
                position = "dodge", color="turquoise4") +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  geom_point() + coord_flip() 

model_output <- tidy(model2)
out_conf1 <- tidy(model2, conf.int = TRUE)
lm_model_out2 <- round_df(out_conf1, digits=2)

p2<-ggplot(lm_model_out2, aes(x=reorder(term, estimate), y=estimate)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), 
                width = 0.2,size  = 1,
                position = "dodge", color="turquoise4") +
  geom_hline(yintercept = 0, color = "red", size = 1) +
  geom_point() + coord_flip() 

v<-ggarrange(p1, p2, ncol=1, nrow=1) 


coff=model1$coefficients
coff1=melt(coff)
Col1 <- rownames(coff1)
Col2 <- coff1[,1]
dt=cbind(Col1,Col2)
#############################################################




#No need
plot_models(model1, model2, std.est = "std2", transform=NULL,show.values = TRUE, vline.color = "black", m.labels = c("Lower", "Upper"), ci.lvl = 0.95)+ylab("Odds Ratio")+scale_color_sjplot()

l=c("red","blue")

?plot_models()

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
