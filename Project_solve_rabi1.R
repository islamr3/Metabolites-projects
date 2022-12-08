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
#install.packages("ggplot2")
library('ggplot2')

#redcap_data_adnin


redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T) 

#Grouping dataset based on quartile method on adenine
#redcap_data_adnin <- redcap_data_adnin %>% mutate(quartile = ntile(adnine_v3y0, 2)) 


eskd <- c("sa_allc_ckdepi_renal1_v3y0","sa_allc_ckdepi_renal1_v5y1","sa_allc_ckdepi_renal1_v7y2","sa_allc_ckdepi_renal1_v9y3",
          "sa_allc_ckdepi_renal1_v11y4","sa_allc_ckdepi_renal1_v13y5","sa_allc_ckdepi_renal1_v15y6",
          "sa_allc_ckdepi_renal1_v17y7","sa_allc_ckdepi_renal1_v19y8","sa_allc_ckdepi_renal1_v21y9",
          "sa_allc_ckdepi_renal1_v23y10")



results <- function(dt,x){
  
  ahh <- as.list(match.call()[-1])
  

  name <- dt %>% 
    mutate(mea = (mean(urine_creatinine_v3y0,na.rm=T))*88.42*1000) %>% 
    mutate(acr=mea/urine_albumin_v3y0) %>% 
    mutate(log_acr = log(acr))%>%
    mutate(adnine = log(adnine_v3y0)) %>% 
    mutate(adnine_n=abs(scale(adnine_v3y0)))%>% 
    mutate(adnine_w=adnine_v3y0-mean(adnine_v3y0,na.rm = TRUE)/sd(adnine_v3y0,na.rm = TRUE))%>%
    select(patientid,ahh$x, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, 
           albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr,adnine,adnine_n,adnine_w) %>%
    #filter(albuminuria_v3y0==1)%>%   # Condition 1
    #filter(albuminuria_v3y0==2)%>%   # Condition 2
    #filter(albuminuria_v3y0!=3)%>%   # Condition 3
    #filter(albuminuria_v3y0!=4)%>%
    #filter(sex_v3y0==2)%>%#           #Condition 3
    #filter(sex_v3y0==1)%>%           # Condition 4
    filter(albuminuria_v3y0==1)%>%   # Condition 4
    #filter(sex_v3y0==1)%>%            # Condition 5
    #filter(albuminuria_v3y0==2)%>%    # Condition 5
    #filter(sex_v3y0==1)%>%            # Condition 6
    #filter(albuminuria_v3y0!=3)%>%   # Condition 6
    #filter(albuminuria_v3y0!=4)%>%   # Condition 6
    #filter(sex_v3y0==2)%>%            #Female
    #filter(albuminuria_v3y0==1)%>%    #Noemal
    #filter(albuminuria_v3y0!=3)%>% 
    #filter(albuminuria_v3y0!=4)%>%    #Micro
    filter(ahh$x!=9) %>%              
    na.omit()
  
  
event <- as.numeric(name[,ahh$x] %in% c(1,2))
surv<-Surv(name$sa_allc_time_esrd_v3y0,event)

redcap_data3 <- cbind(name,event,surv)
  
print(dim(name))

#Adenine+Clinical parameters
logity1 <- glm(event ~ age_integer_v3y0+hemoglobin_a1c_v3y0+map_v3y0+bmi_v3y0+egfr_cric_v3y0+log_acr+adnine+adnine_n, data = redcap_data3, family = binomial)
print(lroc(logity1, graph = T)$auc)
lr1=lroc(logity1, graph = T)

AUCdata1=lr1$diagnostic.table

AUCdata1 <- data.frame(AUCdata1)
#aucFrame=data.frame(AUCdata)
#aucFrame%>%
#  ggplot(aes(x=X1.Specificity, y=Sensitivity),color=metric)+
#  geom_line()



#Clinical parameters
logity2 <- glm(event ~ age_integer_v3y0+hemoglobin_a1c_v3y0+map_v3y0+bmi_v3y0+egfr_cric_v3y0+log_acr, data = redcap_data3, family = binomial)
print(lroc(logity2, graph = T)$auc)


lr2=lroc(logity2, graph = T)


AUCdata2=lr2$diagnostic.table

AUCdata2<- data.frame(AUCdata2)

datafram2 <- data.frame(AUCdata1,AUCdata2)

print(redcap_data3%>%count(event==1))
print(anova(logity1, logity2, test = "LRT"))



#class(datafram2)
#names(datafram2)
#colnames(AUCdata1) <- c("1-spec1","sen1") 
#colnames(AUCdata2) <- c("1-spec2","sen2") 



### Cbind function using for column bind
#rocy <- cbind(AUCdata1$Sensitivity, AUCdata1$X1.Specificity, AUCdata2$Sensitivity, AUCdata2$X1.Specificity)
#class(rocy)

#rocy <- as.data.frame(rocy)

### Roc Curve

#names(rocy)




#so<-c("Cp+Ad"="red", "Cp" = "#377eb8")
so<-c("Cp+Ad"="blue2", "Cp" = "black")


ggplot(datafram2) + 
  geom_line(aes(y = Sensitivity, x = X1.Specificity, color="Cp+Ad")) +
  geom_line(aes(y = Sensitivity.1, x = X1.Specificity.1,color="Cp"))+
  theme(plot.title = element_text(hjust = 0.5))+
  #labs(title = "Compared Two Overlapping Roc Curve",
  labs(title = "Compared Roc Curve",
       x="1-Specificity",
       y="Sensitivity",
       color="Legend")+
  scale_color_manual(values = so)


#ggplotly(p2)


}



par(mfrow=c(2,5))


a <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v5y1")

b <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v7y2")

c <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v9y3")

d <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v11y4")
e <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v13y5")

f <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v15y6")
g <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v17y7")

h <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v19y8")
i <- results(redcap_data_adnin,"sa_allc_ckdepi_renal1_v21y9")
v <- ggarrange(a, b,c,d,e,f,g,h,i, nrow = 2, ncol = 5)
v
#j <- sojib(redcap_data_adnin,"sa_allc_ckdepi_renal1_v23y10")


#v <- ggarrange(a, b,c,d,e,f,g,h, nrow = 2, ncol = 5)

#Baseline data
#base <- sojib(redcap_data_adnin,"sa_allc_ckdepi_renal1_v3y0")

#v <- ggarrange(a, b,c,d,e,f,g,h,i,j, nrow = 2, ncol = 5)
#v
