rm(list=ls())
library(dplyr)
library(ggplot2)
library(rstatix)  
library(ggpubr)  # For easy data-visualization
library(outliers)

getwd()
#### Read the dataset #####

redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)

EGFR<-redcap_data_adnin%>%
  #filter(sex_v3y0==1)%>%
  #filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12)

slope <- numeric()
k <- 1
ra<-NULL;
for (i in 1:nrow(EGFR)) {
  #time <- 1:10
  egfr <- as.numeric(EGFR[i,])
  IDX <- is.na(egfr)
  IDX
  egfr <- na.omit(egfr)
  t <- 1:13
  t1=t[IDX!=1]
  if(length(egfr)>=3)
  {
    model <- lm(t1 ~ egfr)
    slope[k] <- model$coefficients[2]
    k = k+1
  }
  else
  {
    ra <- c(ra,i)
    #ra=i;
  }
}
slope

#Load adnine data
adnine<-redcap_data_adnin%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  filter(!row_number() %in% ra)%>% # Remove adnine data which are less than eGFR less than 3
  dplyr::select(adnine_v3y0)

#Track the nan from adnine
IDX1 <- is.na(adnine)

#Remove Nan based on IDX
adnine1<-adnine[IDX1!=1]
adnine1<-log(adnine1)
slope1<-slope[IDX1!=1]

#z-score normalization of adnine and slope1
x_stand2a <- scale(adnine1)

#outliner detection and remove
x_stand2b <- scale(slope1)

#Making dataframe with Adnine and slope
dt1<-cbind(x_stand2a,x_stand2b)
dt1<-as.data.frame(dt1)
names(dt1)
names(dt1)<-c("Adnine","Slope")
dt2=rm.outlier(dt1,fill=FALSE, opposite=FALSE)

model <- lm(dt2$Slope~ dt2$Adnine)
summary(model)

