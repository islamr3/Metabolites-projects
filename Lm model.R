#####

#setwd(choose.dir())

#####

##### Required Library ##### 

library(dplyr)
library(ggplot2)

#### Read the dataset #####

redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)%>%
  subset()


###################################Male eGFR selection for Normal
EGFR<-redcapegfr%>%
  filter(sex_v3y0==1)%>%
  filter(albuminuria_v3y0==1)%>%
  dplyr::select(egfr_cric_v3y0,egfr_cric_v5y1,egfr_cric_v7y2,egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6,egfr_cric_v17y7,
                egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10,egfr_cric_v25y11,egfr_cric_v27y12,adnine_v3y0)

dim(EGFR)



##### create the required dataset ####

lmd <- redcap_data_adnin %>% 
  select(egfr_cric_v3y0, egfr_cric_v5y1, egfr_cric_v7y2, egfr_cric_v9y3,egfr_cric_v11y4,egfr_cric_v13y5,egfr_cric_v15y6, egfr_cric_v17y7, egfr_cric_v19y8,egfr_cric_v21y9,egfr_cric_v23y10)#%>% na.omit()
#### print the first six value of the new dataset #####
head(lmd)

##### Using lope for the linear model coefficient #####

slope <- numeric()
intercept <- numeric()
for (i in 1:nrow(lmd)) {
  #time <- 1:10
  egfr <- as.numeric(lmd[i,])
  egfr
  #Ramove nan from egfr
  g<-egfr[!is.na(egfr)]
  #define time variable base on egfr length
  time <- 1:length(g)
  model <- lm(time ~ g)
  intercept[i] <- model$coefficients[1]
  slope[i] <- model$coefficients[2]
}
## Print the intercept value 
intercept
### print the slope value
slope
#### Creating the new data set for intercept and slope ####
dt <- cbind(intercept,slope)
dt <- as.data.frame(dt)

head(dt)
#



