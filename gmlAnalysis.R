rm(list=ls())
#Random number seed
set.seed(10079)
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



##Load redcap data
redcap_data_adnin <-  read.csv("RedCapUp2.csv",stringsAsFactors = T)

#Select the normal Albumine for male, Female and male+Female
redData<-redcap_data_adnin%>%
  filter(sex_v3y0==2)%>%
  #filter(albuminuria_v3y0==1)%>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_time_esrd_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0,adnine_v3y0)%>%
  na.omit()

#Estimate log ACR and adnine ACR in units of mg/mmol
redDataR<-redData%>%
  mutate(acr=((urine_albumin_v3y0/urine_creatinine_v3y0)*1000)/88.42) %>%
  #mutate(acr=urine_albumin_v3y0/acr)%>%
  mutate(log_acr = log(acr))%>%
  mutate(adnine = log(adnine_v3y0))%>%
  na.omit()

#Selection of normal albuminuria based on acr group
data<-redDataR%>%
  filter(acr<=3)%>%
  na.omit()


#Create event sumarry
data<-within(data,{
  #Create event variable
  event <- as.numeric(sa_allc_ckdepi_renal1_v3y0 %in% c(1,2)) #############
  time<-sa_allc_time_esrd_v3y0
  #Create a survival vectro
  surv<-Surv(time,event)
})

#Normalize covariates data
dataN<-data%>%
  mutate(age=scale(age_integer_v3y0))%>%
  mutate(ha1c=scale(hemoglobin_a1c_v3y0))%>%
  mutate(map=scale(map_v3y0))%>%
  mutate(bmi=scale(bmi_v3y0))%>%
  mutate(egfr=scale(egfr_cric_v3y0))%>%
  mutate(adnine=scale(adnine))%>%
  na.omit()


#define covariates  
covariates <- c(names(dataN[20:24]),names(dataN[15]),names(dataN[16]))
#covariates <- names(dataN[15:21])

#Load Normalize data into rawdata
rawdata<-dataN%>%
  na.omit()

# treat missing below detection
rawdata[is.na(rawdata)] <- 0.

#treat zeros with half of the minimun
rawdata[covariates] <- lapply(rawdata[covariates], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/2))


univg_formula<- sapply(covariates, function(x) as.formula(paste("event ~", x)))

univg_models<- lapply(univg_formula, function(x) {glm(x, data = rawdata, family = binomial)})


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = rawdata)})

# glm model Extract data 
univg_results <- lapply(univg_models,
                       function(x){ 
                         m<-x
                         x <- summary(x)
                         #print(x)
                         #p.value<-signif(x$coef[8], digits=2)
                         value<-anova(m, test="Chisq")
                         p.value<-signif(value[2,5],digits=2)
                         #wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[2], digits=2);#coeficient beta
                         HR<-signif(exp(beta),digits=2);#exp(beta)
                         cofi<-plogis(confint(m))
                         HR.confint.lower<-signif(cofi[2],digits = 2)
                         HR.confint.upper<-signif(cofi[4],digits=2)
                         #res<-c(beta, HR,HR.confint.lower)
                         HR <- paste0(HR, " (", 
                                     HR.confint.lower, "-", HR.confint.upper, ")")  
                         resg<-c(beta,HR,p.value)
                         names(resg)<-c("beta", "HR (95% CI for HR)","p.value")                        
                         return(resg)
                         return(exp(cbind(coef(x),confint(x))))
                       })
resg <- t(as.data.frame(univg_results, check.names = FALSE))


#Cox model
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = rawdata)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
#Two model comparision
as.data.frame(resg)
as.data.frame(res)


#Next, we compare ESRD survival curves for 2 groups divided by the 50th percentile.

#Visualization - compare pathways dysregulation stratified by 50% percentile

#quantize path data
nstrata = 2
quant_data <- rawdata%>%
  na.omit()


quant_data[covariates] <- apply(rawdata[covariates], MARGIN = 2, 
                                function(x){
                                  cut(x,breaks = quantile(x,probs = seq(0,1,by = 1/nstrata),na.rm = TRUE))
                                }
)

#get survival formulas
uvar_formulas <- sapply(covariates,
                        function(x){
                          as.formula(paste("Surv(time, event)",x, sep=" ~ "))
                        }
)
#fit models
uvar_models <- lapply( uvar_formulas, function(x){surv_fit(x, data = quant_data)})

#plot models
ggsurvplot(uvar_models, conf.int = TRUE, pval = TRUE)

## 
## attr(,"class")
## [1] "list"            "ggsurvplot_list"
#$Lysine

