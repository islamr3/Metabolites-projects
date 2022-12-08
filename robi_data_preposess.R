 
#install.packages("readxl")
library("readxl")

#install.packages("dplyr")
library('dplyr')


CRIC_data_ccombine = read_excel("/Users/mostsheuliakter/Desktop/code/CRIC_AA_Combined_final_Creatinine_Normalized_adenine at nM vs mM.xlsx")%>%
  dplyr::select(patientid, 'Adenine (nM/uM)')
  

redcap_data = read.csv("/Users/mostsheuliakter/Desktop/code/RedCapUp2.csv") %>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_ckdepi_time_renal3_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,urine_creatinine_v3y0,urine_albumin_v3y0) %>%
  mutate(mean_acr=(mean(urine_creatinine_v3y0, na.rm = TRUE))*88.42*1000) %>%
  mutate(acr=mean_acr/urine_albumin_v3y0)%>%
  mutate(log_acr = log(acr))%>%
  dplyr::select(patientid, sa_allc_ckdepi_renal1_v3y0, sa_allc_ckdepi_time_renal3_v3y0, sex_v3y0,age_integer_v3y0, albuminuria_v3y0, hemoglobin_a1c_v3y0, map_v3y0,bmi_v3y0,egfr_cric_v3y0,log_acr)
#redcap_data 
  
redcap_data_adnin <- merge(redcap_data,CRIC_data_ccombine,by="patientid")%>%  # Merge Data in Here
  rename(c(pid=patientid, time=sa_allc_ckdepi_time_renal3_v3y0, sex=sex_v3y0, age=age_integer_v3y0, 
           albuminuria=albuminuria_v3y0, hemoglobin_a1c=hemoglobin_a1c_v3y0, map=map_v3y0, bmi=bmi_v3y0, egfr_cric=egfr_cric_v3y0, adenin= 'Adenine (nM/uM)'))
redcap_data_adnin
dim(redcap_data_adnin)


