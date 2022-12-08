
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
library("read_excel")

##source.file <-"survival_model.r"
##project.id <- "Sharma_653"
##ource_info <- create_source_file_dir(source.description="")

# Program body here

##setwd(paste0(getProjectPath(), '/Programs'))

##################################################################
# Load Data
# 
##################################################################

#RedCap Dataset
print(getwd())
setwd("/Users/rabiul/fwdcodeusedinanalysis")
#redcap_data = read.csv("/Users/rabiul/fwdcodeusedinanalysis/NIHDP3CRIC_DATA_2021-06-25_1359.csv")
#redcap_data = read.csv(paste0(dataDir(),"/NIHDP3CRIC_DATA_2021-06-25_1359.csv")) %>%

redcap_data = read.csv("/Users/rabiul/fwdcodeusedinanalysis/NIHDP3CRIC_DATA_2021-06-25_1359.csv") %>%
  dplyr::select(patientid, sa_allc_ckdepi_time_renal1_v5y1, sa_allc_ckdepi_renal1_v5y1) %>%
  rename_all(function(x) gsub('_v5y1', '', x))

##redcap_data = read.csv("NIHDP3CRIC_DATA_2021-06-25_1359.csv") %>%
##  dplyr::select(patientid, sa_allc_ckdepi_time_renal1_v5y1, sa_allc_ckdepi_renal1_v5y1) %>%
##  rename_all(function(x) gsub('_v5y1', '', x))

#Note: Added variable (09/03/21)
##filename = grep('CRIC_AA', list.files(dataDir()), v=T)

##new_variable = Read(filename, read.fcn = xlsx::read.xlsx, sheetIndex = 1) %>%
##  dplyr::select(1,6) %>%
##  dplyr::rename('adenine' = 2) %>%
##  filter(!duplicated(patientid)) %>%
##  mutate(adenine = as.numeric(adenine))

#Previously used datasets
all_data = Load.branch('all_data.Rdata')
lactate_baseline_long = all_data$three_yr
lactate_baseline = lactate_baseline_long %>%
  filter(year == 0) %>%
  merge(redcap_data, by = 'patientid', all.x  = T)
  #merge(redcap_data, by = 'patientid', all.x  = T) %>%
  #merge(new_variable, by = 'patientid', all.x = T)

##################################################################
# Create data for survival model  
# 
##################################################################

covariates_used = list(c('lactate_tertile','age','bmi','sex','albumin_creatinine_ratio',
                         'hemoglobin_a1c','map','patientid','race_ethnicity_cat',
                         'metformin', 'egfr_ckd_epi'),
                       c('lactate_tertile','age','bmi','sex','albumin_creatinine_ratio',
                         'hemoglobin_a1c','map','patientid','race_ethnicity_cat',
                         'metformin', 'egfr_ckd_epi', 'adenine'),
                       c('lactate_tertile','age','bmi','sex', 'hemoglobin_a1c',
                         'map','patientid','race_ethnicity_cat','metformin'),
                       c('lactate_tertile','age','bmi','sex', 'hemoglobin_a1c',
                         'egfr_ckd_epi', 'map','patientid','race_ethnicity_cat',
                         'metformin', 'albumin_creatinine_ratio'),
                       c('lactate_tertile','age','bmi','sex', 'hemoglobin_a1c',
                         'albumin_creatinine_ratio', 'map','patientid','race_ethnicity_cat',
                         'metformin')) %>%
  set_names(paste0('cbpsdata', 2:6))


propensity_weight_list = covariates_used %>%
  map(function(x){
    
    formula_hold = paste0('lactate_tertile~',
                          paste(setdiff(x, c('lactate_tertile', 'patientid')), collapse = '+'))
    
    data_hold = lactate_baseline %>%
      select(c(x, patientid)) %>%
      na.omit() # OMIT NAN value
    
    cbps_hold = CBPS(formula(formula_hold),
                     data=data_hold,
                     exact = T)
    
    data_hold = data_hold %>%
      mutate(weights = cbps_hold$weights) %>%
      select(patientid, weights)
    
    bal = bal.tab(cbps_hold, un=TRUE)
    bal.df = bal$Balance.Across.Pairs %>%
      tibble::rownames_to_column('Variable') %>%
      select(-Type)
    
    results = list(weights = data_hold, balance_df = bal.df)
    
  }) 

model_data = lactate_baseline %>%
  left_join(propensity_weight_list %>%
              map(function(x) x %>%
                    keep(names(x) == 'weights')) %>%
              flatten() %>%
              set_names(paste0('cbpsdata', 2:6)) %>%
              reduce(left_join, by = 'patientid') %>%
              rename_at(grep('weight', colnames(.)), function(x) paste0('w', 2:6)), by = 'patientid') %>%
  mutate(log_acr = log(albumin_creatinine_ratio),
         log_lactate = log(lactate_normalized)) %>%
  mutate(survtime2 = pmin(sa_allc_ckdepi_time_renal1, 3),
         event2 = ifelse(sa_allc_ckdepi_time_renal1 > 3 | sa_allc_ckdepi_renal1 == 9, 0, 1)) %>%
  mutate(survtime3 = pmin(sa_allc_time_esrd, 3),
         event3 = ifelse(sa_allc_time_esrd > 3 | sa_allc_esrd %in% c(2,9), 0, 1),
         sa_allc_esrd2 = ifelse(sa_allc_esrd %in% c(2, 9), 0, sa_allc_esrd))

mixed_model_data = lactate_baseline_long %>%
  merge(model_data %>%
          select(c('patientid', 'w2')),
        by = 'patientid',
        all.x = T)

#################################
# Write Data
#################################

Write(model_data, 'surv_model_data.Rdata')

##################################################################
# Models  
# 
##################################################################

#################################
# Fit Cox Model with Tertiles
#################################

# Plot KM by Tertile
fit2 = survfit(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1)~lactate_tertile, data = model_data)

# fit survival model with weights
coxmod2 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile,
                weights=w2,
                data=model_data)

model_weights2 = coxmod2 %>%
  gtsummary::tbl_regression(exponentiate = T)

vif_tab = car::vif(coxmod2) %>%
  data.frame() %>%
  rename('GVIF^(1/(2*Df))'=3) %>%
  tibble::rownames_to_column('variable') %>%
  gt(rowname_col = 'variable')

#################################
# Event Summaries
#################################

event_summary1 = model_data %>%
  select(sa_allc_ckdepi_renal1, sa_allc_ckdepi_time_renal1) %>%
  filter(sa_allc_ckdepi_renal1 != 9) %>%
  tbl_summary(by = sa_allc_ckdepi_renal1)

event_summary2 = surv_median(fit2) %>%
  data.frame() %>%
  gt()

#################################
# Table one of tertile model
# variables
#################################

cov_tab1 = model_data %>%
  select(age,bmi,sex,log_acr,hemoglobin_a1c,metformin,map,egfr_ckd_epi,race_ethnicity_cat,lactate_tertile) %>%
  tbl_summary(by = 'lactate_tertile',
              missing = 'no') %>%
  add_p(pvalue_fun = ~ifelse(.x < 0.001, '<0.001', round(.x, 2)))

#################################
# Fit Cox Model with Lactate
# Tertile Only (Unweighted)
#################################

coxmod3 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ lactate_tertile,
                data=model_data)

model3 = coxmod3 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fit Cox Model with Lactate
# Tertile Only (Weighted)
#################################

coxmod4 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ lactate_tertile,
                data=model_data,
                weights = w2)

model_weights4 = coxmod4 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fully Adjusted Model
# (Unweighted)
#################################

# fit survival model with weights
coxmod5 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile,
                data=model_data)

model5 = coxmod5 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fully Adjusted Model
# (Weighted)
#################################

# fit survival model with weights
coxmod52 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+adenine+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile,
                 weights=w2,
                 data=model_data)

model52 = coxmod52 %>%
  gtsummary::tbl_regression(exponentiate = T,
                            pvalue_fun = ~ifelse(.x < 0.001, '<0.001', round(.x, 2)))

#################################
# Cox Model without ACR and
# eGFR (Weighted)
#################################

# fit survival model with weights
coxmod6 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+hemoglobin_a1c+metformin+map+race_ethnicity_cat+lactate_tertile,
                weights=w3,
                data=model_data)

model_weights6 = coxmod6 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Cox Model without Adenine 
# (Weighted)
#################################

# fit survival model with weights
coxmod7 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+log_acr+bmi+sex+hemoglobin_a1c+egfr_ckd_epi+metformin+map+race_ethnicity_cat+lactate_tertile,
                weights=w4,
                data=model_data)

model_weights7 = coxmod7 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Cox Model without eGFR
# (Weighted)
#################################

# fit survival model with weights
coxmod8 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+hemoglobin_a1c+log_acr+metformin+map+race_ethnicity_cat+lactate_tertile,
                weights=w5,
                data=model_data)

model_weights8 = coxmod8 %>%
  gtsummary::tbl_regression(exponentiate = T)

acr_lactate_correlation = data.frame(Correlation = cor(model_data$log_acr, model_data$lactate_normalized,
                                                       method = 'spearman')) %>%
  gt()

#################################
# Fully adjusted Cox-Model
# with weights (ESRD)
#################################

# fit survival model with weights
coxmod9 = coxph(Surv(sa_allc_time_esrd, sa_allc_esrd2) ~ age+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile,
                weights=w2,
                data=model_data)

model_weights9 = coxmod9 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fully adjusted Cox-Model
# adenine variable added
# (09/03/21)
#################################

# fit survival model with weights
coxmod_newvar = coxph(Surv(sa_allc_time_esrd, sa_allc_esrd2) ~ age+log_adenine2+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat,
                      data=model_data %>%
                        mutate(log_adenine2 = log(adenine+1, base = 2)))

coxmod_newvar2 = coxph(Surv(sa_allc_time_esrd, sa_allc_esrd2) ~ age+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat,
                       data=model_data)

model_newvar = coxmod_newvar %>%
  gtsummary::tbl_regression(exponentiate = T)

model_newvar2 = coxmod_newvar2 %>%
  gtsummary::tbl_regression(exponentiate = T)

#Note: Compare to full adjusted model

cindex_tab = coxmod_newvar$concordance %>%
  data.frame() %>%
  rename(vals = 1) %>%
  tibble::rownames_to_column('vars') %>%
  tidyr::pivot_wider(names_from = 'vars', values_from = 'vals') %>%
  mutate(model = 'Full') %>%
  rbind(coxmod_newvar2$concordance %>%
          data.frame() %>%
          rename(vals = 1) %>%
          tibble::rownames_to_column('vars') %>%
          tidyr::pivot_wider(names_from = 'vars', values_from = 'vals') %>%
          mutate(model = 'Adenine Removed')) %>%
  select(model, concordance, std)

##################################################################
# Subset Models
# 
##################################################################

subset_dataset = model_data %>%
  mutate(log_adenine2 = log(adenine+1, base = 2),
         race_ethnicity_cat = as.character(race_ethnicity_cat),
         african_american = ifelse(race_ethnicity_cat == 'black', 'african-american', 'other'),
         hispanic = ifelse(race_ethnicity_cat == 'hispanic', 'hispanic', 'non-hispanic'),
         acr_group = ifelse(is.na(albumin_creatinine_ratio), NA, 'unknown'),
         acr_group = ifelse(albumin_creatinine_ratio <= 30, 'normo', acr_group),
         acr_group = ifelse(albumin_creatinine_ratio > 30 & albumin_creatinine_ratio <= 300, 'micro', acr_group),
         acr_group = ifelse(albumin_creatinine_ratio > 300, 'macro', acr_group),
         acr_group = as.factor(acr_group),
         adenine_cat = ifelse(log_adenine2 < median(log_adenine2, na.rm = T), 'Lower', 'Higher'),
         sex2 = ifelse(sex == 1, 'Male', 'Female'),
         acr_sex_group = paste(acr_group, sex2, sep = '_'),
         log_adenine_std = log_adenine2/sd(log_adenine2, na.rm = T),
         adenine_cat2 = ifelse(log_adenine_std < median(log_adenine_std, na.rm = T), 'Lower', 'Upper'))

subset_groups = c('sex', 'acr_group')

subset_levels = list(c('Male', 'Female'), levels(subset_dataset$acr_group)) %>%
  set_names(subset_groups)

#################################
# Descriptive statistics
#################################

log_adenine_percentiles = subset_dataset$log_adenine_std %>%
  quantile(na.rm = T) %>%
  data.frame() %>%
  rename(Value = 1) %>%
  tibble::rownames_to_column('Percentile') %>%
  gt()

log_adenine_hist = subset_dataset %>%
  ggplot(aes(x = log_adenine_std))+
  geom_histogram()+
  theme_bw()

#################################
# Formula List
#################################

subset_models_covs1 = list('age+adenine_cat2+race_ethnicity_cat+bmi+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi',
                           'age+adenine_cat2+race_ethnicity_cat+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi') %>%
  set_names(subset_groups)

subset_models_covs_list = list(subset_models_covs1) %>%
  set_names(c("adenine_in"))

#################################
# Subset Models
#################################

all_subset_results = subset_models_covs_list %>%
  map(function(formula_list_element){
    
    subset_models_covs = formula_list_element
    
    subset_models = subset_groups %>%
      set_names() %>%
      map(function(nme) subset_dataset %>%
            group_by(eval(parse(text = nme))) %>%
            group_split() %>%
            map(function(data) coxph(formula(paste0('Surv(sa_allc_time_esrd, sa_allc_esrd2)~',
                                                    subset_models_covs[[nme]])), data = data)) %>%
            set_names(subset_levels[[nme]]))
    
    subset_dataset = subset_groups %>%
      set_names() %>%
      map(function(nme) subset_dataset %>%
            group_by(eval(parse(text = nme))) %>%
            group_split() %>%
            set_names(subset_levels[[nme]]))
    
    subset_summaries = subset_groups %>%
      map(function(nme) subset_models[[nme]] %>%
            map(function(coxobj) coxobj %>%
                  tbl_regression(exponentiate = T, include =  2)) %>%
            set_names(subset_levels[[nme]])) %>%
      set_names(subset_groups)
    
    c_index_table = subset_groups %>%
      set_names() %>%
      map(function(nme) subset_models[[nme]] %>%
            map(function(mod) mod$concordance %>%
                  data.frame() %>%
                  dplyr::rename(vals = 1) %>%
                  tibble::rownames_to_column('vars') %>%
                  tidyr::pivot_wider(names_from = 'vars', values_from = 'vals') %>%
                  select(concordance, std)) %>%
            bind_rows() %>%
            data.frame(model = subset_levels[[nme]]) %>%
            dplyr::select(model, everything())) %>%
      bind_rows()
    
    results = list(subset_summaries = subset_summaries,
                   c_index_table = c_index_table,
                   cox_obj = subset_models,
                   subset_dataset = subset_dataset)
    
    return(results)
    
  })

#################################
# Combined C-index Table
#################################

combined_c_index_table = all_subset_results$adenine_in$c_index_table %>%
  rename_at(2:ncol(.), function(x) paste0(x, '_in'))

#################################
# 
#################################

cox_coefs = all_subset_results %>%
  keep(names(all_subset_results) == 'adenine_in') %>%
  map(function(list_object) list_object %>%
        keep(names(list_object) == 'cox_obj')) %>%
  flatten() %>%
  flatten() %>%
  map(function(list_element1) list_element1 %>%
        map(function(list_element2) list_element2 %>%
              summary() %>%
              coef() %>%
              data.frame() %>%
              select(1) %>%
              tibble::rownames_to_column('variable') %>%
              filter(grepl('adenine', variable))) %>%
        bind_rows())

cox_cis = all_subset_results %>%
  keep(names(all_subset_results) == 'adenine_in') %>%
  map(function(list_object) list_object %>%
        keep(names(list_object) == 'cox_obj')) %>%
  flatten() %>%
  flatten() %>%
  map(function(list_element1) list_element1 %>%
        map(function(list_element2) list_element2 %>%
              confint() %>%
              data.frame() %>%
              tibble::rownames_to_column('variable') %>%
              filter(grepl('adenine', variable))) %>%
        bind_rows() %>%
        select(-variable) %>%
        dplyr::rename(lower = 1, upper = 2))

cox_coefs2 = names(cox_coefs) %>%
  set_names() %>%
  map(function(x) cox_coefs[[x]] %>%
        data.frame(levels = subset_levels[[x]]) %>%
        data.frame(cox_cis[[x]]) %>%
        select(variable, levels, everything())) %>%
  bind_rows() %>%
  select(-variable) %>%
  mutate_at(2:ncol(.), function(x) exp(x))

forest_plots = cox_coefs2 %>%
  mutate(levels = factor(levels, levels = subset_levels %>% flatten_chr())) %>%
  ggplot(aes(x = levels, y = coef, ymin = lower, ymax = upper))+
  geom_pointrange()+
  theme_bw()+
  scale_y_log10()+
  coord_flip()+
  labs(x = 'Level', y = 'HR')

##################################################################
# Interaction Models
# 
##################################################################

#################################
# Sex Interaction Model
#################################

sex_interaction_model = model_data %>%
  mutate(log_adenine2 = log(adenine+1, base = 2)) %>%
  (function(data) coxph(Surv(sa_allc_time_esrd, sa_allc_esrd2) ~ age+log_adenine2*sex+race_ethnicity_cat+bmi+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi,data=data)) %>%
  tbl_regression(exponentiate = T, include = 11)

#################################
# Fit Cox Model with acr_group1
# and lactate interaction
# (Weighted)
#################################

# fit survival model with weights
coxmod10 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile*acr_group1,
                 weights=w2,
                 data=model_data)

model_weights10 = coxmod10 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fit Cox Model with acr_group1
# and lactate interaction
# (Unweighted)
#################################

# fit survival model with weights
coxmod11 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile*acr_group1,
                 data=model_data)

model_weights11 = coxmod11 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fit Cox Model with acr_group2
# and lactate interaction
# (Weighted)
#################################

# fit survival model with weights
coxmod12 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile*acr_group2,
                 weights=w2,
                 data=model_data)

model_weights12 = coxmod12 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Fit Cox Model with acr_group1
# and lactate interaction
# (Unweighted)
#################################

# fit survival model with weights
coxmod13 = coxph(Surv(sa_allc_ckdepi_time_renal1, sa_allc_ckdepi_renal1) ~ age+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi+race_ethnicity_cat+lactate_tertile*acr_group2,
                 data=model_data)

model_weights13 = coxmod13 %>%
  gtsummary::tbl_regression(exponentiate = T)

#################################
# Table of mortality at year
# 15
#################################

mortality_15yr = model_data %>%
  select(sa_death) %>%
  table() %>%
  data.frame() %>%
  rename('Status' = 1) %>%
  filter(Status != 9) %>%
  mutate(Status = ifelse(Status == 0, 'Alive', 'Dead')) %>%
  mutate(Total = sum(Freq),
         Pct = round(100*Freq/Total, 2)) %>%
  select(-Total)

model_data %>%
  group_by(sa_death, sa_allc_ckdepi_renal1) %>%
  summarise(N  = n())

#################################
# Stack tables used in analysis
#################################

tbl_list = list(model3, model_weights7, model_weights2)

surv_tbl_stack = tbl_stack(tbl_list,
                           group_header = c('Tertile Only (Unweighted)', 'ACR Excluded Model', 'Fully Adjusted Model'))

##################################################################
# Mixed Models
#
##################################################################

#################################
# Model
#################################

# fit mixed model with weights
mixed_model_weights <- lmer(egfr_ckd_epi ~ age+bmi+sex+log_acr+hemoglobin_a1c+metformin+map+race_ethnicity_cat+lactate_tertile*year+(year|patientid),
                            weights=mixed_model_data$weights,
                            data=mixed_model_data)

mixed_model_coefs = mixed_model_weights %>%
  tbl_regression(pvalue_fun = ~ifelse(.x < 0.001, '<0.001', round(.x, 2)))

##################################################################
# Visualization
#
##################################################################

#################################
# eGFR plots
#################################

set.seed(200)

egfr_plot_data = lactate_baseline_long %>%
  dplyr::select(patientid, year, egfr_ckd_epi) %>%
  merge(model_data %>%
          dplyr::select(patientid, survtime2, event2)) %>%
  #Randomly select 10 subjects
  filter(patientid %in% sample(unique(patientid), replace = F, size = 20))

egfr_plot_data2 = unique(egfr_plot_data$patientid) %>%
  set_names() %>%
  map(function(pid){

    data_hold = egfr_plot_data %>%
      filter(patientid == pid)

    if(unique(data_hold$event2) == 1){

      data_hold = data_hold %>%
        filter(!is.na(egfr_ckd_epi))

      max_year = max(data_hold$year)

      if(max_year == 3){

        data_hold[which(data_hold$year == max_year),]$year = max(data_hold$survtime2)
        max_year = 2

      }

      data_hold = data_hold %>%
        rbind(data.frame(patientid = pid, year = (max_year + 1),
                         egfr_ckd_epi = 0, survtime2 = unique(data_hold$survtime2), event2 = 1))

    }

    return(data_hold)


  }) %>%
  bind_rows()

egfr_plot = egfr_plot_data2 %>%
  mutate(event_color = ifelse(event2 == 1, 'red', 'black')) %>%
  ggplot(aes(x = year, y = egfr_ckd_epi, group = patientid, color = I(event_color)))+
  geom_point()+
  geom_line()+
  theme_bw()+
  labs(x = 'Year', y = 'eGFR')

fail_tab = table((model_data %>%
                    filter(year == 0))$event2) %>%
  data.frame() %>%
  mutate(Var1 = ifelse(Var1 == 0, "Didn't Fail", 'Failure'),
         Proportion = round(Freq/nrow(model_data), 3)) %>%
  gt::gt(rowname_col = 'Var1')

#################################
# KM plots
#################################

tertile_km = ggsurvplot(fit2)$plot +
  scale_color_grey()

tertile_km_hazard = ggsurvplot(fit2)$data.survplot %>%
  mutate(surv_compliment = 1-surv) %>%
  ggplot(aes(x = time, y = surv_compliment, color = strata))+
  geom_step(size = 1, linetype = 1)+
  theme_minimal()+
  scale_color_grey()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top') +
  labs(x='Years From Baseline', y = 'Incidence')

#################################
# KM Plot (Singapore)
#################################

#Plot with same format as Singapore groups
mlabs = paste0('Tertile ', 1:3)
mpal = c('grey80', 'grey40', 'black') %>%
  purrr::set_names(mlabs)
mline = c('solid', 'dashed', 'solid') %>%
  purrr::set_names(mlabs)

tertile_km_singapore = ggsurvplot(fit2, fun = 'event', legend = 'top',
                                  legend.title = 'Lactate:', legend.labs = mlabs,
                                  palette = mpal, linetype = mline,
                                  pval = T, pval.method = T, pval.method.coord = c(0.1, 0.96),
                                  pval.coord=c(0.1, 0.90), pval.size = 5, surv.scale = 'percent', ylim = c(0, 1),
                                  ggtheme = theme_classic2(base_size = 16, base_family = 'Calibri'),
                                  font.family = 'Calibri') +
  labs(x = 'Years From Baseline', y = 'Cumulative incidence')

#################################
# Hazard Ratio Plots
#################################

coxmod_list = list(coxmod3, coxmod7, coxmod2)
coxmod_names = c('Tertile Only',  'Model Excluding ACR',  'Fully Adjusted Model')

coxmod_HR = coxmod_list %>%
  set_names(coxmod_names) %>%
  map(function(x) x %>%
        summary() %>%
        coef() %>%
        data.frame() %>%
        dplyr::select(2) %>%
        rename('HR'=1) %>%
        data.frame(x %>%
                     confint() %>%
                     data.frame() %>%
                     rename('lower' = 1, 'upper' = 2) %>%
                     mutate_all(exp))) %>%
  map(function(data) data %>%
        tibble::rownames_to_column('variable') %>%
        filter(grepl('lactate', variable)))

coxmod_HR2 = 1:length(coxmod_names) %>%
  map(function(i) coxmod_HR[[coxmod_names[i]]] %>%
        mutate(model = coxmod_names[i])) %>%
  set_names(coxmod_names) %>%
  bind_rows() %>%
  select(model, everything()) %>%
  mutate(variable = gsub('lactate_tertile', '', variable))

coxmod_forest = coxmod_HR2 %>%
  ggplot(aes(x = model, y = HR, ymin = lower, ymax = upper, color = variable))+
  geom_pointrange(position = position_dodge(0.5))+
  coord_flip()+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  scale_color_grey()+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.title.y= element_blank()) +
  labs(color = 'Tertile')

#################################
# Effects Plots
#################################

model_effects = allEffects(mixed_model_weights,
                           xlevels = list(year = c(0, 1, 2, 3)))

model_effects = as.data.frame(model_effects)$`lactate_tertile:year` %>%
  mutate(lower_se = fit - se, upper_se = fit + se)

effects_plot = ggplot(data=model_effects, aes(x=year, y=fit, color=lactate_tertile, ymin=lower, ymax=upper))+
  geom_pointrange(position=position_dodge(width=0.15), size=0.8, fatten = 2.5)+
  geom_line(position=position_dodge(width=0.15), size=0.8)+
  scale_color_grey()+
  scale_y_continuous(breaks = seq(30, 45, 5))+
  theme_classic()+
  theme(legend.position = 'bottom')+
  labs(x='Years From Baseline', y='eGFR', color='Lactate Tertile')

effects_plot2 = ggplot(data=model_effects, aes(x=year, y=fit, color=lactate_tertile, ymin=lower_se, ymax=upper_se))+
  geom_pointrange(position=position_dodge(width=0.15), size=0.8, fatten = 2.5)+
  geom_line(position=position_dodge(width=0.15), size=0.8)+
  scale_color_grey()+
  scale_y_continuous(breaks = seq(30, 45, 5))+
  theme_classic()+
  theme(legend.position = 'bottom')+
  labs(x='Years From Baseline', y='eGFR', color='Lactate Tertile')

#################################
# Effects Plots (Singapore)
#################################

rg = ref_grid(mixed_model_weights, at = list(year = 0:3))

emms = emmip(rg,lactate_tertile~year, CIs = TRUE, plotit = FALSE) %>%
  mutate(lactate_tertile = paste0('Tertile ',as.character(lactate_tertile)),
         lactate_tertile = as.factor(lactate_tertile))

Write(emms, 'plot_data.csv')

# Plot figure: means and 95% CIs
effects_plot_singapore1 = emms %>%
  ggplot(aes(x = year, y = yvar, group = lactate_tertile))+
  geom_line(color = "grey30")+
  geom_pointrange(aes(ymin = LCL, ymax = UCL),
                  position = position_dodge(0.1), fill = "white") +
  labs(shape = "Lactate:") +
  scale_shape_manual(values = c(21, 22, 24), labels = c("Tertile 1", "Tertile 2", "Tertile 3"))+
  scale_x_continuous(limits = c(-0.2, 3.2), breaks = c(0, 1, 2, 3)) +
  scale_y_continuous(limits = c(30, 45)) +
  xlab("Years From Baseline") +
  ylab(expression(eGFR~mL/min/1.73~m^{2})) +
  theme_classic() +
  theme(legend.position = "top")

effects_plot_singapore2 = emms %>%
  ggplot(aes(x = year, y = yvar, color = lactate_tertile)) +
  geom_pointrange(aes(ymin = yvar - SE, ymax = yvar + SE), position = position_dodge(0.1)) +
  geom_line(position=position_dodge(width=0.15), size=0.8) +
  scale_color_grey()+
  labs(color = "Lactate:") +
  scale_shape_manual(values = c(21, 22, 24), labels = c("Tertile 1", "Tertile 2", "Tertile 3"))+
  scale_x_continuous(limits = c(-0.2, 3.2), breaks = c(0, 1, 2, 3)) +
  scale_y_continuous(limits = c(30, 45)) +
  xlab("Years From Baseline") +
  ylab(expression(eGFR~mL/min/1.73~m^{2})) +
  theme_classic() +
  theme(legend.position = "top")

#################################
# Save figures in powerpoint
# file
#################################

path = paste0(getProjectPath(), '/Results/', source.file, '/editable_plots.pptx')

doc1 = read_pptx() %>%
  add_slide() %>%
  ph_with(value = dml(ggobj = effects_plot_singapore2),
          location = ph_location_type('body'))

print(doc1, target = path)

##################################################################
# Variable Comparison
#
##################################################################

sa_allc_var_comparisons = with(lactate_baseline, table(sa_allc_ckdepi_renal1, sa_allc_cric_renal1)) %>%
  data.frame() %>%
  gt()

time_death_table = with(lactate_baseline %>%
                          filter(sa_death == 1), table(sa_allc_ckdepi_time_renal1 < time_death)) %>%
  data.frame() %>%
  mutate(Variable = 'sa_allc_ckdepi_time_renal1') %>%
  rbind(with(lactate_baseline %>% filter(sa_death == 1), table(sa_allc_cric_time_renal1 < time_death)) %>%
          data.frame() %>%
          mutate(Variable = 'sa_allc_cric_time_renal1')) %>%
  mutate(Var1 = ifelse(Var1 == TRUE, 'Time Death < Time Event', 'Time Death > Time Event')) %>%
  select(Variable, everything()) %>%
  gt(rowname_col = 'Variable')

##################################################################
# Write Results
#
##################################################################

results_list = list(model_weights2 = model_weights2,
                    tertile_km = tertile_km,
                    tertile_km_hazard = tertile_km_hazard)

Write(results_list, 'SURVIVAL_RES.Rdata')

# End Program Body

dependency.out <- finalize_dependency()
