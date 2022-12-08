
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

source.file <-"survival_model2.r"
project.id <- "Sharma_653"
source_info <- create_source_file_dir(source.description="")

# Program body here

setwd(paste0(getProjectPath(), '/Programs'))

##################################################################
# Load Data
# 
##################################################################
#rawdata <- read.csv("RedCapUp2.csv")#, stringsAsFactors = TRUE)

model_data = Load.branch('surv_model_data.Rdata')

##################################################################
# Functions
#
##################################################################

map_to_tertile = function(x){
  
  quartiles = quantile(x, na.rm = T, probs = seq(0, 1, by = 1/3))
  
  vect_return = rep(NA, length(x)) %>%
    (function(vect) ifelse(x >= quartiles[1] & x < quartiles[2], 'Q1', vect)) %>%
    (function(vect) ifelse(x >= quartiles[2] & x < quartiles[3], 'Q2', vect)) %>%
    (function(vect) ifelse(x >= quartiles[3], 'Q3', vect)) %>% 
    as.factor()
  
  return(vect_return)
}

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
         adenine_tertile = map_to_tertile(log_adenine2),
         sex2 = ifelse(sex == 1, 'Male', 'Female'),
         acr_sex_group = paste(acr_group, sex2, sep = '_'))

subset_groups = c('sex', 'acr_group')

subset_levels = list(c('Male', 'Female'), levels(subset_dataset$acr_group)) %>%
  set_names(subset_groups)

#################################
# Descriptive statistics
#################################

adenine_cat_table1 = subset_dataset %>%
  select(acr_group, sex2, adenine_cat) %>%
  tbl_summary(by = 'adenine_cat')

adenine_cat_table2 = subset_dataset %>%
  select(acr_group, sex2, adenine_tertile) %>%
  tbl_summary(by = 'adenine_tertile')

adenine_cat_table_comb = tbl_merge(list(adenine_cat_table1, adenine_cat_table2),
                                   tab_spanner = c('Binary', 'Tertile'))

#################################
# Formula List
#################################

subset_models_covs1 = list('age+log_adenine2+race_ethnicity_cat+bmi+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi',
                           'age+log_adenine2+race_ethnicity_cat+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi') %>%
  set_names(subset_groups)

subset_models_covs2 = list('age+adenine_cat+race_ethnicity_cat+bmi+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi',
                           'age+adenine_cat+race_ethnicity_cat+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi') %>%
  set_names(subset_groups)

subset_models_covs3 = list('age+adenine_tertile+race_ethnicity_cat+bmi+log_acr+hemoglobin_a1c+metformin+map+egfr_ckd_epi',
                           'age+adenine_tertile+race_ethnicity_cat+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi') %>%
  set_names(subset_groups)

subset_models_covs_list = list(subset_models_covs1, subset_models_covs2, subset_models_covs3) %>%
  set_names(c("adenine_log2", 'adenine_binary', 'adenine_tertile'))


##Rabi update
coxph(formula=Surv(subset_dataset$sa_allc_time_esrd, subset_dataset$sa_allc_esrd2)~
                     age+log_adenine2+bmi+sex+hemoglobin_a1c+map+egfr_ckd_epi+log_acr, data = subset_dataset)



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


subset_summaries_list = all_subset_results %>%
  map(function(x) x %>%
        keep(names(x) == 'subset_summaries') %>%
        flatten()) %>%
  map(function(x) x %>%
        map(function(y) y %>% keep(names(y) != 'macro'))) %>%
  map(function(x) x %>%
        map(function(y) y %>%
              tbl_merge(tab_spanner = names(y))))

sex_subsets_table = subset_summaries_list %>%
  map(function(x) x %>%
        keep(names(x) == 'sex')) %>%
  flatten() %>%
  tbl_stack()

acr_subsets_table = subset_summaries_list %>%
  map(function(x) x %>%
        keep(names(x) == 'acr_group')) %>%
  flatten() %>%
  tbl_stack()

#################################
# Get Kaplan Meier Plots 
#################################

acr_subsets = subset_dataset %>%
  group_split(acr_group) %>%
  set_names(levels(subset_dataset$acr_group))

adjusted_km_plots = levels(subset_dataset$acr_group) %>%
  set_names() %>%
  map(function(nme) coxph(Surv(sa_allc_time_esrd, sa_allc_esrd2)~age+strata(adenine_tertile)+race_ethnicity_cat+bmi+sex+hemoglobin_a1c+metformin+map+egfr_ckd_epi,
                          data = acr_subsets[[nme]]) %>%
        (function(x) ggadjustedcurves(x, data = acr_subsets[[nme]] %>%
                                        data.frame(),
                                      method = 'conditional')+
           labs(title = nme)))

# Save plots as powerpoint

path = paste0(getProjectPath(), '/Results/', source.file, '/editable_plots.pptx')

doc1 = read_pptx() %>%
  add_slide() %>%
  ph_with(value = dml(ggobj = adjusted_km_plots$micro),
          location = ph_location_type('body')) %>%
  add_slide() %>%
  ph_with(value = dml(ggobj = adjusted_km_plots$normo),
          location = ph_location_type('body'))


print(doc1, target = path)

# End Program Body

dependency.out <- finalize_dependency()
