
rm(list=ls())
set.seed( 2011 )

library("adapr")
library('lme4')
library('lmerTest')
library('reshape2')
library('ggplot2')
library('rtf')
library('MatchIt')
library('glmnet')
library('CBPS')
library('dplyr')
library('mitools')
library('effects')
library('gelfondjal')
library("dplyr")
library("magrittr")
library("tidyverse")

source.file <-"preprocess_data.r"
project.id <- "Sharma_653"
source_info <- create_source_file_dir(source.description="Process data")

# Program body here

setwd(paste0(getProjectPath(), '/Programs'))

##################################################################
# Comments for reproducing analysis
#
##################################################################

# Note: run browseURL(dataDir()) to open data directory folder

# Note: Generating HTML report file:
message(paste0('Note: To generate final report place final_results.Rmd file in ',
               paste0(getProjectPath(), '/Programs/Markdown')))

# Note: To obtain final results run the following sequence of scripts:
# preprocess_data.R, mixed_model_update3.R, and final_results.R
# Make sure that final_results.Rmd is in markdown folder before doing this

##################################################################
# Create datasets
#
##################################################################

#################################
# Subsetting variables
#################################

baseline_variables = c('age', 'bmi', 'hemoglobin_a1c', 'map', 'race_ethnicity_cat',
                       'sex', 'albumin_creatinine_ratio', 'metformin')

quantile_variables = c('lactate_normalized', 'egfr_cric_v')

#################################
# Baseline and Full datasets
#################################

long_df = read.csv('/Users/rabiul/fwdcodeusedinanalysis/CRIC-Lactate-Data-Nov2020.csv') %>%
  rename_with(~tolower(.x)) %>%
  mutate(race_ethnicity_cat = plyr::revalue(as.factor(race_ethnicity_cat2), c('1' = 'white', '2' = 'black', '3' = 'hispanic', '4' = 'other')),
         year = (vnum - 3)/2) %>%
  filter(year %in% 0:3) %>%
  rename(patientid = pid, albumin_creatinine_ratio = uacratio, egfr_cric_v = egfr_cric,
         age = age_integer, egfr_ckd_epi = egfr_roche) %>%
  arrange(patientid, year) %>%
  dplyr::select(-race_ethnicity_cat2)

lactate_baseline = long_df %>%
  filter(year == 0) %>%
  dplyr::select(c(baseline_variables, quantile_variables, 'patientid')) %>%
  mutate(lactate_quart = cut(lactate_normalized,
                             breaks = quantile(lactate_normalized,probs=seq(0, 1, by = 0.25),na.rm=T,include.lowest=TRUE)),
         lactate_tertile = cut(lactate_normalized,
                               breaks = quantile(lactate_normalized,probs=seq(0, 1, by = 1/3),na.rm=T,include.lowest=TRUE)),
         lactate_binary = ifelse(lactate_normalized <= median(lactate_normalized), 'low', 'high'),
         egfr_cric_bin = factor(ifelse(egfr_cric_v <= median(egfr_cric_v), 'low', 'high'))) %>%
  mutate(lactate_quart = plyr::revalue(lactate_quart,c('(0.141,12.6]'='Q1', '(12.6,23.4]'='Q2',
                                                       '(23.4,44.9]'='Q3', '(44.9,1.19e+03]'='Q4')),
         lactate_tertile = recode(lactate_tertile, "(0.141,15.9]"='T1', "(15.9,35.3]"='T2',
                                  "(35.3,1.19e+03]"='T3')) %>%
  mutate(acr_group1 = ifelse(albumin_creatinine_ratio <= 500, '0-500', '>500'),
         acr_group1 = ifelse(is.na(albumin_creatinine_ratio), NA, acr_group1),
         acr_group1 = relevel(as.factor(acr_group1), ref = '0-500'),
         acr_group2 = ifelse(albumin_creatinine_ratio <= 500, '0-500', '500-1000'),
         acr_group2 = ifelse(albumin_creatinine_ratio > 1000, '>1000', acr_group2),
         acr_group2 = ifelse(is.na(albumin_creatinine_ratio), NA, acr_group2),
         acr_group2 = factor(acr_group2, levels = c('0-500', '500-1000', '>1000'))) %>%
  dplyr::select(-quantile_variables)

#Final dataset
long_df = long_df %>%
  select(-c(baseline_variables)) %>%
  merge(lactate_baseline, by = 'patientid') %>%
  mutate(log_acr = log(albumin_creatinine_ratio))

#Subsets used in one year and three year models
all_datasets = lapply(list(0:1, 0:3), function(x) long_df %>% filter(year %in% x))
names(all_datasets) = c('one_yr', 'three_yr')

#################################
# Write data
#################################

write(all_datasets,'all_data.Rdata')
write.csv(all_datasets,'~/Desktop/all_data.csv')

# End Program Body

dependency.out <- finalize_dependency()

