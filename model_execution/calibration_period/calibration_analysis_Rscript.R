#last updated sept 29, 2022

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle'), require, character.only=T)


####HYAK OR GITHUB SPECIFIC CODES TO COMMENT/UNCOMMENT####
###HYAK###
# #location where we sent outputs
indir_outputs<-'/gscratch/icrc/cgreene3/calibration_outputs/'
indir_input_params<-'/gscratch/icrc/cgreene3/input_parameters'
indir_target_calibration_estimates<-'/gscratch/icrc/cgreene3/target_calibration_estimates'

#location where we send calib analysis results
outdir<- '/gscratch/icrc/cgreene3/calibration_analysis/'

##GITHUB/LOCAL###
#location where we sent outputs
# library(here)
# indir_outputs<-paste0(here(), '/calibration_outputs_test/run_',run_id)
# indir_target_calibration_estimates<-paste0(here(), '/param_files/target_calibration_estimates')
# indir_input_params<-paste0(here(), '/param_files/input_parameters')
# 
# #location where we send calib analysis results
# outdir<- paste0(here(), '/calibration_analysis/run_',run_id)

setwd(indir_target_calibration_estimates)
HIV_prev_df<-read.csv('KZN_GBD_HIV_prev_rate_calibration_df.csv')%>%
  mutate(sex = as.character(sex))
TB_inc_df<-read.csv('KZN_GBD_TB_inc_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))
TB_mort_df<-read.csv('KZN_GBD_TB_mort_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))

setwd(indir_input_params)
sim_calibration_ref_df<-read.csv('calibration_sets_df.csv')

##combine simulation estimates##
setwd(indir_outputs)
all_files_in_outdir <- list.files(pattern="calib_metrics*")

outputs_combined_df<-read.csv(all_files_in_outdir[1])
for (i in 2:length(all_files_in_outdir)){
  print(all_files_in_outdir[i])
  temp<-read.csv(all_files_in_outdir[i])
  outputs_combined_df<-rbind(outputs_combined_df, temp)
}

missing_sims_df<-data.frame(sim_id = 1:100000)%>%
  left_join(outputs_combined_df, by = c('sim_id'))%>%
  filter(is.na(year))

outputs_combined_df_reshape<-melt(outputs_combined_df, id = c("year", "sim_id"))
outputs_combined_df_reshape<-outputs_combined_df_reshape%>%
  mutate(sex = if_else(grepl('female', variable), 'Female', 'Male'),
         TB_HIV_coinfection = if_else(grepl('neg', variable), 'no', 'yes'))

outputs_combined_df_reshape<-outputs_combined_df_reshape%>%
  filter(grepl('ppl', variable) == FALSE)%>%
  filter(grepl('pop', variable) == FALSE)%>%
  mutate(variable = as.character(variable))

outputs_TB_mort<-outputs_combined_df_reshape%>%
  filter(grepl('mort', variable) == TRUE)%>%
  mutate(variable = as.character(variable),
         sex = as.character(sex),
         TB_HIV_coinfection = as.character(TB_HIV_coinfection))%>%
  left_join(TB_mort_df, by = c('year', 'sex', 'TB_HIV_coinfection'))%>%
  mutate(measure = 'TB_mortality')%>%
  mutate(in_confidence_interval = if_else(((value < upper_rate) & (value > lower_rate)),
                                     1, 0))%>%
  #mutate(HIV_prev_only_CI = 0)%>%
  #mutate(chi_sq = ((value - val_mort_rate)**2)/val_mort_rate)%>%
  #mutate(TB_only_chi_sq = ((value - val_mort_rate)**2)/val_mort_rate)%>%
  #mutate(TB_mort_only_chi_sq = ((value - val_mort_rate)**2)/val_mort_rate)%>%
  #mutate(TB_inc_HIV_prev_only_CI = 0)%>%
  select(c('sim_id', 'year', 'variable', 'sex', 'measure', 'TB_HIV_coinfection',
           'value', 'in_confidence_interval', 'upper_rate', 'lower_rate'))
  
outputs_TB_inc<-outputs_combined_df_reshape%>%
  filter(grepl('inc', variable) == TRUE)%>%
  mutate(value = value)%>%
  left_join(TB_inc_df, by = c('year', 'sex', 'TB_HIV_coinfection'))%>%
  mutate(in_confidence_interval = if_else(((value < upper_rate) & (value > lower_rate)),
                                     1, 0))%>%
  mutate(measure = 'TB_incidence')%>%
  #mutate(HIV_prev_only_CI = 0)%>%
  #mutate(chi_sq = ((value - val_inc_rate)**2)/val_inc_rate)%>%
  #mutate(TB_only_chi_sq = ((value - val_inc_rate)**2)/val_inc_rate)%>%
  #mutate(TB_mort_only_chi_sq = 0)%>%
  #mutate(TB_inc_HIV_prev_only_CI = if_else(((value < upper_inc_rate) & (value > lower_inc_rate)),
  #                                               1, 0))%>%
  select(c('sim_id', 'year', 'variable', 'sex', 'measure', 'TB_HIV_coinfection',
           'value', 'in_confidence_interval', 'upper_rate', 'lower_rate'))

outputs_HIV_prev<-outputs_combined_df_reshape%>%
  select(-c('TB_HIV_coinfection'))%>%
  filter(grepl('hiv_prev', variable) == TRUE)%>%
  mutate(value = value)%>%
  left_join(HIV_prev_df, by = c('year', 'sex'))%>%
  mutate(in_confidence_interval = if_else(((value < upper_rate) & (value > lower_rate)),
                                     1, 0))%>%
  mutate(measure = 'HIV_prevalence')%>%
  mutate(TB_HIV_coinfection = 'NA')%>%
  #mutate(HIV_prev_only_CI = if_else(((value < upper_prev_rate) & (value > lower_prev_rate)),
  #                                  1, 0))%>%
  #mutate(chi_sq = ((value - val_prev_rate)**2)/val_prev_rate)%>%
  #mutate(TB_only_chi_sq = 0)%>%
  #mutate(TB_mort_only_chi_sq = 0)%>%
  #mutate(TB_inc_HIV_prev_only_CI = if_else(((value < upper_prev_rate) & (value > lower_prev_rate)),
  #                                         1, 0))%>%
  select(c('sim_id', 'year', 'variable', 'sex', 'measure', 'TB_HIV_coinfection',
           'value', 'in_confidence_interval', 'upper_rate', 'lower_rate'))

##combine and analyze
results_df<-rbind(outputs_TB_inc, outputs_TB_mort, outputs_HIV_prev)

results_df_summarised_selection_options<-results_df%>%
  group_by(sim_id)%>%
  summarise(total_in_confidence = sum(in_confidence_interval)#,
            #total_in_confidence_HIV_prev = sum(HIV_prev_only_CI),
            #total_TB_inc_HIV_prev_only_CI = sum(TB_inc_HIV_prev_only_CI),
            #total_TB_inc_only_CI = total_TB_inc_HIV_prev_only_CI-total_in_confidence_HIV_prev,
            #total_TB_mort_only_CI = total_in_confidence-total_TB_inc_HIV_prev_only_CI,
            #total_chi_sq_TB_only = sum(TB_only_chi_sq),
            #total_chi_sq_TB_mort_only = sum(TB_mort_only_chi_sq)
            )%>%
  left_join(sim_calibration_ref_df, by = c('sim_id'))

#lr_total_confidence <- lm(formula = total_in_confidence ~ . - sim_id, 
#                          data = results_df_summarised_selection_options%>%
#                            select(-c('total_in_confidence_HIV_prev',
#                                      'total_chi_sq_TB_only', 
#                                      'total_TB_inc_HIV_prev_only_CI', 
#                                      'total_chi_sq_TB_mort_only')))

#lr_TB_chi_sq_hiv_prev <-lm(formula = total_chi_sq_TB_only ~ . - sim_id, 
#                           data = results_df_summarised_selection_options%>%
#                             select(-c('total_in_confidence_HIV_prev',
#                                       'total_in_confidence', 
#                                       'total_TB_inc_HIV_prev_only_CI', 
#                                       'total_chi_sq_TB_mort_only')))

#if you want to look at the distribution of each parameter for sets that hit 90% of targets
almost_best_calibration_sets_ref_df<-results_df_summarised_selection_options%>%
  filter(total_in_confidence >= 19)

#accepted parameter sets
accepted_calibration_sets_ref_df<-results_df_summarised_selection_options%>%
  filter(total_in_confidence == 20)

#to summarise metrics of accepted points
accepted_calibration_metrics_df<-results_df%>%
  filter(sim_id %in% accepted_calibration_sets_ref_df$sim_id)

#to evaluate metrics that are missed by sets that hit at least 90% of metrics
almost_calibration_metrics_df<-results_df%>%
  filter(sim_id %in% almost_best_calibration_sets_ref_df$sim_id)

setwd(outdir)
write.csv(accepted_calibration_sets_ref_df, 'accepted_calibration_sets_ref_df.csv', row.names = FALSE)
write.csv(almost_best_calibration_sets_ref_df, 'almost_best_calibration_sets_ref_df.csv', row.names = FALSE)
write.csv(almost_calibration_metrics_df, 'almost_calibration_metrics_df.csv', row.names = FALSE)
write.csv(accepted_calibration_metrics_df, 'accepted_calibration_metrics_df.csv', row.names = FALSE)
write.csv(missing_sims_df, 'missing_sims_df.csv', row.names = FALSE)

setwd(indir_input_params)
write.csv(missing_sims_df, 'missing_sims_df.csv', row.names = FALSE)


#sink("lm_total_in_confidence.txt")
#print(summary(lr_total_confidence))
#sink()
#sink("lm_lr_TB_chi_sq_hiv_prev.txt")
#print(summary(lr_TB_chi_sq_hiv_prev))
#sink()
