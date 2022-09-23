#draws second samples fitting to beta distribution

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle'), require, character.only=T)

indir<-paste0(here(), '/calibration_analysis/LHS_1')

outdir <- paste0(here(),'/param_files/input_parameters')

setwd(indir)
best_calib_param_sets_df<-read.csv('best_results_df_summarised.csv')%>%
  filter(total_in_confidence == 20)
best_calib_param_sets_df<-best_calib_param_sets_df[10:ncol(best_calib_param_sets_df)]

setwd(outdir)
model_params_df<-read_excel('KZN_SA_model_parameters.xlsx')%>%
  select(c('model_matched_param', 'value', 'min', 'max', 'calibrated', 'TB_compartment', 
           'DR_compartment', 'HIV_compartment', 'G_compartment', 'notation'))

model_params_df$model_matched_param<-str_replace_all(model_params_df$model_matched_param, c("," = "."))


best_calib_param_sets_df_melt <- melt(best_calib_param_sets_df)
best_calib_param_sets_df_summ<-best_calib_param_sets_df_melt%>%
  group_by(variable)%>%
  summarise(min_tightened_norm = as.double(min(value)),
            max_tightened_norm = as.double(max(value)))%>%
  left_join(model_params_df, by = c('variable' = 'model_matched_param'))%>%
  mutate(min_tightened = as.double(((min_tightened_norm)*(max - min))+min))%>%
  mutate(max_tightened = as.double(((max_tightened_norm)*(max - min))+min))%>%
  select(c('variable', 'min_tightened', 'max_tightened'))


model_params_df_2<-model_params_df%>%
  left_join(best_calib_param_sets_df_summ, by = c('model_matched_param' = 'variable'))%>%
  select(c('model_matched_param', 'value', 'min_tightened', 'max_tightened', 'calibrated', 
           'TB_compartment', 
           'DR_compartment', 'HIV_compartment', 'G_compartment', 'notation'))

colnames(model_params_df_2)<-c('model_matched_param', 'value', 'min', 'max', 'calibrated',
                               'TB_compartment', 
                               'DR_compartment', 'HIV_compartment', 'G_compartment', 'notation')

setwd(outdir)
write.csv(model_params_df_2, 'KZN_SA_model_parameters_round_2.csv', row.names = FALSE)
