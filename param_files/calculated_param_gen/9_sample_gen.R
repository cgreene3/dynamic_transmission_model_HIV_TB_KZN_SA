#clean workspace
rm(list = ls())
gc()


library(readxl)
library(here)
library(dplyr)
library(ExtDist)
library(lhs)
library(ggplot2)
library(stringr)
library(reshape2)

n_samples <- 100000
sample_df <- data.frame(matrix(nrow = n_samples, ncol = 0))

set.seed(as.integer(1)) 

#update parameter file
outdir <- paste0(here(),'/param_files/input_parameters')
indir_params <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data')

setwd(indir_params)
model_params_df<-read_excel('KZN_SA_model_parameters.xlsx')%>%
  select(c('model_matched_param', 'min', 'max', 'value', 'calibrated'))

model_params_df$model_matched_param<-str_replace_all(model_params_df$model_matched_param, c("," = "."))

n_calib_params<-sum(model_params_df$calibrated == 'yes')
calib_param_names<-model_params_df%>%
  filter(calibrated == 'yes')%>%
  select(c('model_matched_param'))
calib_param_names<-as.vector(t(calib_param_names))

raw_lhs_df <- as.data.frame(randomLHS(n_samples, n_calib_params))
colnames(raw_lhs_df)<-calib_param_names

for (mp in model_params_df$model_matched_param){
  if(model_params_df[model_params_df$model_matched_param == mp, 
                     c('calibrated')] == "yes"){
    
    max_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                              c('max')])
    min_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                              c('min')])
    
    sample_df_temp<-as.data.frame(qunif(unlist(raw_lhs_df%>%select(mp)),
                    min_temp, max_temp))
    colnames(sample_df_temp)<-mp
    sample_df<-cbind(sample_df, sample_df_temp)
  } else {
    value_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                         c('value')])
    sample_df_temp<-as.data.frame(rep(value_temp, times = n_samples))
    colnames(sample_df_temp)<-mp
    sample_df<-cbind(sample_df, sample_df_temp)
    
  }
}
sample_df$sim_id <-1:nrow(sample_df)

setwd(outdir)
write.csv(sample_df, 'calibration_sets_df.csv', row.names = FALSE)

