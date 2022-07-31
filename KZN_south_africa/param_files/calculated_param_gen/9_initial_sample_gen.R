#draws initial samples using latin hypercube

#clean workspace
rm(list = ls())
gc()

library(lhs)
library(readxl)
library(here)
library(dplyr)

set.seed(1)
n_samples <- 10

#update parameter file
indir_outdir <- paste0(here(),'/param_files/input_parameters')

setwd(indir_outdir)
model_params_df<-read_excel('KZN_SA_model_parameters.xlsx')%>%
  select(c('model_matched_param', 'min', 'max', 'calibrated'))

n_calib_params<-sum(model_params_df$calibrated == 'yes')
calib_param_names<-model_params_df%>%
  filter(calibrated == 'yes')%>%
  select(c('model_matched_param'))
calib_param_names<-as.vector(t(calib_param_names))
raw_lhs_df <- as.data.frame(randomLHS(n_samples, n_calib_params))
colnames(raw_lhs_df)<-calib_param_names
raw_lhs_df$sim_id<-1:nrow(raw_lhs_df)

write.csv(raw_lhs_df, 'sim_calibration_ref_df.csv', row.names = FALSE)
