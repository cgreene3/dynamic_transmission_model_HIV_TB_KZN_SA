#draws initial samples using latin hypercube

#clean workspace
rm(list = ls())
gc()

library(lhs)
library(readxl)
library(here)
library(dplyr)

n_samples <- 100000

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


set.seed(1)
raw_lhs_df_1 <- as.data.frame(randomLHS(n_samples, n_calib_params))
colnames(raw_lhs_df_1)<-calib_param_names
raw_lhs_df_1$sim_id<-1:nrow(raw_lhs_df_1)

set.seed(2)
raw_lhs_df_2 <- as.data.frame(randomLHS(n_samples, n_calib_params))
colnames(raw_lhs_df_2)<-calib_param_names
raw_lhs_df_2$sim_id<-1:nrow(raw_lhs_df_2)

write.csv(raw_lhs_df_1, 'sim_calibration_ref_df_LHS_1.csv', row.names = FALSE)
write.csv(raw_lhs_df_2, 'sim_calibration_ref_df_LHS_2.csv', row.names = FALSE)


