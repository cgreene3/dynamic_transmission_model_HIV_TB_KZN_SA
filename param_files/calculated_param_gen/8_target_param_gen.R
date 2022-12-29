#Estimating active TB mortalities/TB incidence/HIV prev per 100K to calibrate to
#Uses GBD
#KZN population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa
indir <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data/GBD')
outdir <- paste0(here(),'/param_files/target_calibration_estimates')

##########read in pop files#################
#data pull date 02/21/2022
calib_years <-c(2005, 2017)

setwd(indir)

#get popestimates
pop_df<-read.csv('pop_estimates_15_59.csv')

#input mortality rate calibration targets
TB_mort_num_df<-read.csv('disease_mort_num.csv')

#TB incidence rate calibration targets
TB_inc_num_df<-read.csv('tb_incidence_num.csv')

#HIV prevalence calibration targets
HIV_prev_num_df<-read.csv('hiv_prev_num.csv')

#TB mort calibration estimates 2005 and 2017
TB_mort_calibration_df<-TB_mort_num_df%>%
  filter(cause != "HIV/AIDS resulting in other diseases")%>%
  mutate(TB_HIV_coinfection = if_else(grepl('HIV', cause, ignore.case = TRUE), 'yes', 'no'))%>%
  group_by(year, sex, TB_HIV_coinfection)%>%
  summarise(val_mort_num = sum(val),
            upper_mort_num = sum(upper),
            lower_mort_num = sum(lower))%>%
  left_join(pop_df, by = c('year', 'sex'))%>%
  mutate(val_rate = (val_mort_num/expected_total_pop)*100000,
         upper_rate = (upper_mort_num/expected_total_pop)*100000,
         lower_rate = (lower_mort_num/expected_total_pop)*100000)%>%
  select(c('year', 'sex', 'TB_HIV_coinfection', 
           'val_rate', 'upper_rate', 'lower_rate'))

#TB incidence calibration estimates
TB_inc_calibration_df<-TB_inc_num_df%>%
  mutate(TB_HIV_coinfection = if_else(grepl('HIV', cause, ignore.case = TRUE), 'yes', 'no'))%>%
  group_by(year, sex, TB_HIV_coinfection)%>%
  summarise(val_inc_num = sum(val),
            upper_inc_num = sum(upper),
            lower_inc_num = sum(lower))%>%
  left_join(pop_df, by = c('year', 'sex'))%>%
  mutate(val_rate = (val_inc_num/expected_total_pop)*100000,
         upper_rate = (upper_inc_num/expected_total_pop)*100000,
         lower_rate = (lower_inc_num/expected_total_pop)*100000)%>%
  select(c('year', 'sex', 'TB_HIV_coinfection', 
           'val_rate', 'upper_rate', 'lower_rate'))

#HIV prev calibration estimates 2005 and 2017
HIV_prev_calibration_df<-HIV_prev_num_df%>%
  group_by(year, sex)%>%
  summarise(val_prev_num = sum(val),
            upper_prev_num = sum(upper),
            lower_prev_num = sum(lower))%>%
  left_join(pop_df, by = c('year', 'sex'))%>%
  mutate(val_rate = (val_prev_num/expected_total_pop)*100000,
         upper_rate = (upper_prev_num/expected_total_pop)*100000,
         lower_rate = (lower_prev_num/expected_total_pop)*100000)%>%
  select(c('year', 'sex', 
           'val_rate', 'upper_rate', 'lower_rate'))

setwd(outdir)
write.csv(TB_mort_calibration_df, 'KZN_GBD_TB_mort_rate_calibration_df.csv', row.names = FALSE)
write.csv(TB_inc_calibration_df, 'KZN_GBD_TB_inc_rate_calibration_df.csv', row.names = FALSE)
write.csv(HIV_prev_calibration_df, 'KZN_GBD_HIV_prev_rate_calibration_df.csv', row.names = FALSE)

