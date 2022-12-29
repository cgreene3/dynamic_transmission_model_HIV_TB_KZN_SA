#calculated IPT initiations based on linear scale up in warm up (calibration period)
#in 2018 differentiated by program

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#'normalr'

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB/KZN_south_africa for here to work
indir<-paste0(here(), '/param_files/calculated_param_gen/raw_input_data/DO_ART')

#update parameter file
outdir <- paste0(here(),'/param_files/input_parameters')

setwd(indir)
IPT_scenario_df<-read.csv('IPT_initiation_scenario_ref.csv')


#####add in baseline scenario###
#linear scale up of art coverage in calibration period#


#for reps
n_gender = 2
n_years = 28
n_on_off_art_cats = 2

#IPT does not start until 2005
n_years_warmup_without_ipt<-15
n_years_warmup_with_ipt<-13

scale_up_male_yes<-IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 1,
         on_art == 'yes')

scale_up_male_yes<-(scale_up_male_yes$ipt_init_perc)/n_years_warmup_with_ipt

scale_up_male_no<-IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 1,
         on_art == 'no')

scale_up_male_no<-(scale_up_male_no$ipt_init_perc)/n_years_warmup_with_ipt

scale_up_female_yes<-IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 2,
         on_art == 'yes')

scale_up_female_yes<-(scale_up_female_yes$ipt_init_perc)/n_years_warmup_with_ipt

scale_up_female_no<-IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 2,
         on_art == 'no')

scale_up_female_no<-(scale_up_female_no$ipt_init_perc)/n_years_warmup_with_ipt

###create a dataframe for IPT initiation percent calculations in warmup period###
year <- rep(1990:2017, times = n_gender*n_on_off_art_cats)

gender_name<-c(rep("male", times = n_years*n_on_off_art_cats),
          rep("female", times = n_years*n_on_off_art_cats))

gender_id<-c(rep(1, times = n_years*n_on_off_art_cats),
             rep(2, times = n_years*n_on_off_art_cats))

on_art<-c(rep('yes', times = n_years),
          rep('no', times = n_years),
          rep('yes', times = n_years),
          rep('no', times = n_years))

scale_up<-c(rep(0, times = n_years_warmup_without_ipt),
            rep(scale_up_male_yes, times = n_years_warmup_with_ipt),
            rep(0, times = n_years_warmup_without_ipt),
            rep(scale_up_male_no, times = n_years_warmup_with_ipt),
            rep(0, times = n_years_warmup_without_ipt),
            rep(scale_up_female_yes, times = n_years_warmup_with_ipt),
            rep(0, times = n_years_warmup_without_ipt),
            rep(scale_up_female_no, times = n_years_warmup_with_ipt))

ipt_initiation_df<-data.frame(year,
                              gender_name,
                              gender_id,
                              on_art,
                              scale_up)

ipt_initiation_df<-ipt_initiation_df%>%
  mutate(ipt_init_perc = (year-2004)*scale_up)%>%
  select(-c('scale_up'))

#baseline (program 1) in calib period
ipt_initiation_df$program_id<-rep(1, times = nrow(ipt_initiation_df))

###add in program info

###program 1 on art males###
program_1_IPT_coverage_male_on_art = IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 1,
         on_art == 'yes')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'male', 1, 'yes', 
                           program_1_IPT_coverage_male_on_art$ipt_init_perc, 
                           1))

###program 1 on art females###
program_1_IPT_coverage_female_on_art = IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 2,
         on_art == 'yes')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'female', 2, 'yes', 
                           program_1_IPT_coverage_female_on_art$ipt_init_perc, 
                           1))

###program 1 not on art males###
program_1_IPT_coverage_male_not_on_art = IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 1,
         on_art == 'no')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'male', 1, 'no', 
                           program_1_IPT_coverage_male_not_on_art$ipt_init_perc, 
                           1))

###program 1 not on art females###
program_1_IPT_coverage_female_not_on_art = IPT_scenario_df%>%
  filter(program_id == 1,
         G_SET == 2,
         on_art == 'no')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'female', 2, 'no', 
                           program_1_IPT_coverage_male_not_on_art$ipt_init_perc, 
                           1))

###program 2
###program 2 on art males###
program_2_IPT_coverage_male_on_art = IPT_scenario_df%>%
  filter(program_id == 2,
         G_SET == 1,
         on_art == 'yes')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'male', 1, 'yes', 
                           program_2_IPT_coverage_male_on_art$ipt_init_perc, 
                           2))

###program 2 on art females###
program_2_IPT_coverage_female_on_art = IPT_scenario_df%>%
  filter(program_id == 2,
         G_SET == 2,
         on_art == 'yes')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'female', 2, 'yes', 
                           program_2_IPT_coverage_female_on_art$ipt_init_perc, 
                           2))

###program 2 not on art males###
program_2_IPT_coverage_male_not_on_art = IPT_scenario_df%>%
  filter(program_id == 2,
         G_SET == 1,
         on_art == 'no')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'male', 1, 'no', 
                           program_2_IPT_coverage_male_not_on_art$ipt_init_perc, 
                           2))

###program 2 not on art females###
program_2_IPT_coverage_female_not_on_art = IPT_scenario_df%>%
  filter(program_id == 2,
         G_SET == 2,
         on_art == 'no')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'female', 2, 'no', 
                           program_2_IPT_coverage_male_not_on_art$ipt_init_perc, 
                           2))

###program 3
###program 3 on art males###
program_3_IPT_coverage_male_on_art = IPT_scenario_df%>%
  filter(program_id == 3,
         G_SET == 1,
         on_art == 'yes')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'male', 1, 'yes', 
                           program_3_IPT_coverage_male_on_art$ipt_init_perc, 
                           3))

###program 3 on art females###
program_3_IPT_coverage_female_on_art = IPT_scenario_df%>%
  filter(program_id == 3,
         G_SET == 2,
         on_art == 'yes')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'female', 2, 'yes', 
                           program_3_IPT_coverage_female_on_art$ipt_init_perc, 
                           3))

###program 3 not on art males###
program_3_IPT_coverage_male_not_on_art = IPT_scenario_df%>%
  filter(program_id == 3,
         G_SET == 1,
         on_art == 'no')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'male', 1, 'no', 
                           program_3_IPT_coverage_male_not_on_art$ipt_init_perc, 
                           3))

###program 3 not on art females###
program_3_IPT_coverage_female_not_on_art = IPT_scenario_df%>%
  filter(program_id == 3,
         G_SET == 2,
         on_art == 'no')

ipt_initiation_df<-rbind(ipt_initiation_df,
                         c(2018, 'female', 2, 'no', 
                           program_3_IPT_coverage_male_not_on_art$ipt_init_perc, 
                           3))

#now assume all on IPT also on ART
ipt_initiation_df<-ipt_initiation_df%>%
  filter(on_art == "yes")%>%
  select(-c(on_art))

setwd(outdir)
write.csv(ipt_initiation_df, 'ipt_initiation_df.csv')
