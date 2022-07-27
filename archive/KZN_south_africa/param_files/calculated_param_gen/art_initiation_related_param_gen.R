#calculate parameters used with model states to generate art initiation estimates
#ART coverage
#% of compartment 2 (CD4 < 200) eligible to initiate ART

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr', 'normalr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
#projections from HIV-HPV Cara's model
hiv_data_indir<-paste0(here(), '/param_files/calculated_param_gen/input_data/')

#update parameter file
outdir <- paste0(here(),'/param_files/')

setwd(hiv_data_indir)
hiv_input_df<-read.csv('hiv_input_gen_data_1980_2028.csv')
art_scenario_ref_df<-read.csv('art_coverage_scenario_ref.csv')

art_initiation_related_param_df<-hiv_input_df%>%
  #before 2010 not eligible to anyone in hiv compartment 2
  mutate(percent_hiv_compartment_2_eligible = if_else(year <= 2010, 0, 
  #between 2011-2015 some PLWH from hiv compartment 2 become eligible                                                     
                                                      if_else(year <=2015,
                                                              (PLHIV_2a_CD4_350_less/PLHIV_2_CD4_200_more),
  #in 2016 everyone becomes eligible                                                            
                                                              1)))%>%
  select(c('year', 'gender_name', 'gender_id', 'year_gender',
           'percent_hiv_compartment_2_eligible'))


#ART does not start until 2004
n_years_warmup_without_art<-14
n_years_warmup_with_art<-14

scale_up_male<-art_scenario_ref_df%>%
  filter(POLICY_ID == 1,
         G_SET == 1)

scale_up_male<-(scale_up_male$art_coverage)/n_years_warmup_with_art

scale_up_female<-art_scenario_ref_df%>%
  filter(POLICY_ID == 1,
         G_SET == 2)

scale_up_female<-(scale_up_female$art_coverage)/n_years_warmup_with_art


art_initiation_related_param_df<-art_initiation_related_param_df%>%
  filter(year <= 2018)%>%
  mutate(art_coverage = if_else(year <= 2004, 0,#before 2004 ART coverage did not exist
                                     #calibration period scale up to art coverage
                                     if_else(year<=2017 & gender_id == 1, 
                                             ((year-2003)*scale_up_male),
                                             if_else(year<=2017 & gender_id == 2,
                                                     ((year-2003)*scale_up_female),
                                                     if_else(year==2018 & gender_id == 1, 
                                                             ((n_years_warmup_with_art)*scale_up_male),
                                                             ((n_years_warmup_with_art)*scale_up_female))))))

art_initiation_related_param_df$POLICY_ID<-rep(1, times = nrow(art_initiation_related_param_df))

policy_2_art_coverage_male<-art_scenario_ref_df%>%
  filter(POLICY_ID == 2,
         G_SET == 1)

policy_2_art_coverage_male<-policy_2_art_coverage_male$art_coverage

#add 2018 for evaluation period for other policies
art_initiation_related_param_df<-rbind(art_initiation_related_param_df,
                                       c(2018, 'male', 1, '2018_male', 1, policy_2_art_coverage_male, 2))

policy_3_art_coverage_male<-art_scenario_ref_df%>%
  filter(POLICY_ID == 3,
         G_SET == 1)

policy_3_art_coverage_male<-policy_3_art_coverage_male$art_coverage

#add 2018 for evaluation period for other policies
art_initiation_related_param_df<-rbind(art_initiation_related_param_df,
                                       c(2018, 'male', 1, '2018_male', 1, policy_3_art_coverage_male, 3))



policy_2_art_coverage_female<-art_scenario_ref_df%>%
  filter(POLICY_ID == 2,
         G_SET == 2)

policy_2_art_coverage_female<-policy_2_art_coverage_female$art_coverage

#add 2018 for evaluation period for other policies
art_initiation_related_param_df<-rbind(art_initiation_related_param_df,
                                       c(2018, 'female', 2, '2018_female', 1, policy_2_art_coverage_female, 2))

policy_3_art_coverage_female<-art_scenario_ref_df%>%
  filter(POLICY_ID == 3,
         G_SET == 2)

policy_3_art_coverage_female<-policy_3_art_coverage_female$art_coverage

#add 2018 for evaluation period for other policies
art_initiation_related_param_df<-rbind(art_initiation_related_param_df,
                                       c(2018, 'female', 1, '2018_female', 1, policy_3_art_coverage_female, 3))
  
  
setwd(outdir)
write.csv(art_initiation_related_param_df, 'art_initiation_related_param_df.csv')
  
