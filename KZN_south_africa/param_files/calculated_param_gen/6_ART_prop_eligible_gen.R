#calculate parameters used with model states to generate art initiation estimates
#% of compartment 2 (CD4 < 200) eligible to initiate ART

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 
         'ggplot2', 'stringr', 'normalr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work

#projections from HIV-HPV Cara's model
indir<-paste0(here(), '/param_files/calculated_param_gen/input_data/HPV_HIV_model')

#update parameter file
outdir <- paste0(here(),'/param_files/')

setwd(indir)
hiv_input_df<-read.csv('ART_prop_eligible_1980_2028_HPV_HIV.csv')

prop_eligible_df<-hiv_input_df%>%
  #before 2010 not eligible to anyone in hiv compartment 2
  mutate(percent_hiv_compartment_2_eligible = if_else(year <= 2010, 0, 
                                                      #between 2011-2015 some PLWH from hiv compartment 2 become eligible                                                     
                                                      if_else(year <=2015,
                                                              (PLHIV_2a_CD4_350_less/PLHIV_2_CD4_200_more),
                                                              #in 2016 everyone becomes eligible                                                            
                                                              1)))%>%
  select(c('year', 'gender_name', 'gender_id', 'year_gender',
           'percent_hiv_compartment_2_eligible'))

prop_eligible_df$max<-prop_eligible_df$percent_hiv_compartment_2_eligible*1.25
prop_eligible_df$min<-prop_eligible_df$percent_hiv_compartment_2_eligible*.75


setwd(outdir)
write.csv(prop_eligible_df, 'ART_prop_eligible_df.csv', row.names = FALSE)
