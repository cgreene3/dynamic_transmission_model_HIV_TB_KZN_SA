#calculate parameters used with model states to generate art initiation estimates
#ART coverage

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#, 'normalr'

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
indir<-paste0(here(), '/param_files/calculated_param_gen/raw_input_data/DO_ART')

#update parameter file
outdir <- paste0(here(),'/param_files/input_parameters')

setwd(indir)
art_scenario_ref_df<-read.csv('art_coverage_scenario_ref.csv')

#ART does not start until 2004
n_years_warmup_without_art<-14
n_years_warmup_with_art<-14

scale_up_male<-art_scenario_ref_df%>%
  filter(program_id == 1,
         G_SET == 1)

scale_up_male<-(scale_up_male$art_coverage)/n_years_warmup_with_art

scale_up_female<-art_scenario_ref_df%>%
  filter(program_id == 1,
         G_SET == 2)

scale_up_female<-(scale_up_female$art_coverage)/n_years_warmup_with_art

year <- c(1980:2028, 1980:2028)
sex <- c(rep('Male', times = length(1980:2028)), rep('Female', times = length(1980:2028)))

art_coverage_df<-data.frame(year, sex)

art_coverage_df<-art_coverage_df%>%
  filter(year <= 2018)%>%
  mutate(art_coverage = if_else(year < 2004, 0,#before 2004 ART coverage did not exist
                                #calibration period scale up to art coverage
                                if_else(year<=2017 & sex == 'Male', 
                                        ((year-2003)*scale_up_male),
                                        if_else(year<=2017 & sex == 'Female',
                                                ((year-2003)*scale_up_female),
                                                if_else(year==2018 & sex == 'Male', 
                                                        ((n_years_warmup_with_art)*scale_up_male),
                                                        ((n_years_warmup_with_art)*scale_up_female))))))

art_coverage_df$program_id<-rep(1, times = nrow(art_coverage_df))


#add in information from other policies

#p2 males
policy_2_art_coverage_male<-art_scenario_ref_df%>%
  filter(program_id == 2,
         G_SET == 1)

art_coverage_df<-rbind(art_coverage_df, 
                       data.frame(year = 2018, sex = 'Male', 
                                  art_coverage = policy_2_art_coverage_male$art_coverage,
                                  program_id = 2))

#p2 females
policy_2_art_coverage_female<-art_scenario_ref_df%>%
  filter(program_id == 2,
         G_SET == 2)

art_coverage_df<-rbind(art_coverage_df, 
                       data.frame(year = 2018, sex = 'Female', 
                                  art_coverage = policy_2_art_coverage_female$art_coverage,
                                  program_id = 2))

#p3 males
policy_3_art_coverage_male<-art_scenario_ref_df%>%
  filter(program_id == 3,
         G_SET == 1)

art_coverage_df<-rbind(art_coverage_df, 
                       data.frame(year = 2018, sex = 'Male', 
                                  art_coverage = policy_3_art_coverage_male$art_coverage,
                                  program_id = 3))

#p3 females
policy_3_art_coverage_female<-art_scenario_ref_df%>%
  filter(program_id == 3,
         G_SET == 2)

art_coverage_df<-rbind(art_coverage_df, 
                       data.frame(year = 2018, sex = 'Female', 
                                  art_coverage = policy_2_art_coverage_female$art_coverage,
                                  program_id = 3))


setwd(outdir)
write.csv(art_coverage_df, 'art_coverage_df.csv', row.names = FALSE)

