#Estimates populations every year by gender in KZN
#Uses GBD

#for 15 - 59 (all ages in model) used for mortality and hiv incidnece estimate
#15-19 (for calculating aging into model, births perc)

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa
indir_outdir <- paste0(here(), '/param_files/calculated_param_gen/input_data/GBD/')

setwd(indir_outdir)

#get pop estimates by taking 1/rate * num of deaths (for all cause)
pop_df<-read.csv('all_cause_mort_num_rate_pulled_2_22_22.csv')%>%
  select(-c('location', 'cause', 'upper', 'lower', 'measure'))
pop_df<-dcast(pop_df, sex+age+year ~ metric, mean, 
              value.var = 'val') #note base population estimates are the same 

##all ages in model
pop_df_all_ages<-pop_df%>%
  mutate(pop_val = (Number/Rate))%>%
  group_by(year, sex)%>%
  summarise(expected_total_pop = sum(pop_val)*100000)

write.csv(pop_df_all_ages, 'pop_estimates_all_ages.csv', row.names = FALSE)

##15-19
pop_df_15_19<-pop_df%>%
  filter(age == "15 to 19")%>%
  mutate(pop_val = (Number/Rate))%>%
  group_by(year, sex)%>%
  summarise(expected_total_pop = sum(pop_val)*100000)

write.csv(pop_df_15_19, 'pop_estimates_15_19.csv', row.names = FALSE)

