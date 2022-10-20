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
pop_df_15_59<-read.csv('all_cause_mort_num_rate.csv')%>%
  select(-c('location', 'cause', 'upper', 'lower', 'measure'))
pop_df_15_19<-read.csv('all_cause_mort_num_rate.csv')%>%
  select(-c('location', 'cause', 'upper', 'lower', 'measure'))%>%
  filter(age == "15-19 years")
pop_df_15_59<-dcast(pop_df_15_59, sex+age+year ~ metric, mean, 
              value.var = 'val') #note base population estimates are the same
pop_df_15_19<-dcast(pop_df_15_19, sex+age+year ~ metric, mean, 
              value.var = 'val') #note base population estimates are the same

##all ages in model
pop_df_15_59<-pop_df_15_59%>%
  mutate(Rate = Rate/100000)%>% #is per 100,000
  mutate(pop_val = Number/Rate)%>%
  group_by(year, sex)%>%
  summarise(expected_total_pop = sum(pop_val))

write.csv(pop_df_15_59, 'pop_estimates_15_59.csv', row.names = FALSE)

##15_19

pop_df_15_19<-pop_df_15_19%>%
  mutate(pop_val = (Number/Rate))%>%
  group_by(year, sex)%>%
  summarise(expected_total_pop = sum(pop_val)*100000)

write.csv(pop_df_15_19, 'pop_estimates_15_19.csv', row.names = FALSE)

