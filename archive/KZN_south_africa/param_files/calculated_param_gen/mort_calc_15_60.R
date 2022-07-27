#Estimating background mortality rate from GBD
#Uses GBD 2019 - Number of Deaths
#KZN population estimates
#updated Jan 20th

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/calculated_param_gen/input_data')
outdir <- paste0(here(),'/param_files')

##########read in pop files#################
#data pull date 02/28/2021

setwd(indir) 
rate_df<-read.csv('KZN_rate_IHME_deaths_1990_2017_Jan20.csv')
num_df<-read.csv('KZN_num_IHME_deaths_1990_2017_Jan20.csv')

#get pop estimates by taking 1/rate * num of deaths (for all cause)
pop_df<-num_df%>%
  filter(cause_id == 294)%>%
  left_join(rate_df, by = c('age_id', 'sex_id', 'year'))%>%
  select('age_id', 'sex_id', 'year', 'val.x', 'val.y')%>% #val.x = num, val.y = rate
  mutate(pop_val = (val.x * (1/val.y))*100000)

#group populations over all age groups
pop_df<-pop_df%>%
  group_by(sex_id, year)%>%
  summarise(expected_total_pop = sum(pop_val))

setwd(indir)
write.csv(pop_df, 'pop_df.csv', row.names = FALSE)
#pop_df<-read.csv('pop_df.csv')

#prevalence and mortality data frames consider the following disease states
#HIV/AIDS - causeid 298
#DS TB - cause id 934
#MDR TB - cause ids 946, 947

##########read in mortality data##########

#pulls all cause, TB, and HIV deaths
#subtracts out TB and HIV related deaths
mort_num_df<-num_df%>%
  mutate(disease_cat = if_else(cause_name == 'All causes', 'base', 'disease'))%>%
  group_by(sex_id, year, disease_cat)%>%
  summarise(expected_disease_mort = sum(val))%>%
  ungroup()%>%
  mutate(base_mort_calcs = if_else(disease_cat == 'base', 
                                   expected_disease_mort, 
                                   -expected_disease_mort))%>% 
  group_by(sex_id, year)%>%
  summarise(expected_non_disease_mort_total = sum(base_mort_calcs))


mort_rate_df<-pop_df%>%
  left_join(mort_num_df, by = c('sex_id', 'year'))%>%
  mutate(mort_rate = expected_non_disease_mort_total/expected_total_pop,
         mort_rate_per_100K = mort_rate*100000)


setwd(outdir)
write.csv(mort_rate_df, 'KZN_mort_rate_df.csv', row.names = FALSE)

