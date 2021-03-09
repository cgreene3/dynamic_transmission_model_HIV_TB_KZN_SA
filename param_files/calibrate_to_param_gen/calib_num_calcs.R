#Estimating active TB mortalities per 100K to calibrate to
#Uses GBD 2019 - Number of Deaths and Scales to 100K
#SA population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/calibrate_to_param_gen')
indir_pop <- paste0(here(),'/param_files/GBD_pop_estimates')
outdir <- paste0(here(),'/param_files')

#to name plots according to date generated
data_gen_date<-Sys.Date()
data_gen_date<-str_replace_all(data_gen_date, '-', '_')

##########read in pop files#################
#data pull date 02/28/2021
#put GBD estimates into seperate folder so can read in all

#commented to just read in pop_df for quicker read - rerun in need updated GBD pop estimates

setwd(indir_pop) 
#pop_files = list.files(pattern="*.CSV")

#pop_df<-data.frame()

#lapply(pop_files, function(file){
#  temp_df <-read.csv(file)%>%
#    filter(location_id == 196,
#           age_group_id >= 8, #between 15
#           age_group_id <=15, #less than 60
#           sex_id != 3) #only females and males
#  pop_df <<- rbind(pop_df, temp_df)
#})

#group populations over all age groups
#pop_df<-pop_df%>%
#  group_by(sex_id, sex_name, year_id)%>%
#  summarise(expected_total_pop = sum(val))%>%
#  rename(year = year_id)

#write_csv(pop_df, 'pop_df.csv')

#read pop
pop_df<-read.csv('pop_df.csv')

#group pop df by gender since rate will be per 100k people
pop_df<-pop_df%>%
  group_by(year)%>%
  summarise(expected_total_pop = sum(expected_total_pop))

#read in TB mort estimates
setwd(indir)
total_mort_df<-read.csv('IHME_TB_deaths_1990_2017_Mar4.csv')
total_mort_df<-total_mort_df%>%
  filter(cause_id != 300)%>%
  mutate(co_infection = if_else(cause_id %in% c(946, 947, 934), 'TB_only', 'HIV/TB_coinfection'),
         calibration_group = paste0(co_infection, '_', sex_name))

total_mort_calibration_params<-total_mort_df%>%
  group_by(year, sex_name, calibration_group)%>%
  summarise(min_calib_total = sum(lower),
            max_calib_total = sum(upper),
            expected_calib_total = sum(val))%>%
  mutate(sex_name = tolower(sex_name))%>%
  left_join(pop_df, by = c('year'))%>%
  mutate(min_percent = min_calib_total/expected_total_pop,
         max_percent =max_calib_total/expected_total_pop,
         expected_percent =expected_calib_total/expected_total_pop,
         min_rate = min_percent*100000,
         max_rate = max_percent*100000,
         expected_rate = expected_percent*100000)%>%
  select(c('year', 'sex_name', 'calibration_group', 'min_rate', 'max_rate', 'expected_rate'))

setwd(outdir)
write.csv(total_mort_calibration_params, 'calibration_rates_df.csv', row.names = FALSE)

