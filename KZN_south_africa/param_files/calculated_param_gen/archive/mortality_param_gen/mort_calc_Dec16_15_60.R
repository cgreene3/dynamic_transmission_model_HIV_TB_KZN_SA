#Estimating background mortality rate from GBD
#Uses GBD 2019 - Number of Deaths
#SA population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#to name plots according to date generated
data_gen_date<-'2021_12_16'
#data_gen_date<-str_replace_all(data_gen_date, '-', '_')

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/mortality_param_gen')
indir_pop <- paste0(here(),'/param_files/GBD_pop_estimates')
outdir <- paste0(here(),'/param_files/calibration_code_results/', data_gen_date)
#will give warning if wd already exists
#ignore warning!
#dir.create(file.path(outdir))

##########read in pop files#################
#data pull date 02/28/2021
#put GBD estimates into seperate folder so can read in all

#commented to just read in pop_df for quicker read - rerun in need updated GBD pop estimates

setwd(indir_pop) 
pop_files = list.files(pattern="*.CSV")

pop_df<-data.frame()

lapply(pop_files, function(file){
  temp_df <-read.csv(file)%>%
    filter(location_id == 196,
           age_group_id >= 8, #between 15
           age_group_id <=16, #less than 60
           sex_id != 3) #only females and males
  pop_df <<- rbind(pop_df, temp_df)
})

#group populations over all age groups
pop_df<-pop_df%>%
  group_by(sex_id, sex_name, year_id)%>%
  summarise(expected_total_pop = sum(val))%>%
  rename(year = year_id)

write.csv(pop_df, 'pop_df.csv')
#pop_df<-read.csv('pop_df.csv')

#prevalence and mortality data frames consider the following disease states
#HIV/AIDS - causeid 298
#DS TB - cause id 934
#MDR TB - cause ids 946, 947

#########read in prevalence data#############
setwd(indir)

#read in population data
#prev_df<-read.csv('IHME_prev_1990_2017_Feb28.csv')
#prev_df<-prev_df%>%
#  group_by(sex_id, year)%>%
#  summarise(expected_disease_pop = sum(val))
  
#use prevalence data to calculate non disease mort
#non_disease_mort_df<-pop_df%>%
#  left_join(prev_df, by = c('sex_id', 'year'))%>%
#  mutate(expected_non_disease_pop = expected_total_pop-expected_disease_pop)


###changed from -expected_disease_mort that subtracts out disease pop

##########read in mortality data##########
mort_num_df<-read.csv('IHME_deaths_1990_2017_Feb28.csv')%>%
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
  left_join(mort_df, by = c('sex_id', 'year'))%>%
  mutate(mort_rate = expected_non_disease_mort_total/expected_total_pop,
         mort_rate_per_100K = mort_rate*100000)


setwd(outdir)
write.csv(mort_rate_df, 'mort_rate_df.csv', row.names = FALSE)

