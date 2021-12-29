#Calibration Analysis Dec 21th

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops', 'data.table'), require, character.only=T)

#############Set in directory and out directory###########

start_eval_date <- '2021_12_21' #set to date started calibration
calib_analysis_desc<- "" #'_W_M'

###########Make sure epi_model_HIV_TB/SouthAfrica_KZN.Rproj is open, otherwise will need to change wd manually########
indir_ref_data<-paste0(here(),'/param_files/')

#set calib outputs to local file on computer (otherwise files will overwhelm github)
indir_state_prog<-paste0("~/Documents/academic_posttt_2020/HIV_TB/calib_outdir/", start_eval_date)
outdir_mse_dfs<-paste0(here(),'/param_files/calibration_code_results/', start_eval_date, calib_analysis_desc)
outdir_analysis <- paste0(here(),'/param_files/calibration_code_results/', start_eval_date, 
                          calib_analysis_desc, '/analysis')

#######read in calibration ref df#########
setwd(indir_ref_data)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')

# #####Graphs that describe model states overtime####
HIV_calibration_df<-read.csv('hiv_prev_calibration_df.csv')

ART_df<-read.csv('art_initiation_related_param_df.csv')%>%
  filter(POLICY_ID == 1)

#GBD calibration rates in general param file folder
GBD_calibration_df<-read.csv('TB_mort_calibration_rates_df.csv')
#match calibration group with model output calibration groups
GBD_calibration_df$calibration_group<-if_else(GBD_calibration_df$calibration_group == "TB_only_Female",
                                              "HIV_neg_female",
                                              if_else(GBD_calibration_df$calibration_group == "TB_only_Male",
                                                      "HIV_neg_male",
                                                      if_else(GBD_calibration_df$calibration_group == "HIV/TB_coinfection_Male",
                                                              "HIV_pos_male", "HIV_pos_female")))


GBD_calibration_df$G_compartment <-if_else(GBD_calibration_df$sex_name == 'male',
                                           1, 2)

#read in mse df top selected from analysis
setwd(outdir_mse_dfs)
mse_df_top<-read.csv('mse_df_top.csv')

########pull incidence and mort stats########
setwd(indir_state_prog)

top_df_all_info_df<-data.frame()

lapply(unique(mse_df_top$sim_id), function(sim_id_itr){ 
  file_name<-paste0('out_df_sim_id_', sim_id_itr, '.csv')
  print(sim_id_itr)
  temp_data <- fread(file_name) 
  temp_data$sim_id <- rep(sim_id_itr, times = nrow(temp_data))
  #for each iteration, bind the new data to the building dataset
  top_df_all_info_df <<- rbind(top_df_all_info_df, temp_data) 
})

state_prog_df <-reshape2::melt(top_df_all_info_df, 
                               id.vars = c("time", "year",
                                           'sim_id'))

state_prog_df <- cbind(state_prog_df,
                       data.frame(do.call('rbind',
                                          strsplit(as.character(state_prog_df$variable),
                                                   '_',fixed=TRUE))))

names(state_prog_df)[names(state_prog_df) == "X2"] <- "TB_compartment"
names(state_prog_df)[names(state_prog_df) == "X3"] <- "DR_compartment"
names(state_prog_df)[names(state_prog_df) == "X4"] <- "HIV_compartment"
names(state_prog_df)[names(state_prog_df) == "X5"] <- "G_compartment"
state_prog_df<-subset(state_prog_df, select = -X1)
state_prog_df<-subset(state_prog_df, select = -variable)

#pull incidence stats from 2017
incidence_stats_df_4cats<-state_prog_df%>%
  filter(TB_compartment == 'incidence')%>%
  mutate(calibration_group = paste0('HIV_', HIV_compartment, '_', G_compartment))%>%
  select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
  group_by(year, sim_id, calibration_group)%>%
  summarise(cum_incidence = max(value))%>%
  as.data.frame()%>%
  group_by(calibration_group, sim_id)%>%
  mutate(incidence_in_year = cum_incidence - lag(cum_incidence))%>%
  ungroup()%>%
  mutate(incidence_in_year = if_else(year == 1990, cum_incidence, incidence_in_year))%>%
  filter(year == 2017)%>%
  group_by(calibration_group)%>%
  summarise(max_incidence = max(incidence_in_year),
            min_incidence = min(incidence_in_year),
            mean_incidence = mean(incidence_in_year))

incidence_stats_df<-state_prog_df%>%
  filter(TB_compartment == 'incidence')%>%
  mutate(calibration_group = paste0('HIV_', HIV_compartment, '_', G_compartment))%>%
  select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
  group_by(year, sim_id, calibration_group)%>%
  summarise(cum_incidence = max(value))%>%
  as.data.frame()%>%
  group_by(calibration_group, sim_id)%>%
  mutate(incidence_in_year = cum_incidence - lag(cum_incidence))%>%
  ungroup()%>%
  mutate(incidence_in_year = if_else(year == 1990, cum_incidence, incidence_in_year))%>%
  filter(year == 2017)%>%
  group_by(sim_id)%>%
  summarise(incidence_in_year_total = sum(incidence_in_year))%>%
  summarise(max_incidence = max(incidence_in_year_total),
            min_incidence = min(incidence_in_year_total),
            mean_incidence = mean(incidence_in_year_total))

setwd(outdir_mse_dfs)
write.csv(incidence_stats_df, 'incidence_stats_df.csv', row.names = FALSE)
write.csv(incidence_stats_df_4cats, 'incidence_stats_df_4cats.csv', row.names = FALSE)

#pull prev stats from 2017
prev_stats_df<-state_prog_df%>%
  filter(TB_compartment == 6)%>%
  group_by(time, year, sim_id)%>%
  summarise(prev = sum(value))%>%
  group_by(year, sim_id)%>%
  summarise(prev = mean(prev))%>%
  group_by(year)%>%
  summarise(max_prev = max(prev),
            min_prev = min(prev),
            mean_prev = mean(prev))%>%
  filter(year == 2017)

write.csv(prev_stats_df, 'prev_stats_df.csv', row.names = FALSE)

#pull mort stats from 2017
mort_stats_df_4cats<-state_prog_df%>%
  filter(TB_compartment == 'mort')%>%
  mutate(calibration_group = paste0('HIV_', HIV_compartment, '_', G_compartment))%>%
  select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
  group_by(year, sim_id, calibration_group)%>%
  summarise(cum_mort = max(value))%>%
  as.data.frame()%>%
  group_by(calibration_group, sim_id)%>%
  mutate(mort_in_year = cum_mort - lag(cum_mort))%>%
  ungroup()%>%
  mutate(mort_in_year = if_else(year == 1990, cum_mort, mort_in_year))%>%
  filter(year == 2017)%>%
  group_by(calibration_group)%>%
  summarise(max_mort_model = max(mort_in_year),
            min_mort_model = min(mort_in_year),
            mean_mort_model = mean(mort_in_year))

mort_stats_df<-state_prog_df%>%
  filter(TB_compartment == 'mort')%>%
  mutate(calibration_group = paste0('HIV_', HIV_compartment, '_', G_compartment))%>%
  select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
  group_by(year, sim_id, calibration_group)%>%
  summarise(cum_mort = max(value))%>%
  as.data.frame()%>%
  group_by(calibration_group, sim_id)%>%
  mutate(mort_in_year = cum_mort - lag(cum_mort))%>%
  ungroup()%>%
  mutate(mort_in_year = if_else(year == 1990, cum_mort, mort_in_year))%>%
  filter(year == 2017)%>%
  group_by(sim_id)%>%
  summarise(mort_in_year_total = sum(mort_in_year))%>%
  summarise(max_mort = max(mort_in_year_total),
            min_mort = min(mort_in_year_total),
            mean_mort = mean(mort_in_year_total))

write.csv(mort_stats_df, 'mort_stats_df.csv', row.names = FALSE)
write.csv(mort_stats_df_4cats, 'mort_stats_df_4cats.csv', row.names = FALSE)


#compare mort per 100K males/females (as in graphs)...above as per 100K people
total_pop_in_gender_df<-state_prog_df%>%
  filter(G_compartment %in% c(1,2))%>%
  group_by(time, year, G_compartment, sim_id)%>%
  summarise(total_gender_pop = sum(value))%>%
  ungroup()%>%
  group_by(year, sim_id, G_compartment)%>%
  summarise(total_gender_pop = median(total_gender_pop))

total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)

mort_calibration_df<-state_prog_df%>%
  filter(TB_compartment == 'mort')%>%
  mutate(calibration_group = paste0('HIV_', HIV_compartment, '_', G_compartment))%>%
  select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
  group_by(year, sim_id, calibration_group)%>%
  summarise(cum_mort = max(value))%>%
  as.data.frame()%>%
  group_by(calibration_group, sim_id)%>%
  mutate(mort_in_year = cum_mort - lag(cum_mort))%>%
  ungroup()%>%
  mutate(mort_in_year = if_else(year == 1990, cum_mort, mort_in_year))%>%
  filter(year == 2017)%>%
  mutate(G_compartment = if_else(grepl('female', calibration_group), '2', '1'))%>%
  left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
  mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
  as.data.frame()%>%
  group_by(calibration_group)%>%
  summarise(max_mort_model_per100K_gender = max(mort_rate_per_100K),
            min_mort_model_per100K_gender = min(mort_rate_per_100K),
            mean_mort_model_per100K_gender = mean(mort_rate_per_100K))%>%
  left_join(GBD_calibration_df%>%
              filter(year == 2017)%>%
            select(min_rate, max_rate, expected_rate, calibration_group),
            by = c('calibration_group'))

write.csv(mort_calibration_df, 'mort_stats_df_4cats_per100K_gender_GBD_compare.csv', row.names = FALSE)

