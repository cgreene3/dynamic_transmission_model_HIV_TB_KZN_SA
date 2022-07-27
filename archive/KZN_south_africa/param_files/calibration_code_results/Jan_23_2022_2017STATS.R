#Calibration Analysis Dec 21th

#clean workspace
rm(list = ls())
gc()

#############Set in directory and out directory###########
start_eval_date <- '2022_01_21' #set to date started calibration
calib_analysis_desc<- "" #'_W_M'

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops', 'data.table'), require, character.only=T)

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
GBD_calibration_df<-read.csv('KZN_TB_mort_calibration_rates_df.csv')
#match calibration group with model output calibration groups
GBD_calibration_df$calibration_group<-if_else(GBD_calibration_df$calibration_group == "TB_only_Female",
                                              "HIV_neg_female",
                                              if_else(GBD_calibration_df$calibration_group == "TB_only_Male",
                                                      "HIV_neg_male",
                                                      if_else(GBD_calibration_df$calibration_group == "HIV/TB_coinfection_Male",
                                                              "HIV_pos_male", "HIV_pos_female")))


names(GBD_calibration_df)[names(GBD_calibration_df) == 'sex_id'] <- "G_compartment"

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

#compare rates per 100K males/females
total_pop_in_gender_df<-state_prog_df%>%
  filter(G_compartment %in% c(1,2))%>%
  group_by(time, year, G_compartment, sim_id)%>%
  summarise(total_gender_pop = sum(value))%>%
  ungroup()%>%
  group_by(year, sim_id, G_compartment)%>%
  summarise(total_gender_pop = median(total_gender_pop))%>%
  mutate(gender_name = if_else(G_compartment == 1, 'male', 'female'))

total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)

#######pull incidence stats from 2017##########

#for each category (males HIV-, females HIV-, males HIV+, females HIV+)
#by sim
incidence_stats_df_4cats_by_sim<-state_prog_df%>%
  filter(TB_compartment == 'incidence')%>%
  mutate(calibration_group = paste0('HIV_', HIV_compartment, '_', G_compartment))%>%
  select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
  group_by(year, sim_id, calibration_group)%>%
  summarise(cum_incidence = max(value))%>%
  as.data.frame()%>%
  group_by(calibration_group, sim_id)%>%
  mutate(incidence_in_year = cum_incidence - lag(cum_incidence))%>%
  ungroup()%>%
  mutate(incidence_in_year_per_100K_people = if_else(year == 1990, cum_incidence, incidence_in_year))%>%
  filter(year == 2017)%>%
  mutate(gender_name = if_else(grepl('female', calibration_group), 'female', 'male'))%>%
  left_join(total_pop_in_gender_df, by = c('gender_name', 'sim_id', 'year'))%>%
  mutate(incidence_in_year_per_100K_gender = (incidence_in_year_per_100K_people/total_gender_pop)*100000)


#overall
#get max/min sims
incidence_stats_overall_df_max_min<-incidence_stats_df_4cats_by_sim%>%
  group_by(sim_id)%>%
  summarise(incidence_in_year_per_100K_people = sum(incidence_in_year_per_100K_people))%>%
  slice(which.max(incidence_in_year_per_100K_people),
        which.min(incidence_in_year_per_100K_people))%>%
  mutate(range_indicator = if_else(rank(incidence_in_year_per_100K_people)==1, 'min', 'max'))

incidence_stats_df_4cats_overall<-incidence_stats_df_4cats_by_sim%>%
  filter(sim_id %in% unique(incidence_stats_overall_df_max_min$sim_id))%>%
  left_join(incidence_stats_overall_df_max_min%>%select(-c('incidence_in_year_per_100K_people')), by = c('sim_id'))%>%
  select(c('calibration_group', 'incidence_in_year_per_100K_people', 'range_indicator'))%>%
  rbind(., incidence_stats_df_4cats_by_sim%>%
          group_by(calibration_group)%>%
          summarise(incidence_in_year_per_100K_people = mean(incidence_in_year_per_100K_people))%>%
          mutate(range_indicator = 'mean'))

incidence_stats_overall_df<-incidence_stats_df_4cats_overall%>%
  group_by(range_indicator)%>%
  summarise(incidence_in_year_per_100K_people = sum(incidence_in_year_per_100K_people))
  
#per 100K, gender
#get max/min sims
incidence_stats_gender_df_max_min<-incidence_stats_df_4cats_by_sim%>%
  filter(gender_name == 'male')%>%
  group_by(sim_id)%>%
  summarise(incidence_in_year_per_100K_gender = sum(incidence_in_year_per_100K_gender))%>%
  slice(which.max(incidence_in_year_per_100K_gender),
        which.min(incidence_in_year_per_100K_gender))%>%
  mutate(gender_name = 'male',
         range_indicator = if_else(rank(incidence_in_year_per_100K_gender)==1, 'min', 'max'))%>%
  rbind(., incidence_stats_df_4cats_by_sim%>%
          filter(gender_name == 'female')%>%
          group_by(sim_id)%>%
          summarise(incidence_in_year_per_100K_gender = sum(incidence_in_year_per_100K_gender))%>%
          slice(which.max(incidence_in_year_per_100K_gender),
                which.min(incidence_in_year_per_100K_gender))%>%
          mutate(gender_name = 'female',
                 range_indicator = if_else(rank(incidence_in_year_per_100K_gender)==1, 'min', 'max')))%>%
  mutate(sim_id_gender_name = paste0(sim_id, gender_name))

incidence_stats_df_4cats_gender<-incidence_stats_df_4cats_by_sim%>%
  mutate(sim_id_gender_name = paste0(sim_id, gender_name))%>%
  filter(sim_id_gender_name %in% unique(incidence_stats_gender_df_max_min$sim_id_gender_name))%>%
  left_join(incidence_stats_gender_df_max_min%>%select(-c('incidence_in_year_per_100K_gender')), by = c('sim_id_gender_name', 'gender_name'))%>%
  select(c('calibration_group', 'gender_name', 'incidence_in_year_per_100K_gender', 'range_indicator'))%>%
  rbind(., incidence_stats_df_4cats_by_sim%>%
          group_by(calibration_group)%>%
          summarise(incidence_in_year_per_100K_gender = mean(incidence_in_year_per_100K_gender))%>%
          mutate(range_indicator = 'mean',
                 gender_name = if_else(grepl('female', calibration_group), 'female', 'male')))

incidence_stats_gender_df<-incidence_stats_df_4cats_gender%>%
  group_by(range_indicator, gender_name)%>%
  summarise(incidence_in_year_per_100K_gender = sum(incidence_in_year_per_100K_gender))

# 
# incidence_stats_sentence<-paste0("In the baseline scenario, in the year prior to the start of the intervention in 2017, ", 
# " we estimated an annual incidence of active TB disease of ", 
# prettyNum(incidence_stats_df$mean_incidence_per_100K_people[1], big.mark = ",", scientific = FALSE),
# " [", prettyNum(incidence_stats_df$min_incidence_per_100K_people[1], big.mark = ",", scientific = FALSE),
# "; ", prettyNum(incidence_stats_df$max_incidence_per_100K_people[1], big.mark = ",", scientific = FALSE), "] ",
# "including ", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_pos_female")%>%select(mean_incidence_per_100K_people),
#                         big.mark = ",", scientific = FALSE),
# " [", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_pos_female")%>%select(min_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE),
# "; ", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_pos_female")%>%select(max_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE), "] ",
# "among women living with HIV, ",
# prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_neg_female")%>%select(mean_incidence_per_100K_people),
#           big.mark = ",", scientific = FALSE),
# " [", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_neg_female")%>%select(min_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE),
# "; ", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_neg_female")%>%select(max_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE), "] ",
# "among women without HIV, ", 
# prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_pos_male")%>%select(mean_incidence_per_100K_people),
#                         big.mark = ",", scientific = FALSE),
# " [", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_pos_male")%>%select(min_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE),
# "; ", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_pos_male")%>%select(max_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE), "] ",
# "among men living with HIV, ",
# prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_neg_male")%>%select(mean_incidence_per_100K_people),
#           big.mark = ",", scientific = FALSE),
# " [", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_neg_male")%>%select(min_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE),
# "; ", prettyNum(incidence_stats_df_4cats%>%filter(calibration_group == "HIV_neg_male")%>%select(max_incidence_per_100K_people),
#                 big.mark = ",", scientific = FALSE), "] ",
# "among men without HIV.")

setwd(outdir_mse_dfs)

#make easier comparison to GBD
GBD_compare_df<-data.frame(matrix(ncol = length(c('metric', 'gender', incidence_stats_overall_df$range_indicator)),
                        nrow = 0))
colnames(GBD_compare_df) <- c('metric', 'gender', incidence_stats_overall_df$range_indicator)
GBD_compare_df[nrow(GBD_compare_df) + 1,] = c('incidence', 'both', 
                          round(incidence_stats_overall_df$incidence_in_year_per_100K_people))
temp<-incidence_stats_gender_df%>%
  filter(gender_name == 'female')
GBD_compare_df[nrow(GBD_compare_df) + 1,] = c('incidence', 'female', 
                          round(temp$incidence_in_year_per_100K_gender))
temp<-incidence_stats_gender_df%>%
  filter(gender_name == 'male')
GBD_compare_df[nrow(GBD_compare_df) + 1,] = c('incidence', 'male', 
                          round(temp$incidence_in_year_per_100K_gender))

GBD_compare_df$category<-rep('all_categories', times = nrow(GBD_compare_df))

GBD_compare_df<-GBD_compare_df%>%select(c('metric', 'gender', 'category', 'mean', 'max', 'min'))

incidence_stats_df_4cats_overall$gender = 'both'
incidence_stats_df_4cats_overall$metric = 'incidence'
names(incidence_stats_df_4cats_overall)[names(incidence_stats_df_4cats_overall) 
                                        == 'calibration_group']<-'category'
incidence_stats_df_4cats_overall$incidence_in_year_per_100K_people<-round(incidence_stats_df_4cats_overall$incidence_in_year_per_100K_people)

temp<-reshape2::dcast(incidence_stats_df_4cats_overall, 
            metric + gender + category ~ range_indicator, value.var = "incidence_in_year_per_100K_people")%>%
  select(c('metric', 'gender', 'category', 'mean', 'max', 'min'))

GBD_compare_df<-rbind(GBD_compare_df, temp)

names(incidence_stats_df_4cats_gender)[names(incidence_stats_df_4cats_gender) 
                                       == 'gender_name']<-'gender'
incidence_stats_df_4cats_gender$metric = 'incidence'
names(incidence_stats_df_4cats_gender)[names(incidence_stats_df_4cats_gender) 
                                        == 'calibration_group']<-'category'
incidence_stats_df_4cats_gender$incidence_in_year_per_100K_gender<-
  round(incidence_stats_df_4cats_gender$incidence_in_year_per_100K_gender)

temp<-reshape2::dcast(incidence_stats_df_4cats_gender, 
                      metric + gender + category ~ range_indicator, value.var = "incidence_in_year_per_100K_gender")%>%
  select(c('metric', 'gender', 'category', 'mean', 'max', 'min'))

GBD_compare_df<-rbind(GBD_compare_df, temp)

# write.csv(incidence_stats_overall_df, 'incidence_stats_overall_df.csv', row.names = FALSE)
# write.csv(incidence_stats_df_4cats_overall, 'incidence_stats_df_4cats_overall.csv', row.names = FALSE)
# write.csv(incidence_stats_gender_df, 'incidence_stats_gender_df.csv', row.names = FALSE)
# write.csv(incidence_stats_df_4cats_gender, 'incidence_stats_df_4cats_gender.csv', row.names = FALSE)

#if want to automatically write sentence to file
# write.table(incidence_stats_sentence, file = "incidence_stats_sentence.txt",
#             col.names = FALSE, row.names = FALSE)


#######pull prevalence stats from 2017##########

prev_stats_df_4_cats_by_sim_df<-state_prog_df%>%
  filter(TB_compartment == 6,
         year == 2017,
         time %% 1 == 0.5)%>% #take prev values from middle of the year
  mutate(group = paste0(if_else(HIV_compartment == 1, "HIV_neg", "HIV_pos"),
                        if_else(G_compartment == 1, "_male", "_female")))%>%
  group_by(sim_id, group)%>%
  summarise(prevalence_per_100K_people = sum(value))%>%
  mutate(gender_name = if_else(grepl('female', group), 'female', 'male'))%>%
  left_join(total_pop_in_gender_df%>%filter(year == 2017), by = c('sim_id', 'gender_name'))%>%
  mutate(prevalence_per_100K_gender = (prevalence_per_100K_people/total_gender_pop)*100000)

#per 100K, overall
#get max/min sims
prev_stats_overall_max_min <-prev_stats_df_4_cats_by_sim_df%>%
  group_by(sim_id)%>%
  summarise(prevalence_per_100K_people = sum(prevalence_per_100K_people))%>%
  slice(which.max(prevalence_per_100K_people),
        which.min(prevalence_per_100K_people))%>%
  mutate(range_indicator = if_else(rank(prevalence_per_100K_people)==1, 'min', 'max'))

prev_stats_df_4cats_overall<-prev_stats_df_4_cats_by_sim_df%>%
  filter(sim_id %in% unique(prev_stats_overall_max_min$sim_id))%>%
  left_join(prev_stats_overall_max_min%>%select(-c('prevalence_per_100K_people')), by = c('sim_id'))%>%
  select(c('group', 'prevalence_per_100K_people', 'range_indicator'))%>%
  rbind(., prev_stats_df_4_cats_by_sim_df%>%
          group_by(group)%>%
          summarise(prevalence_per_100K_people = mean(prevalence_per_100K_people))%>%
          mutate(range_indicator = 'mean'))

prev_stats_overall_df<-prev_stats_df_4cats_overall%>%
  group_by(range_indicator)%>%
  summarise(prevalence_per_100K_people = sum(prevalence_per_100K_people))

#per 100K, gender
#get max/min sims
prev_stats_gender_df_max_min<-prev_stats_df_4_cats_by_sim_df%>%
  filter(gender_name == 'male')%>%
  group_by(sim_id)%>%
  summarise(prevalence_per_100K_gender = sum(prevalence_per_100K_gender))%>%
  slice(which.max(prevalence_per_100K_gender),
        which.min(prevalence_per_100K_gender))%>%
  mutate(gender_name = 'male',
         range_indicator = if_else(rank(prevalence_per_100K_gender)==1, 'min', 'max'))%>%
  rbind(., prev_stats_df_4_cats_by_sim_df%>%
          filter(gender_name == 'female')%>%
          group_by(sim_id)%>%
          summarise(prevalence_per_100K_gender = sum(prevalence_per_100K_gender))%>%
          slice(which.max(prevalence_per_100K_gender),
                which.min(prevalence_per_100K_gender))%>%
          mutate(gender_name = 'female',
                 range_indicator = if_else(rank(prevalence_per_100K_gender)==1, 'min', 'max')))%>%
  mutate(sim_id_gender_name = paste0(sim_id, gender_name))

prev_stats_df_4cats_gender<-prev_stats_df_4_cats_by_sim_df%>%
  mutate(sim_id_gender_name = paste0(sim_id, gender_name))%>%
  filter(sim_id_gender_name %in% unique(prev_stats_gender_df_max_min$sim_id_gender_name))%>%
  left_join(prev_stats_gender_df_max_min%>%select(-c('prevalence_per_100K_gender')), by = c('sim_id_gender_name', 'gender_name'))%>%
  select(c('group', 'gender_name', 'prevalence_per_100K_gender', 'range_indicator'))%>%
  rbind(., prev_stats_df_4_cats_by_sim_df%>%
          group_by(group)%>%
          summarise(prevalence_per_100K_gender = mean(prevalence_per_100K_gender))%>%
          mutate(range_indicator = 'mean',
                 gender_name = if_else(grepl('female', group), 'female', 'male')))

prev_stats_gender_df<-prev_stats_df_4cats_gender%>%
  group_by(range_indicator, gender_name)%>%
  summarise(prevalence_per_100K_gender = sum(prevalence_per_100K_gender))

setwd(outdir_mse_dfs)

#make easier comparison to GBD
temp<-prev_stats_overall_df%>%
  mutate(order_fun = if_else(range_indicator == 'mean', 1, if_else(range_indicator == 'max', 2, 3)))%>%
  arrange(order_fun)

GBD_compare_df[nrow(GBD_compare_df) + 1,] = c('prevalence', 'both', 'all_categories',
                                              round(temp$prevalence_per_100K_people))
temp<-prev_stats_gender_df%>%
  filter(gender_name == 'female')%>%
  mutate(order_fun = if_else(range_indicator == 'mean', 1, if_else(range_indicator == 'max', 2, 3)))%>%
  arrange(order_fun)
  
GBD_compare_df[nrow(GBD_compare_df) + 1,] = c('prevalence', 'female', 'all_categories',
                                              round(temp$prevalence_per_100K_gender))
temp<-prev_stats_gender_df%>%
  filter(gender_name == 'male')%>%
  mutate(order_fun = if_else(range_indicator == 'mean', 1, if_else(range_indicator == 'max', 2, 3)))%>%
  arrange(order_fun)

GBD_compare_df[nrow(GBD_compare_df) + 1,] = c('prevalence', 'male', 'all_categories',
                                              round(temp$prevalence_per_100K_gender))


prev_stats_df_4cats_overall$gender = 'both'
prev_stats_df_4cats_overall$metric = 'prevalence'
names(prev_stats_df_4cats_overall)[names(prev_stats_df_4cats_overall) 
                                        == 'group']<-'category'
prev_stats_df_4cats_overall$prevalence_per_100K_people<-round(prev_stats_df_4cats_overall$prevalence_per_100K_people)

temp<-reshape2::dcast(prev_stats_df_4cats_overall, 
                      metric + gender + category ~ range_indicator, value.var = "prevalence_per_100K_people")%>%
  select(c('metric', 'gender', 'category', 'mean', 'max', 'min'))

GBD_compare_df<-rbind(GBD_compare_df, temp)

names(prev_stats_df_4cats_gender)[names(prev_stats_df_4cats_gender) 
                                       == 'gender_name']<-'gender'
prev_stats_df_4cats_gender$metric = 'prevalence'
names(prev_stats_df_4cats_gender)[names(prev_stats_df_4cats_gender) 
                                       == 'group']<-'category'
prev_stats_df_4cats_gender$prevalence_per_100K_gender<-
  round(prev_stats_df_4cats_gender$prevalence_per_100K_gender)

temp<-reshape2::dcast(prev_stats_df_4cats_gender, 
                      metric + gender + category ~ range_indicator, value.var = "prevalence_per_100K_gender")%>%
  select(c('metric', 'gender', 'category', 'mean', 'max', 'min'))

GBD_compare_df<-rbind(GBD_compare_df, temp)


setwd(outdir_mse_dfs)
write.csv(GBD_compare_df, 'GBD_compare_df.csv', row.names = FALSE)

# write.csv(prev_stats_overall_df, 'prev_stats_overall_df.csv', row.names = FALSE)
# write.csv(prev_stats_df_4cats_overall, 'prev_stats_df_4cats_overall.csv', row.names = FALSE)
# write.csv(prev_stats_gender_df, 'prev_stats_gender_df.csv', row.names = FALSE)
# write.csv(prev_stats_df_4cats_gender, 'prev_stats_df_4cats_gender.csv', row.names = FALSE)

##########pull mort stats from 2017##########
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

setwd(outdir_mse_dfs)
write.csv(mort_stats_df, 'mort_stats_df.csv', row.names = FALSE)
write.csv(mort_stats_df_4cats, 'mort_stats_df_4cats.csv', row.names = FALSE)


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

setwd(outdir_mse_dfs)
write.csv(mort_calibration_df, 'mort_stats_df_4cats_per100K_gender_GBD_compare.csv', row.names = FALSE)

