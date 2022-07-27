#Estimating active TB mortalities per 100K to calibrate to
#Uses GBD 2019 - Number of Deaths and Scales to 100K
#KZN population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa
indir <- paste0(here(),'/param_files/calculated_param_gen/input_data')
outdir <- paste0(here(),'/param_files')

##########read in pop files#################
#data pull date 01/20/2022

setwd(indir) 
mort_rate_df<-read.csv('KZN_rate_IHME_deaths_1990_2017_Jan20.csv')
all_num_mort_df<-read.csv('KZN_num_IHME_deaths_1990_2017_Jan20.csv') #all HIV/TB deaths
TB_mort_num_df<-read.csv('KZN_TB_num_IHME_deaths_1990_2017_Jan20.csv')

TB_prev_num_df<-read.csv('KZN_TB_num_IHME_inc_prev_2017_Jan20.csv')%>%
  filter(measure_name == 'Prevalence')
TB_inc_num_df<-read.csv('KZN_TB_num_IHME_inc_prev_2017_Jan20.csv')%>%
  filter(measure_name == 'Incidence')

HIV_inc_num_df<-read.csv('GBD_KZN_num_hiv_inc.csv')

#get pop estimates by taking 1/rate * num of deaths (for all cause)
pop_df<-all_num_mort_df%>%
  filter(cause_id == 294)%>%
  left_join(mort_rate_df, by = c('age_id', 'sex_id', 'year'))%>%
  select('age_id', 'sex_id', 'year', 'val.x', 'val.y')%>% #val.x = num, val.y = rate
  mutate(pop_val = (val.x * (1/val.y)))

#group populations over all age groups
pop_df<-pop_df%>%
  group_by(sex_id, year)%>%
  summarise(expected_total_pop = sum(pop_val))

setwd(indir)
write.csv(pop_df, 'pop_df.csv', row.names = FALSE)

#TB mort estimates 1990-2017
total_mort_df<-TB_mort_num_df%>%
  filter(cause_id != 294)%>%
  mutate(co_infection = if_else(cause_id %in% c(946, 947, 934), 'TB_only', 'HIV/TB_coinfection'),
         calibration_group = paste0(co_infection, '_', sex_name))

total_mort_calibration_params<-total_mort_df%>%
  group_by(year, sex_id, calibration_group)%>%
  summarise(min_calib_total = sum(lower),
            max_calib_total = sum(upper),
            expected_calib_total = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex_id'))%>%
  mutate(min_rate = min_calib_total/expected_total_pop,
         max_rate =max_calib_total/expected_total_pop,
         expected_rate =expected_calib_total/expected_total_pop)%>%
  select(c('year', 'sex_id', 'calibration_group', 'min_rate', 'max_rate', 'expected_rate'))


#read in TB inc estimates
total_inc_df<-TB_inc_num_df%>%
  filter(cause_id != 294)%>%
  mutate(co_infection = if_else(cause_id %in% c(946, 947, 934), 'TB_only', 'HIV/TB_coinfection'),
         calibration_group = paste0(co_infection, '_', sex_name))

total_inc_df<-total_inc_df%>%
  group_by(year, sex_id, calibration_group)%>%
  summarise(min_calib_total = sum(lower),
            max_calib_total = sum(upper),
            expected_calib_total = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex_id'))%>%
  mutate(min_rate = min_calib_total/expected_total_pop,
         max_rate =max_calib_total/expected_total_pop,
         expected_rate =expected_calib_total/expected_total_pop)%>%
  select(c('year', 'sex_id', 'calibration_group', 'min_rate', 'max_rate', 'expected_rate'))


#read in TB prev estimates
total_prev_df<-TB_prev_num_df%>%
  filter(cause_id != 294)%>%
  mutate(co_infection = if_else(cause_id %in% c(946, 947, 934), 'TB_only', 'HIV/TB_coinfection'),
         calibration_group = paste0(co_infection, '_', sex_name))

total_prev_df<-total_prev_df%>%
  group_by(year, sex_id, calibration_group)%>%
  summarise(min_calib_total = sum(lower),
            max_calib_total = sum(upper),
            expected_calib_total = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex_id'))%>%
  mutate(min_rate = min_calib_total/expected_total_pop,
         max_rate =max_calib_total/expected_total_pop,
         expected_rate =expected_calib_total/expected_total_pop)%>%
  select(c('year', 'sex_id', 'calibration_group', 'min_rate', 'max_rate', 'expected_rate'))


#HIV incidence

names(HIV_inc_num_df)[names(HIV_inc_num_df) == 'sex'] <- "sex_id"

inc_rate_df<-HIV_inc_num_df%>%
  group_by(year, sex_id)%>%
  summarise(min_total = sum(lower),
            max_total = sum(upper),
            mean_total = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex_id'))%>%
  mutate(min_rate = min_total/(expected_total_pop*100000),
         max_rate =max_total/(expected_total_pop*100000),
         expected_rate =mean_total/(expected_total_pop*100000))

setwd(outdir)
write.csv(inc_rate_df, 'GBD_KZN_hiv_incidence_df.csv')


#if want to visualize mort estimates
#setwd(....)

# for (g in unique(total_mort_calibration_params$sex_name)){
#   
#   file_name <- paste0('TB_moratlity_g_compartment_', g, '.png')
#   
#   df_temp <- total_mort_calibration_params%>%filter(sex_name == g)
#   df1<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[1])
#   df2<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[2])
#   
#   graph_temp <- ggplot(df_temp,aes(x = year, y = expected_rate, color = calibration_group))+
#     geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), inherit.aes = FALSE,fill = "plum2")+
#     geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), inherit.aes = FALSE,fill = "lightgreen")+
#     geom_line(aes(colour = calibration_group), size = 1)+
#     labs(title = paste0('Deaths, rate per 100K, ', g))+
#     #scale_y_continuous(name="rate", breaks=seq(from = 0, to = 300, by = 20))+
#     scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
#     scale_color_manual(values=c('purple', 'green'))
#   
#   png(file_name)
#   print(graph_temp)
#   dev.off()
# }

setwd(outdir)
write.csv(total_mort_calibration_params, 'KZN_TB_mort_calibration_rates_df.csv', row.names = FALSE)
write.csv(total_inc_df, 'KZN_TB_inc_calibration_rates_df.csv', row.names = FALSE)

