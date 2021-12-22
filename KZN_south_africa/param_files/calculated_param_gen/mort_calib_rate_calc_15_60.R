#Estimating active TB mortalities per 100K to calibrate to
#Uses GBD 2019 - Number of Deaths and Scales to 100K
#SA population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa
indir <- paste0(here(),'/param_files/calculated_param_gen/input_data')
indir_pop <- paste0(indir, '/GBD_pop_estimates')
outdir <- paste0(here(),'/param_files')

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

#write.csv(pop_df, 'pop_df.csv')

#read pop
#pop_df<-read.csv('pop_df.csv')

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
  left_join(pop_df, by = c('year', 'sex_name'))%>%
  mutate(min_percent = min_calib_total/expected_total_pop,
         max_percent =max_calib_total/expected_total_pop,
         expected_percent =expected_calib_total/expected_total_pop,
         min_rate = min_percent*100000,
         max_rate = max_percent*100000,
         expected_rate = expected_percent*100000)%>%
  select(c('year', 'sex_name', 'calibration_group', 'min_rate', 'max_rate', 'expected_rate'))

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
write.csv(total_mort_calibration_params, 'TB_mort_calibration_rates_df.csv', row.names = FALSE)

