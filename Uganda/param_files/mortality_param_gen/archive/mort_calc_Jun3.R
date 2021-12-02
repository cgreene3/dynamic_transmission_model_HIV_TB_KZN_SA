#Estimating background mortality rate from GBD
#Uses GBD 2019 - Number of Deaths
#SA population estimates

#clean workspace
rm(list = setdiff(ls(), 'df_input_data_all'))
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#to name plots according to date generated
data_gen_date<-Sys.Date()
data_gen_date<-str_replace_all(data_gen_date, '-', '_')

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/mortality_param_gen')
indir_pop <- paste0(here(),'/param_files/GBD_pop_estimates')
outdir <- paste0(here(),'/param_files/calib_ref_data/', data_gen_date)
#will give warning if wd already exists
#ignore warning!
dir.create(file.path(outdir))


setwd(indir)

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
pop_df<-read.csv('pop_df.csv')

#prevalence and mortality data frames consider the following disease states
#HIV/AIDS - causeid 298
#DS TB - cause id 934
#MDR TB - cause ids 946, 947

#########read in prevalence data#############
setwd(indir)
prev_df<-read.csv('IHME_prev_1990_2017_Feb28.csv')
prev_df<-prev_df%>%
  group_by(sex_id, year)%>%
  summarise(expected_disease_pop = sum(val))
  
#use prevalence data to calculate non disease mort
non_disease_mort_df<-pop_df%>%
  left_join(prev_df, by = c('sex_id', 'year'))%>%
  mutate(expected_non_disease_pop = expected_total_pop-expected_disease_pop)

##########read in mortality data##########
mort_df<-read.csv('IHME_deaths_1990_2017_Feb28.csv')%>%
  mutate(disease_cat = if_else(cause_name == 'All causes', 'base', 'disease'))%>%
  group_by(sex_id, year, disease_cat)%>%
  summarise(expected_disease_mort = sum(val))%>%
  ungroup()%>%
  mutate(base_mort_calcs = if_else(disease_cat == 'base', 
                                   expected_disease_mort, 
                                   -expected_disease_mort))%>%
  group_by(sex_id, year)%>%
  summarise(expected_non_disease_mort_total = sum(base_mort_calcs))

non_disease_mort_df<-non_disease_mort_df%>%
  left_join(mort_df, by = c('sex_id', 'year'))%>%
  mutate(expected_non_disease_mort_rate = expected_non_disease_mort_total/expected_non_disease_pop)

#manipulate for mortality calcs
mort_graph_IHME<-read.csv('IHME_deaths_1990_2017_Feb28.csv')%>%
  mutate(Disease_Category = if_else(cause_id %in% c(947, 949, 934), 'Tuberculosis',
                               if_else(cause_id == 298, 'HIV/AIDS',
                                       'All Causes')))%>%
  group_by(sex_id, year, Disease_Category)%>%
  summarise(expected_mort = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex_id'))%>%
  mutate(per_100K = (expected_mort/expected_total_pop)*100000)%>%
  select(c('sex_id', 'year', 'Disease_Category', 'sex_name', 'per_100K'))

non_disease_mort_df2<-non_disease_mort_df%>%
  mutate(Disease_Category = 'Baseline',
         per_100K = expected_non_disease_mort_rate*100000)%>%
  select(c('sex_id', 'year', 'Disease_Category', 'sex_name', 'per_100K'))
# 
# mort_graph_IHME<-rbind(mort_graph_IHME, non_disease_mort_df2)
# 
# 
# p1<-ggplot2::ggplot(data = mort_graph_IHME%>%filter(sex_id == 1), 
#                 aes(x = year, y = per_100K))+
#   geom_line(aes(group = Disease_Category, color = Disease_Category))+
#   #ggtitle('male')+
#   ylab('mortality rate, per 100K males')+
#   scale_x_continuous(breaks=seq(1990,2017,3))+
#   theme(axis.title = element_text(size = 14),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12))+
#   theme_calc()+ scale_colour_calc()+ theme(
#     # Hide panel borders and remove grid lines
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Change axis line
#     axis.line = element_line(colour = "black"),
#     legend.position = "bottom"
#   ) +
#   guides(color=guide_legend(nrow=1, byrow = TRUE))
# 
# 
# p2<-ggplot2::ggplot(data = mort_graph_IHME%>%filter(sex_id == 2), 
#                     aes(x = year, y = per_100K))+
#   geom_line(aes(group = Disease_Category, color = Disease_Category))+
#   #ggtitle('male')+
#   ylab('mortality rate, per 100K females')+
#   scale_x_continuous(breaks=seq(1990,2017,3))+
#   theme(axis.title = element_text(size = 14),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12))+
#   theme_calc()+ scale_colour_calc()+ theme(
#     # Hide panel borders and remove grid lines
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     # Change axis line
#     axis.line = element_line(colour = "black"),
#     legend.position = "bottom"
#   ) +
#   guides(color=guide_legend(nrow=1, byrow = TRUE))
# 
# 
# non_disease_mort_graph<-ggplot(non_disease_mort_df, aes(x = year, 
#                                 y = expected_non_disease_mort_rate,
#                                 group = sex_name)) + 
#   geom_point() + 
#   geom_line(aes(color = sex_name))+
#   labs(title = "Non-disease specific mortality",
#        y = 'Non-disease specific mortality rate',
#        colour = "Gender")+
#   ylim(0,.02)
# 
# png(paste0('non_disease_mort_graph_',data_gen_date,'.png'), width=450,height=350,res=100)
# print(non_disease_mort_graph)
# dev.off()
# 
# #write data to read into epi model

#calculate base death rates by gender
df_input_data <- non_disease_mort_df

TB_SET <- 1:8
HIV_SET <- 1:4
#calculate disease specific mort
df_input_data<-do.call("rbind", replicate(length(TB_SET)*length(HIV_SET), 
                                          df_input_data, 
                                          simplify = FALSE))

df_input_data<-df_input_data%>%
  mutate(TB_compartment = 0,
         HIV_compartment = 0)%>%
  rename(G_compartment = sex_id)

df_input_data$year<-as.integer(df_input_data$year)
df_input_data<-df_input_data%>%
    arrange(year)%>%
    arrange(G_compartment)

counter <-1

for (yr in unique(df_input_data$year)){
    for (g in 1:2){
        for (t in TB_SET){
            for (h in HIV_SET){
                df_input_data$TB_compartment[counter]<-t
                df_input_data$HIV_compartment[counter]<-h
                counter<-counter+1
            }
        }
    }
}



#adjusted mortality
HIV_2_increase <- c(10)
HIV_3_increase <- c(35, 50)
HIV_4_increase <- c(1.2, 1.5)

#L-H
#TB (HIV neg) – 14, 18
#TB, HIV CD4>200 – 20, 25
#TB, HIV CD4<200 – 30, 40
#TB, HIV on ART - 17, 22
#L-M-H
#TB (HIV neg) - 10, 15, 20
#TB, HIV CD4>200 – 15, 22, 30
#TB, HIV CD4<200 - 30, 40, 50
#TB, HIV on ART - 12, 18, 24


TB_Active_HIV_1<-c(10,15)
TB_Active_HIV_2<-c(22,30) 
TB_Active_HIV_3<-c(50,70)
TB_Active_HIV_4<-c(15,20,25)

df_input_data_all<-df_input_data%>%
  filter(HIV_compartment == 1,
         TB_compartment != 6)

#no HIV no TB adj
df_input_data_all$mort_adj = rep(1, times = nrow(df_input_data_all))
df_input_data_all$level = rep(0, times = nrow(df_input_data_all))
df_input_data_all$level_ref = rep('base', times = nrow(df_input_data_all))

#hiv ONLY adjustments
for (h2 in 1:length(HIV_2_increase)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 2,
           TB_compartment != 6)%>%
    mutate(mort_adj = HIV_2_increase[h2],
           level = h2,
           level_ref = paste0('hiv_2_only_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}

for (h3 in 1:length(HIV_3_increase)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 3,
           TB_compartment != 6)%>%
    mutate(mort_adj = HIV_3_increase[h3],
           level = h3,
           level_ref = paste0('hiv_3_only_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}

for (h4 in 1:length(HIV_4_increase)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 4,
           TB_compartment != 6)%>%
    mutate(mort_adj = HIV_4_increase[h4],
           level = h4,
           level_ref = paste0('hiv_4_only_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}


#TB adjustments
for (tb1 in 1:length(TB_Active_HIV_1)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 1,
           TB_compartment == 6)%>%
    mutate(mort_adj = TB_Active_HIV_1[tb1],
           level = tb1,
           level_ref = paste0('TB_only_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}

for (tb2 in 1:length(TB_Active_HIV_2)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 2,
           TB_compartment == 6)%>%
    mutate(mort_adj = TB_Active_HIV_2[tb2],
           level = tb2,
           level_ref = paste0('TB_HIV2_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}

for (tb3 in 1:length(TB_Active_HIV_3)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 3,
           TB_compartment == 6)%>%
    mutate(mort_adj = TB_Active_HIV_3[tb3],
           level = tb3,
           level_ref = paste0('TB_HIV3_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}

for (tb4 in 1:length(TB_Active_HIV_4)){
  df_input_data_temp<-df_input_data%>%
    filter(HIV_compartment == 4,
           TB_compartment == 6)%>%
    mutate(mort_adj = TB_Active_HIV_4[tb4],
           level = tb4,
           level_ref = paste0('TB_HIV4_level_', level))
  
  df_input_data_all<-rbind(df_input_data_temp, df_input_data_all)
}

df_input_data_all$mort_rate<-df_input_data_all$expected_non_disease_mort_rate*
  df_input_data_all$mort_adj


df_input_data_all<-df_input_data_all%>%
  dplyr::select(c('year', 'TB_compartment',
                  'HIV_compartment', 'G_compartment', 
                  'mort_adj', 'level', 'level_ref', 
                  'mort_rate'))

setwd(outdir)
write.csv(df_input_data_all, 'mort_df.csv', row.names = FALSE)

