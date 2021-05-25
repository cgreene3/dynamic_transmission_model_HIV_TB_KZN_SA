#Calibration Analysis April 29

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops'), require, character.only=T)


#############Set in directory and out directory###########
###########Make sure epi_model_HIV_TB.Rproj is open, otherwise will need to change wd manually########
date <- gsub('-', '_', Sys.Date())
indir <- paste0(here(),'/model_outputs/calibration/', '2021_05_12')
#paste0(here(),'/model_outputs/calibration/', date)
outdir <- paste0(here(),'/model_outputs/calibration/', '2021_05_12', '/analysis')
#paste0(here(),'/model_outputs/calibration/', date, '/analysis')

#will give warning if wd already exists
#ignore warning!
dir.create(file.path(outdir))

#file_list <- list.files(path=indir)
file_list<-c('out_df_sim_id_2916.csv',
             'out_df_sim_id_2841.csv',
             'out_df_sim_id_2860.csv',
             'out_df_sim_id_2891.csv',
             'out_df_sim_id_2890.csv')


#file_list<-file_list[file_list != "analysis"]
#file_list<-file_list[file_list != "sim_calibration_ref_df.csv"]
setwd(indir)
#sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')
mort_calib_TB_only_test <- c('L', 'M', 'H')
mort_calib_TB_HIV_CD4More_test <- c('L', 'M', 'H')
mort_calib_TB_HIV_CD4Less_test <- c('L', 'M', 'H')
mort_calib_TB_HIV_ART_test <- c('L', 'M', 'H')

#betas to test
beta_1_test<-seq(from = 10, to = 12, by = .4)
beta_2_test<-seq(from = 10, to = 12, by = .4)

#create a sim id for each combination
sim_calib_id<-c()
b1_calib_id<-c()
b2_calib_id<-c()
m1_calib_id<-c()
m2_calib_id<-c()
m3_calib_id<-c()
m4_calib_id<-c()

sim_id_temp<-1

for (b1 in beta_1_test){
  for (b2 in beta_2_test){
    for (m1 in mort_calib_TB_only_test){
      for (m2 in mort_calib_TB_HIV_CD4More_test){
        for (m3 in mort_calib_TB_HIV_CD4Less_test){
          for (m4 in mort_calib_TB_HIV_ART_test){
            sim_calib_id<-c(sim_calib_id, sim_id_temp)
            b1_calib_id<-c(b1_calib_id, b1)
            b2_calib_id<-c(b2_calib_id, b2)
            m1_calib_id<-c(m1_calib_id, paste0('1-', m1)) #to id proper
            m2_calib_id<-c(m2_calib_id, paste0('2-', m2)) #calibration params
            m3_calib_id<-c(m3_calib_id, paste0('3-', m3))
            m4_calib_id<-c(m4_calib_id, paste0('4-', m4))
            sim_id_temp<-sim_id_temp+1
          }
        }
      }
    }
  }
}

sim_calibration_ref_df<-data.frame(sim_calib_id, b1_calib_id, b2_calib_id,
                                   m1_calib_id, m2_calib_id, m3_calib_id, m4_calib_id)



calibration_df_all<-data.frame()

for (i in 1:length(file_list)){
  temp_data <- read.csv(file_list[i]) 
  calibration_df_all <- rbind(calibration_df_all, temp_data) #for each iteration, bind the new data to the building dataset
}


# #####Graphs that describe model states overtime####
calibration_df_all<-subset(calibration_df_all, select = -X)

#melt out all df for easy manipulation
state_prog_df <-reshape2::melt(calibration_df_all, id.vars = c("time", "year",
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

#####Calibration Calculations######
total_pop_in_gender_df<-state_prog_df%>%
  group_by(time, year, G_compartment, sim_id)%>%
  summarise(total_gender_pop = sum(value))%>%
  ungroup()%>%
  group_by(year, sim_id, G_compartment)%>%
  summarise(total_gender_pop = median(total_gender_pop))


total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  
mort_calibration_df<-calibration_df_all%>%
  select(c('year', 'sim_id', 'time', 'HIV_neg_male', 'HIV_neg_female', 'HIV_pos_male', 'HIV_pos_female'))%>%
  group_by(year, sim_id)%>%
  summarise(cum_mort_hiv_neg_male = max(HIV_neg_male),
         cum_mort_hiv_neg_female= max(HIV_neg_female),
         cum_mort_hiv_pos_male= max(HIV_pos_male),
         cum_mort_hiv_pos_female= max(HIV_pos_female))%>%
  as.data.frame()
# 
mort_calibration_df<-reshape2::melt(mort_calibration_df, id = c('year', 'sim_id'))
mort_calibration_df2<-mort_calibration_df%>%
  group_by(variable, sim_id)%>%
  mutate(mort_in_year = value - lag(value))%>%
  ungroup()%>%
  mutate(mort_in_year = if_else(year == 1990, value, mort_in_year))%>%
  filter(year < 2018)%>%
  mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
  left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
  mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
  as.data.frame()

mort_calibration_df<-mort_calibration_df2

setwd(paste0(here(), '/param_files'))
GBD_calibration_df<-read.csv('calibration_rates_df.csv')

#####overlap model outputs with calibration#####
setwd(indir)

mort_calibration_df_all<-mort_calibration_df[,c('year', 'variable', 'sim_id', 'mort_rate_per_100K')]

calibration_group_name<-c()
#
for (n in 1:nrow(mort_calibration_df_all)){
   if(mort_calibration_df_all$variable[n] ==
      'cum_mort_hiv_neg_male'){
     calibration_group_name<-c(calibration_group_name, 'TB_only_Male')
   } else if (mort_calibration_df_all$variable[n] ==
              'cum_mort_hiv_neg_female'){
     calibration_group_name<-c(calibration_group_name, 'TB_only_Female')
   } else if (mort_calibration_df_all$variable[n] ==
              'cum_mort_hiv_pos_female'){
     calibration_group_name<-c(calibration_group_name, 'HIV/TB_coinfection_Female')
   } else {
     calibration_group_name<-c(calibration_group_name, 'HIV/TB_coinfection_Male')
   }
 }
#
mort_calibration_df_all$calibration_group<-calibration_group_name
mort_calibration_df_all<-mort_calibration_df_all%>%
   rename(model_output = mort_rate_per_100K)%>%
   select(c('year', 'calibration_group', 'model_output', 'sim_id'))%>%
   left_join(GBD_calibration_df, by = c('year','calibration_group'))%>%
   filter(year < 2017)

setwd(outdir)

for (s in unique(mort_calibration_df_all$sim_id)){
  mort_calibration_df_temp<-mort_calibration_df_all%>%
    filter(sim_id == s)
  for (g in unique(mort_calibration_df_all$sex_name)){

  
  file_name <- paste0('TB_moratality_calibration', g, 'simid_', s, '.png')
  
  df_temp <- mort_calibration_df_temp%>%filter(sex_name == g)
  df_temp_model<-df_temp
  df_temp_model$calibration_group<-paste0('Model--', 
                                          df_temp_model$calibration_group)
  df_temp_GBD<-df_temp
  df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--', 
                                        df_temp_GBD$calibration_group)
  df1<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[1])
  df2<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[2])
  
  graph_temp <- ggplot()+
    geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE,fill = "plum2")+
    geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE,fill = "darkseagreen1")+
    geom_line(data = df_temp_model, aes(x = year, y = model_output, 
                                        color = calibration_group),
              linetype="dashed", size = 1.2)+
    geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate, 
                                      color = calibration_group))+
    scale_color_manual(values=c('green', 'purple', 'darkgreen', 'darkorchid4'))+
    labs(title = paste0('Deaths, rate per 100K, ', g))+
    scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
    scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
    theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  png(file_name, width = 480, height = 480)
  print(graph_temp)
  dev.off()
  }
}



for (g in unique(mort_calibration_df_all$sex_name)){
  
  
  file_name <- paste0('TB_moratality_calibration', g, 'overlap.png')
  
  df_temp <- mort_calibration_df_all%>%filter(sex_name == g)
  df_temp_model<-df_temp
  df_temp_model$calibration_group<-paste0('Model--', 
                                          df_temp_model$calibration_group)
  df_temp_GBD<-df_temp
  df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--', 
                                        df_temp_GBD$calibration_group)
  df1<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[1])
  df2<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[2])
  
  graph_temp <- ggplot()+
    geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE,fill = "plum2")+
    geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE,fill = "darkseagreen1")+
    geom_line(data = df_temp_model, aes(x = year, y = model_output, 
                                        color = calibration_group),
              linetype="dashed", size = 1.2)+
    geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate, 
                                      color = calibration_group))+
    scale_color_manual(values=c('green', 'purple', 'darkgreen', 'darkorchid4'))+
    labs(title = paste0('Deaths, rate per 100K, ', g))+
    scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
    scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
    theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  png(file_name, width = 480, height = 480)
  print(graph_temp)
  dev.off()
}



mse_df<-mort_calibration_df_all%>%
  mutate(diff = (model_output - expected_rate),
         mse = diff^2)%>%
  group_by(sim_id, calibration_group)%>%
  summarise(calibration_group_mse = sum(mse))%>%
  ungroup()

mse_df<-dcast(mse_df, sim_id ~ calibration_group)
mse_df$total_mse <-rowSums(mse_df[,2:5])

mse_df<-mse_df%>%
  mutate(rank_total = rank(total_mse),
         rank_1 = rank(`HIV/TB_coinfection_Male`),
         rank_2 = rank(`HIV/TB_coinfection_Female`),
         rank_3 = rank(TB_only_Male),
         rank_4 = rank(TB_only_Female))%>%
  left_join(sim_calibration_ref_df%>%
              rename(sim_id = sim_calib_id), 
            by = c('sim_id'))%>%
  select(-c(X))%>%
  as.data.frame()

mse_df2<-mse_df%>%select(-c(rank, sim_id))

library(randomForest)

mse.rf <- randomForest(mse_df2$total_mse ~ mse_df2$b1_calib_id+
                         mse_df2$b2_calib_id+
                         as.factor(mse_df2$m1_calib_id)+
                         as.factor(mse_df2$m2_calib_id)+
                         as.factor(mse_df2$m3_calib_id)+
                         as.factor(mse_df2$m4_calib_id))

varImpPlot(mse.rf)

mse_df_rank<-mse_df%>%
  filter(rank < nrow(mse_df)*.1)

mse_df_rank_betas<-subset(mse_df_rank, select = c(b1_calib_id, b2_calib_id, total_mse))
mse_df_rank_betas<-melt(mse_df_rank_betas)
mse_df_rank_betas<-mse_df_rank_betas%>%
  mutate(variable = if_else(variable == 'b1_calib_id', 'males',
                            'females'))%>%
  group_by(variable, value)%>%
  summarise(count = n())

ggplot(mse_df_rank_betas, aes(y=b1_calib_id, x=b2_calib_id, size = total_mse)) + 
  #geom_bar(position="dodge", stat="identity")+
  geom_point()+
  labs(x = 'effective contact rate (females)',
       y = 'effective contact rate (males)',
       size = 'total mean \nsquared error')
  #scale_fill_viridis(discrete = T)+
  labs(title = "total mse of top 10% of beta values \n for males and females")


plot(mse.rf)


mse_df_rank_mus<-subset(mse_df_rank, 
                        select = c(m1_calib_id, m2_calib_id, 
                                   m3_calib_id, m4_calib_id))
mse_df_rank_mus<-array(data = mse_df_rank_mus)
mse_df_rank_mus<-mse_df_rank_mus[1:length(mse_df_rank_mus)]



G_SET<-1:2

for (g in unique(mort_calibration_df_all$sex_name)){

  file_name <- paste0('TB_moratality_calibration', g, '.png')

  df_temp <- mort_calibration_df_all%>%filter(sex_name == g, sim_id == 1)
  df_temp_model<-df_temp
  df_temp_model$calibration_group<-paste0('Model--',
                                          df_temp_model$calibration_group)
  df_temp_GBD<-df_temp
  df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--',
                                        df_temp_GBD$calibration_group)
  df1<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[1])
  df2<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[2])

  graph_temp <- ggplot()+
    geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate),
                inherit.aes = FALSE,fill = "plum2")+
    geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate),
                inherit.aes = FALSE,fill = "darkseagreen1")+
    geom_line(data = df_temp_model, aes(x = year, y = model_output,
                  color = calibration_group),
                  linetype="dashed", size = 1.2)+
    geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate,
                  color = calibration_group))+
    scale_color_manual(values=c('green', 'purple', 'darkgreen', 'darkorchid4'))+
    labs(title = paste0('Deaths, rate per 100K, ', g))+
    scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
    scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
    theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

  png(file_name, width = 480, height = 480)
  print(graph_temp)
  dev.off()
}




#melt out all df for easy manipulation
state_prog_df <-reshape2::melt(melted_df, id.vars = c("time", "year",
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

#####Calibration Calculations######
total_pop_in_gender_df<-state_prog_df%>%
  group_by(time, year, G_compartment, sim_id)%>%
  summarise(total_gender_pop = sum(value))%>%
  ungroup()%>%
  group_by(year, sim_id, G_compartment)%>%
  summarise(total_gender_pop = median(total_gender_pop))




# 
# 
# 
# 
# 
# 
# 
#   
#   
#   #data.frame(diff(as.matrix(mort_df$cum_mort_hiv_pos_female)))
# # 
# # HIV_overtime<-out_df_melt%>%
# #   mutate(HIV_temp = if_else(HIV_compartment > 1, value, 0))%>%
# #   group_by(time, G_compartment)%>%
# #   mutate(total_in_g_compartment = sum(value))%>%
# #   group_by(time, G_compartment, total_in_g_compartment)%>%
# #   summarise(total_hiv_pos = sum(HIV_temp))%>%
# #   mutate(hiv_prev = total_hiv_pos/total_in_g_compartment)
# # 
# # #png("HIV_prev_overtime_reactivation.png")
# # #hivplot1 <- ggplot(data = HIV_overtime, 
# #                    mapping = aes(x = time, y = hiv_prev, color = G_compartment))+
# #   geom_line()+
# #   lims(y = c(0,1), x = c(0,27))+
# #   labs(title = 'REACTIVATION (based on remote rates) - HIV Prev Over Time')
# # print(hivplot1)
# # dev.off()
# # 
# # 
# # ART_coverage_overtime<-out_df_melt%>%
# #   group_by(time, HIV_compartment, G_compartment)%>%
# #   summarise(total = sum(value))%>%
# #   filter(HIV_compartment > 1)%>%
# #   group_by(time, G_compartment)%>%
# #   mutate(total_PLHIV = sum(total))%>%
# #   filter(HIV_compartment == 4)%>%
# #   mutate(ART_coverage = total/total_PLHIV)
# # 
# # png("ART_coverage_overtime_reactivation.png")
# # artcov_plot <- ggplot(data = ART_coverage_overtime, 
# #                       mapping = aes(x = time, y = ART_coverage, color = G_compartment))+
# #   geom_line()
# # print(artcov_plot)
# # dev.off()
# # 
# # 
# # 
# # 
# # # 
# # #   
# # # 
# # # 
# # # 
# # # 
# # # #hiv compartments overtime test to compare with caras coverages
# # # out_df_by_HIV<-out_df_melt%>%
# # #   group_by(HIV_compartment, G_compartment, time)%>%
# # #   summarise(total_hiv_pop = sum(value))
# # # 
# # # #testing hiv prevalence
# # # ggplot(data = out_df_by_HIV, 
# # #        mapping = aes(x = time, y = total_pop, fill = HIV_compartment))+
# # #   geom_area()
# # # 
# # # 
# # # #filter only active TB compartments, since that is what we are calibrating to
# # # out_df_TB_active<-out_df_melt%>%
# # #   filter(TB_compartment == 6)%>%
# # #   mutate(calibration_group_name = if_else((HIV_compartment == 1 & G_compartment == 1),
# # #                                      'TB_only_Male',
# # #                                      if_else((HIV_compartment != 1 & G_compartment == 1),
# # #                                              'HIV/TB_coinfection_Male',
# # #                                              if_else((HIV_compartment == 1 & G_compartment == 2),
# # #                                              'TB_only_Female',
# # #                                              'HIV/TB_coinfection_Female'))),
# # #          calibration_group_id = if_else((HIV_compartment == 1 & G_compartment == 1),
# # #                                      1,
# # #                                      if_else((HIV_compartment != 1 & G_compartment == 1),
# # #                                              2,
# # #                                              if_else((HIV_compartment == 1 & G_compartment == 2),
# # #                                                      3,
# # #                                                      4))))
# # # 
# # # #so I do not need to call on mort_param_func too many times
# # # out_df_TB_active<-out_df_TB_active%>%
# # #   arrange(year)
# # # 
# # # #testing TB prev overtime
# # # test<-out_df_TB_active%>%
# # #   group_by(G_compartment, HIV_compartment, time)%>%
# # #   summarise(value = sum(value))%>%
# # #   mutate(ID = paste0('HIV_', HIV_compartment, '_G_', G_compartment))
# # # 
# # # ggplot(data = test, mapping = aes(x = time, y = value, fill = ID))+
# # #   geom_area()
# # # 
# # # 
# # # TB_grouping_test_df<-out_df_melt%>%
# # #   group_by(TB_compartment, time)%>%
# # #   summarise(total_pop = sum(value))
# # # 
# # # 
# # # ggplot(data = TB_grouping_test_df%>%filter(TB_compartment == 6), mapping = aes(x = time, y = total_pop,
# # #                                                  fill = TB_compartment))+
# # #   geom_area()
# # # 
# # # counter <-1
# # # mort_est <- rep(0, times = nrow(out_df_TB_active))
# # # 
# # # for (yr in start_yr:end_yr){
# # #   temp<-out_df_TB_active%>%
# # #     filter(year == yr)
# # #   mu_t_h_g<-mort_param_func(yr)
# # #   for (row in 1:nrow(temp)){
# # #     hiv <- as.integer(temp[row, 'HIV_compartment'])
# # #     gender <-as.integer(temp[row, 'G_compartment'])
# # #     pop<-as.double(temp[row, 'value'])
# # #     mort_rate<-(mu_t_h_g[6,hiv,gender]*(1/12))
# # #     
# # #     mort_est[counter]<-(mort_rate*pop)
# # #     counter <- counter + 1
# # #   }
# # # }
# # # 
# # # out_df_TB_active$mort_est <- mort_est
# # # 
# # # #group and combine data for calibration
# # # setwd(indir)
# # # calibration_rates_df<-read.csv('calibration_rates_df.csv')
# # # 
# # # calibration_df<-out_df_TB_active%>%
# # #   group_by(calibration_group_id, calibration_group_name, 
# # #            year)%>%
# # #   summarise(model_rate = sum(mort_est))%>%
# # #   left_join(calibration_rates_df, by = c('calibration_group_name', 'year'))%>%
# # #   filter(year < 2017)
# # #   #mutate(time = time+1990) #so that time graphs go from 1990 - 2017
# # # 
# # # calibration_df$within<-if_else((calibration_df$model_rate<=calibration_df$max_rate)&(calibration_df$model_rate>=calibration_df$min_rate),
# # #                                1, 0)
# # # calibration_df$diff <-calibration_df$model_rate-calibration_df$expected_rate
# # # calibration_df$mse <-(calibration_df$diff)^2
# # # 
# # # 
# # # 





