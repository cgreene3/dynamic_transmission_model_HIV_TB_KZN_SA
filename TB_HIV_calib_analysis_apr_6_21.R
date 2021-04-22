#Calibration Analysis April 6

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops'), require, character.only=T)

#######Create calibration ref df#########
#mort calibration test
#(1 = low, 2 = medium, 3 = high)
#TB_only, TB_HIV_CD4More, TB_HIV_CD4Less, TB_HIV_ART
mort_calib_TB_only_test <- 1:3
mort_calib_TB_HIV_CD4More_test <- 1:3
mort_calib_TB_HIV_CD4Less_test <- 1:3
mort_calib_TB_HIV_ART_test <- 1:3
  
#betas to test
beta_1_test<-seq(from = 4, to = 6, by = .4)
beta_2_test<-seq(from = 4, to = 6, by = .4)

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
            m1_calib_id<-c(m1_calib_id, m1)
            m2_calib_id<-c(m2_calib_id, m2)
            m3_calib_id<-c(m3_calib_id, m3)
            m4_calib_id<-c(m4_calib_id, m4)
            sim_id_temp<-sim_id_temp+1
          }
        }
      }
    }
  }
}

sim_calibration_ref_df<-data.frame(sim_calib_id, b1_calib_id, b2_calib_id,
                                   m1_calib_id, m2_calib_id, m3_calib_id, m4_calib_id)

#remove everything but df
rm(list=setdiff(ls(), "sim_calibration_ref_df"))

#############Set in directory and out directory###########
###########Make sure epi_model_HIV_TB.Rproj is open, otherwise will need to change wd manually########
date <- gsub('-', '_', Sys.Date())
indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/model_outputs/calibration/', date)

#will give warning if wd already exists
#ignore warning!
dir.create(file.path(outdir))
setwd(outdir)


# #####Graphs that describe model states overtime####
# state_prog_df<-out_df%>%
#   select(-c(calibration_mort_states))
# 
# #melt out all df for easy manipulation
# state_prog_df <-melt(state_prog_df, id.vars = c("time", "year",
#                                                'sim_id', "beta_1", "beta_2"))
# state_prog_df <- cbind(state_prog_df, 
#                      data.frame(do.call('rbind', 
#                                         strsplit(as.character(state_prog_df$variable),
#                                                  '_',fixed=TRUE))))
# names(state_prog_df)[names(state_prog_df) == "X2"] <- "TB_compartment"
# names(state_prog_df)[names(state_prog_df) == "X3"] <- "DR_compartment"
# names(state_prog_df)[names(state_prog_df) == "X4"] <- "HIV_compartment"
# names(state_prog_df)[names(state_prog_df) == "X5"] <- "G_compartment"
# state_prog_df<-state_prog_df%>%select(-c('X1'))
# state_prog_df$month <- round(((state_prog_df$time%%1)*(12)+1),0)
# 
# outdir_state_prog <- paste0(here(),'/model_outputs/state_prog/date')
# setwd(outdir_state_prog)
# 
# ######TB states########
# TB_all_overtime<-state_prog_df%>%
#   group_by(TB_compartment, time)%>%
#   summarise(total_in_compartment = sum(value))%>%
#   filter(TB_compartment != 'pop')%>%
#   mutate(time = time+start_yr)
# 
# tb_all_graph<-ggplot(data = TB_all_overtime, 
#                      mapping = aes(x = time, 
#                                    y = total_in_compartment, 
#                                    color = TB_compartment))+
#   geom_line()+
#   labs(title = 'All TB compartments overtime')
# 
# png("tb_all_overtime_after_warmup.png")
# print(tb_all_graph)
# dev.off()
# 
# TB_recovered_overtime<-state_prog_df%>%
#   filter(TB_compartment)%>%
#   group_by(HIV_compartment, time)%>%
#   summarise(total_in_compartment = sum(value))%>%
#   mutate(time = time+start_yr)
# 
# ggplot(data = TB_recovered_overtime, 
#        mapping = aes(x = time, 
#                      y = total_in_compartment, 
#                      color = HIV_compartment))+
#   geom_line()+
#   labs(title = 'TB compartments overtime')
# 
# tb_active_graph_summarised<-ggplot()+
#   +
#   labs(title = 'Active TB overtime')
# 
# png("tb_active_overtime_after_warmup.png")
# print(tb_active_graph_summarised)
# dev.off()
# 
# TB_prev_overtime_gender_HIV<-state_prog_df%>%
#   group_by(time, G_compartment)%>%
#   mutate(total_g_compartment = sum(value))%>%
#   ungroup()%>%
#   filter(TB_compartment == 6)%>%
#   group_by(time, HIV_compartment, G_compartment)%>%
#   mutate(total_pop_in_compartment = sum(value),
#          time = time+1990,
#          percent_pop_in_compartment = 
#            total_pop_in_compartment/total_g_compartment,
#          rate = percent_pop_in_compartment*100000)%>%
#   ungroup()%>%
#   group_by(year, HIV_compartment, G_compartment)%>%
#   mutate()
# 
# for (g in G_SET){
#   
#   file_name <- paste0('ActiveTB_g_compartment_', g, '.png')
#   
#   gender_name = if_else(g == 1, 'Males', 'Females')
#   
#   graph_temp <-ggplot() +
#     geom_line(TB_prev_overtime_gender_HIV%>%filter(G_compartment == g),
#               mapping = aes(x = time, y = rate, 
#                             color = HIV_compartment), size = 1)+
#     labs(title = paste0('TB Prevalence, rate per 100K,\nby HIV compartment for ', gender_name))+
#     scale_y_continuous(name=paste0("Rate per 100k, ", gender_name), limits=c(0, 600))+
#     scale_x_continuous(name = 'Time', limits = c(1990, 2017), breaks=(seq(1990, 2017, 5)))+
#     scale_color_manual(values=c("red", "#56B4E9", "purple", "green4"))+
#     #theme(legend.position='none', text = element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     #      panel.background = element_blank(), axis.line = element_line(colour = "black"))
#   
#   png(file_name)
#   print(graph_temp)
#   dev.off()
# }
# 
# 
# #####Calibration Calculations######
# total_pop_in_gender_df<-state_prog_df%>%
#   group_by(time, year, G_compartment)%>%
#   summarise(total_gender_pop = sum(value))%>%
#   ungroup()%>%
#   group_by(year, G_compartment)%>%
#   summarise(total_gender_pop = median(total_gender_pop))
# 
# 
# mort_calibration_df<-out_df%>%
#   select(c('year', 'time', 'HIV_neg_male', 'HIV_neg_female', 'HIV_pos_male', 'HIV_pos_female'))%>%
#   group_by(year)%>%
#   summarise(cum_mort_hiv_neg_male = max(HIV_neg_male),
#          cum_mort_hiv_neg_female= max(HIV_neg_female),
#          cum_mort_hiv_pos_male= max(HIV_pos_male),
#          cum_mort_hiv_pos_female= max(HIV_pos_female))
# 
# mort_calibration_df<-melt(mort_calibration_df, id = c('year'))
# mort_calibration_df<-mort_calibration_df%>%
#   group_by(variable)%>%
#   mutate(mort_in_year = value - lag(value))%>%
#   ungroup()%>%
#   mutate(mort_in_year = if_else(year == 1990, value, mort_in_year))%>%
#   filter(year < 2018)%>%
#   mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
#   left_join(total_pop_in_gender_df, by = c('G_compartment', 'year'))%>%
#   mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))
# 
# for (g in G_SET){
#   
#   file_name <- paste0('TB_moratlity_g_compartment_', g, '.png')
#   
#   gender = if_else(g == 1, 'Males', 'Females')
#   
#   graph_temp <- ggplot(data = mort_calibration_df%>%filter(G_compartment == g), 
#                        mapping = aes(x = year, y = mort_rate_per_100K))+
#     geom_line(aes(colour = variable), size = 1)+
#     labs(title = paste0('Deaths, rate per 100K, ', gender))+
#     #scale_y_continuous(name="rate", breaks=seq(from = 0, to = 300, by = 20))+
#     scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
#     scale_color_manual(values=c("green", 'purple'))+
#     ylim(0, 600)
#   
#   png(file_name)
#   print(graph_temp)
#   dev.off()
# }
# 
# #####overlap model outputs with calibration#####
# setwd(indir)
# GBD_calibration_df<-read.csv('calibration_rates_df.csv')
# mort_calibration_df_all<-mort_calibration_df%>%
#   select(c('year', 'variable', 'mort_rate_per_100K'))
# 
# calibration_group_name<-c()
# 
# for (n in 1:nrow(mort_calibration_df_all)){
#   if(mort_calibration_df_all$variable[n] == 
#      'cum_mort_hiv_neg_male'){
#     calibration_group_name<-c(calibration_group_name, 'TB_only_Male')
#   } else if (mort_calibration_df_all$variable[n] == 
#              'cum_mort_hiv_neg_female'){
#     calibration_group_name<-c(calibration_group_name, 'TB_only_Female')
#   } else if (mort_calibration_df_all$variable[n] == 
#              'cum_mort_hiv_pos_female'){
#     calibration_group_name<-c(calibration_group_name, 'HIV/TB_coinfection_Female')
#   } else {
#     calibration_group_name<-c(calibration_group_name, 'HIV/TB_coinfection_Male')
#   }
# }
# 
# mort_calibration_df_all$calibration_group<-calibration_group_name
# mort_calibration_df_all<-mort_calibration_df_all%>%
#   rename(model_output = mort_rate_per_100K)%>%
#   select(c('year', 'calibration_group', 'model_output'))%>%
#   left_join(GBD_calibration_df, by = c('year','calibration_group'))%>%
#   filter(year < 2017)
# 
# setwd(outdir_state_prog)
# for (g in unique(mort_calibration_df_all$sex_name)){
#   
#   file_name <- paste0('TB_moratality_calibration', g, '.png')
#   
#   df_temp <- mort_calibration_df_all%>%filter(sex_name == g)
#   df_temp_model<-df_temp
#   df_temp_model$calibration_group<-paste0('Model--', 
#                                           df_temp_model$calibration_group)
#   df_temp_GBD<-df_temp
#   df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--', 
#                                         df_temp_GBD$calibration_group)
#   df1<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[1])
#   df2<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[2])
#   
#   graph_temp <- ggplot()+
#     geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), 
#                 inherit.aes = FALSE,fill = "plum2")+
#     geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), 
#                 inherit.aes = FALSE,fill = "darkseagreen1")+
#     geom_line(data = df_temp_model, aes(x = year, y = model_output, 
#                   color = calibration_group),
#                   linetype="dashed", size = 1.2)+
#     geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate, 
#                   color = calibration_group))+
#     scale_color_manual(values=c('green', 'purple', 'darkgreen', 'darkorchid4'))+
#     labs(title = paste0('Deaths, rate per 100K, ', g))+
#     scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
#     scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
#     theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"))
#   
#   png(file_name, width = 480, height = 480)
#   print(graph_temp)
#   dev.off()
# }
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
