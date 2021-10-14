#Calibration Analysis July 21th to 2028

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
###########Make sure epi_model_HIV_TB.Rproj is open, otherwise will need to change wd manually########
start_eval_date <- '2021_07_21' #may need to change manually 
desc<-'/past_2018_changing_births'
#if analysis date is later that start date (running calib prog)
description <- '1990_2028_ANALYSIS_changing_births'

indir_ref_data<-paste0(here(),'/param_files/calibration_code_results/', start_eval_date)

#set calib outputs to local file on computer (otherwise files will overwhelm github)
indir_state_prog<-paste0("~/Documents/academic_posttt_2020/HIV_TB/calib_outdir/", start_eval_date, desc)

outdir_analysis <- paste0(here(),'/param_files/calibration_code_results/', start_eval_date, '/analysis/', description)
dir.create(paste0(here(),'/param_files/calibration_code_results/', start_eval_date, '/analysis/'))
dir.create(file.path(outdir_analysis))

#######Create calibration ref df#########
setwd(indir_ref_data)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')

# #####Graphs that describe model states overtime####

setwd(paste0(here(), '/param_files'))
HIV_calibration_df<-read.csv('hiv_transmission_df_v4.csv')
HIV_calibration_df<-HIV_calibration_df%>%
  mutate(G_compartment = if_else(gender == 'Females', 2, 1))%>%
  dplyr::select('year', 'G_compartment', 'hiv_prevalence')

HIV_ART_cov_df<-read.csv('hiv_transmission_df_v4.csv')
HIV_ART_cov_df<-HIV_ART_cov_df%>%
  mutate(G_compartment = if_else(gender == 'Females', 2, 1))%>%
  dplyr::select('year', 'G_compartment', 'art_coverage')%>%
  rename(art_coverage_HPV_est = art_coverage)


#GBD_calibration_df<-read.csv('calibration_rates_df.csv')
#match calibration group with model output calibration groups
#GBD_calibration_df$calibration_group<-if_else(GBD_calibration_df$calibration_group == "TB_only_Female",
#                                              "HIV_neg_female",
#                                              if_else(GBD_calibration_df$calibration_group == "TB_only_Male",
#                                                      "HIV_neg_male",
#                                                      if_else(GBD_calibration_df$calibration_group == "HIV/TB_coinfection_Male",
#                                                              "HIV_pos_male", "HIV_pos_female")))


#GBD_calibration_df$G_compartment <-if_else(GBD_calibration_df$sex_name == 'male',
#                                           1, 2)

#####Evaluating Top MSE#####

####TOP GRAPHS FUNCTION#####
top_graphs_func<-function(mse_df_top){
  
  calibration_df_all_top<-data.frame()
  state_prog_df<-data.frame()
  total_pop_in_gender_df<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in mse_df_top$sim_id){
  
    file_name<-paste0('out_df_sim_id_', sim_id_itr, '_to_2028.csv')
    print(sim_id_itr)
    temp_data <- fread(file_name) 
    temp_data$sim_id <- rep(sim_id_itr, times = nrow(temp_data))
    calibration_df_all_top <- rbind(calibration_df_all_top, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  state_prog_df <-reshape2::melt(calibration_df_all_top, 
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
  
  total_pop_in_gender_df<-state_prog_df%>%
    filter(G_compartment %in% c(1,2))%>%
    group_by(time, year, G_compartment, sim_id)%>%
    summarise(total_gender_pop = sum(value))%>%
    ungroup()%>%
    group_by(year, sim_id, G_compartment)%>%
    summarise(total_gender_pop = median(total_gender_pop))
  
  total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  ###hiv prev calcs####
  for (sim_id_itr in mse_df_top$sim_id){
    setwd(outdir_best_analysis_HIV_PREV)
    file_name_hiv_prev_graph<-paste0('hiv_prev_sim_id_', sim_id_itr, '.png')
    HIV_prev_df<-state_prog_df%>%
      filter(sim_id == sim_id_itr)%>%
      filter(HIV_compartment %in% c(2,3,4))%>%
      group_by(year, time, G_compartment)%>%
      summarise(total_prev = sum(value))%>%
      group_by(year, G_compartment)%>%
      summarise(total_prev = median(total_prev))%>%
      left_join(total_pop_in_gender_df%>%filter(sim_id == sim_id_itr), by = c('year', 'G_compartment'))%>%
      mutate(HIV_prev_TB_model_estimates_per_100K = total_prev*(100000/total_gender_pop))%>%
      mutate(G_compartment = as.integer(G_compartment))%>%
      left_join(HIV_calibration_df, by = c('year', 'G_compartment'))%>%
      mutate(HIV_prev_HPV_model_estimates_per_100K = hiv_prevalence*100000)%>%
      select(c('year', 'G_compartment', 'HIV_prev_TB_model_estimates_per_100K',
               'HIV_prev_HPV_model_estimates_per_100K'))
    
    
    HIV_prev_df_melt <-reshape2::melt(HIV_prev_df, 
                                   id.vars = c("year",
                                               'G_compartment'))
    
    HIV_prev_df_melt$calibration_group<-(paste0('Gender_', HIV_prev_df_melt$G_compartment, '_',
                                     HIV_prev_df_melt$variable))
    
    
    hiv_prev_graph<-ggplot()+
      geom_line(data = HIV_prev_df_melt,
                aes(x = year, 
                    y = value,
                    group = calibration_group,
                    color = calibration_group,
                    linetype = calibration_group))+
      scale_linetype_manual(values=c("dashed", "solid", "dashed", "solid"))+
      scale_color_manual(values=c('darkgreen','darkgreen', 'purple', 'purple'))+
      ggtitle(paste0('hiv_prevalence_sim_id_', sim_id_itr))+
      ylab('HIV_prev_per_100K')
    
    png(file_name_hiv_prev_graph, width = 480, height = 480)
    print(hiv_prev_graph)
    dev.off()
    
    setwd(outdir_best_analysis_ART)
    file_name_ART_coverage_graph<-paste0('ART_coverage_sim_id_', sim_id_itr, '.png')
    
    
    ART_coverage_df<-state_prog_df%>%
      filter(sim_id == sim_id_itr)%>%
      filter(HIV_compartment %in% c(2,3,4))%>%
      mutate(temp_ART = if_else(HIV_compartment == 4, value, 0))%>%
      group_by(year, time, G_compartment)%>%
      summarise(total_prev = sum(value),
                total_art = sum(temp_ART))%>%
      mutate(ART_coverage = total_art/total_prev)%>%
      group_by(year, G_compartment)%>%
      summarise(ART_coverage_TB_model_estimates = median(ART_coverage))%>%
      mutate(G_compartment = as.integer(G_compartment))%>%
      left_join(HIV_ART_cov_df, by = c('year', 'G_compartment'))
    
    ART_coverage_df_melt<-reshape2::melt(ART_coverage_df, 
                                      id.vars = c("year",
                                                  'G_compartment'))
    
    ART_coverage_df_melt$calibration_group<-(paste0('Gender_', ART_coverage_df_melt$G_compartment, '_',
                                                    ART_coverage_df_melt$variable))
    
    ART_coverage_graph<-ggplot()+
      geom_line(data = ART_coverage_df_melt%>%filter(year < 2027),
                aes(x = year, 
                    y = value,
                    group = calibration_group,
                    color = calibration_group,
                    linetype = calibration_group))+
      scale_linetype_manual(values=c("dashed", "solid", "dashed", "solid"))+
      scale_color_manual(values=c('darkgreen','darkgreen', 'purple', 'purple'))+
      ggtitle(paste0('ART_coverage_sim_id_', sim_id_itr))+
      ylab('ART coverage (percent)')
    
    png(file_name_ART_coverage_graph, width = 480, height = 480)
    print(ART_coverage_graph)
    dev.off()
    
  }
  
  ###mort calcs###
  setwd(outdir_best_analysis_TB_MORT)
  mort_calibration_df<-state_prog_df%>%
    filter(TB_compartment == 'mort')%>%
    mutate(calibration_group = paste0(DR_compartment, '_', HIV_compartment, '_', G_compartment))%>%
    select(c('year', 'sim_id', 'time', 'calibration_group', 'value'))%>%
    group_by(year, sim_id, calibration_group)%>%
    summarise(cum_mort = max(value))%>%
    as.data.frame()
  
  #mort_calibration_df<-reshape2::melt(mort_calibration_df, id = c('year', 'sim_id'))
  mort_calibration_df<-mort_calibration_df%>%
    group_by(calibration_group, sim_id)%>%
    mutate(mort_in_year = cum_mort - lag(cum_mort))%>%
    ungroup()%>%
    mutate(mort_in_year = if_else(year == 1990, cum_mort, mort_in_year))%>%
    filter(year < 2028)%>%
    mutate(G_compartment = if_else(grepl('female', calibration_group), '2', '1'))%>%
    left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()
  
  mort_calibration_df_all<-mort_calibration_df[,c('sim_id',
                                                  'year', 
                                                  'calibration_group',
                                                  'G_compartment',
                                                  'mort_rate_per_100K')]
  
  for (s in mse_df_top$sim_id){
    
    mort_calibration_df_all_temp<-mort_calibration_df_all%>%
      filter(sim_id == s)
    
    #graph most promising policies
    for (g in unique(mort_calibration_df_all_temp$G_compartment)){
      
      file_name <- paste0('TB_moratality_calibration_g_', g, '_simid_', s, '.png')
      
      df_temp <- mort_calibration_df_all_temp%>%filter(G_compartment == g)
      df_temp_model<-df_temp
      df_temp_model$calibration_group<-paste0('Model--', 
                                              df_temp_model$calibration_group)
      #df_temp_GBD<-GBD_calibration_df%>%filter(G_compartment == g)
      #df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--', 
      #                                      df_temp_GBD$calibration_group)
      #df1<-df_temp_GBD%>%filter(calibration_group == unique(df_temp_GBD$calibration_group)[2])
      #df2<-df_temp_GBD%>%filter(calibration_group == unique(df_temp_GBD$calibration_group)[1])
    
      gender_temp<-if_else(g == 1, 'Males', 'Females')
      
      graph_temp <- ggplot()+
        #geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), 
        #            inherit.aes = FALSE,fill = "plum2")+
        #geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), 
        #            inherit.aes = FALSE,fill = "darkseagreen1")+
        geom_line(data = df_temp_model, aes(x = year, y = mort_rate_per_100K, 
                                            color = calibration_group),
                  linetype="dashed", size = 1.2)+
        #geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate, 
        #                                  color = calibration_group))+
        scale_color_manual(values=c('darkorchid4', 'darkgreen'))+
        #scale_color_manual(values=c('purple', 'green', 'darkorchid4', 'darkgreen'))+
        labs(title = paste0('Deaths, rate per 100K, ', gender_temp))+
        scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
        scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
        theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      png(file_name, width = 480, height = 480)
      print(graph_temp)
      dev.off()
    }
  }
}

#find best based on TB mort and HIV prev
setwd(paste0(here(),
             "/param_files/calibration_code_results/2021_07_21/",
             "analysis/"))

best_mse_df<-read.csv('mse_df_top.csv')

outdir_best_analysis<-paste0(outdir_analysis, '/tb_mort_hiv_prev_best')
outdir_best_analysis_ART<-paste0(outdir_best_analysis, '/ART_COVERAGE')
outdir_best_analysis_HIV_PREV<-paste0(outdir_best_analysis, '/HIV_PREV')
outdir_best_analysis_TB_MORT<-paste0(outdir_best_analysis, '/TB_MORT')

dir.create(file.path(outdir_best_analysis))
dir.create(file.path(outdir_best_analysis_ART))
dir.create(file.path(outdir_best_analysis_HIV_PREV))
dir.create(file.path(outdir_best_analysis_TB_MORT))

top_graphs_func(best_mse_df)

