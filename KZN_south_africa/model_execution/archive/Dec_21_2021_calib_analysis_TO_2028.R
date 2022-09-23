#Calibration Analysis Jan 21 2022 to 2028

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
start_eval_date <- '2022_01_21/' 
indir_ref_data<-paste0(here(),'/param_files/')
indir_mse_dfs<-paste0(here(),'/param_files/calibration_code_results/', start_eval_date)

#set calib outputs to local file on computer (otherwise files will overwhelm github)
indir_state_prog<-paste0("~/Documents/academic_posttt_2020/HIV_TB/calib_outdir/", start_eval_date, 'to_2028')

outdir_analysis <- paste0(here(),'/param_files/calibration_code_results/', start_eval_date, 'analysis/1990_2028_ANALYSIS')
dir.create(file.path(outdir_analysis))

#######Create calibration ref df#########
setwd(indir_ref_data)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')
setwd(indir_mse_dfs)
mse_df_top<-read.csv('mse_df_top.csv')

setwd(indir_ref_data)
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

#####Evaluating Top MSE#####

####TOP GRAPHS FUNCTION#####
top_graphs_func<-function(mse_df_top){
  
  calibration_df_all_top<-data.frame()
  state_prog_df<-data.frame()
  total_pop_in_gender_df<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in mse_df_top$sim_id){
  
    file_name<-paste0('out_df_sim_id_', sim_id_itr, '.csv')
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
    
    
    
    ####ART COVERAGE###
    setwd(outdir_best_analysis_ART)
    file_name_ART_coverage_graph<-paste0('ART_coverage_sim_id_', sim_id_itr, '.png')
    
    past_2018_art_cov<-ART_df%>%
      filter(year == 2018)
    past_2018_art_cov<-past_2018_art_cov$art_coverage
    
    
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
      left_join(ART_df%>%select('year', 'G_compartment','art_coverage'), 
                by = c('year', 'G_compartment'))%>%
      mutate(art_coverage = if_else(year > 2018, past_2018_art_cov[G_compartment],
                                    art_coverage))
    
    ART_coverage_df_melt<-reshape2::melt(ART_coverage_df, 
                                      id.vars = c("year",
                                                  'G_compartment'))
    
    ART_coverage_df_melt$calibration_group<-(paste0('gender_', ART_coverage_df_melt$G_compartment, '_',
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
        labs(title = paste0('Tuberculosis death rate per 100K ', gender_temp))+
        scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2028, by = 5))+
        scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
        theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      png(file_name, width = 480, height = 480)
      print(graph_temp)
      dev.off()
    }
  }
}

outdir_best_analysis_ART<-paste0(outdir_analysis, '/ART_COVERAGE')
outdir_best_analysis_HIV_PREV<-paste0(outdir_analysis, '/HIV_PREV')
outdir_best_analysis_TB_MORT<-paste0(outdir_analysis, '/TB_MORT')
outdir_best_analysis_TB_MORT_GBD_OVERLAY<-paste0(outdir_analysis, '/TB_MORT_GBD_OVERLAY')

dir.create(file.path(outdir_best_analysis_ART))
dir.create(file.path(outdir_best_analysis_HIV_PREV))
dir.create(file.path(outdir_best_analysis_TB_MORT))
dir.create(file.path(outdir_best_analysis_TB_MORT_GBD_OVERLAY))

top_graphs_func(mse_df_top)

####merge and visualize outputs
vis_all_top_graphs_func<-function(mse_df_top){
  
  calibration_df_all_top<-data.frame()
  state_prog_df<-data.frame()
  total_pop_in_gender_df<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in mse_df_top$sim_id){
    
    file_name<-paste0('out_df_sim_id_', sim_id_itr, '.csv')
    print(sim_id_itr)
    temp_data <- fread(file_name) 
    temp_data$sim_id <- rep(sim_id_itr, times = nrow(temp_data))
    #for each iteration, bind the new data to the building dataset
    calibration_df_all_top <- rbind(calibration_df_all_top, temp_data) 
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
  
  #hiv prev overlay
  HIV_calibration_df_model_proj<-state_prog_df%>%
    filter(!TB_compartment %in% c('mort', 'incidence'),
           HIV_compartment > 1)%>%
    group_by(time, year, sim_id, G_compartment)%>%
    summarise(total_prev = sum(value))%>%
    group_by(sim_id, year, G_compartment)%>%
    summarise(total_prev = median(total_prev))%>%
    left_join(total_pop_in_gender_df, by = c('sim_id', 'year', 'G_compartment'))%>%
    mutate(HIV_prev_TB_model_estimates_per_100K = total_prev*(100000/total_gender_pop))%>%
    mutate(G_compartment = as.integer(G_compartment))%>%
    select(c('year', 'G_compartment', 'HIV_prev_TB_model_estimates_per_100K'))%>%
    mutate(categories_data_vis = paste0('TB-HIV model projections for ',
                                        if_else(G_compartment == 1, 'males', 'females')))
    
  hiv_prev_calibration_df_vis<-HIV_calibration_df_model_proj%>%
    group_by(categories_data_vis, year, G_compartment)%>%
    summarise(max_rate = max(HIV_prev_TB_model_estimates_per_100K)/100000,
              min_rate = min(HIV_prev_TB_model_estimates_per_100K)/100000,
              hiv_prevalence = mean(HIV_prev_TB_model_estimates_per_100K)/100000)
  
  #https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
  colors_for_hiv_graph = palette.colors(palette = "Accent")
  
  HIV_calibration_df$categories_data_vis <- paste0('HPV-HIV model projections for ', HIV_calibration_df$gender_name, 's')
  
  hiv_prev_all_data_vis <- rbind(hiv_prev_calibration_df_vis%>%
                                   select(c('year', 'categories_data_vis', 'hiv_prevalence')), 
                                 HIV_calibration_df%>%select(c('year', 
                                                               'categories_data_vis',
                                                               'hiv_prevalence')))
  
  
  hiv_prev_graph_temp <- ggplot()+
    geom_ribbon(data = hiv_prev_calibration_df_vis%>%
                  filter(G_compartment == 1), aes(x = year, ymin = min_rate, ymax = max_rate),
                inherit.aes = FALSE, fill = colors_for_hiv_graph[1])+
    geom_ribbon(data = hiv_prev_calibration_df_vis%>%
                  filter(G_compartment == 2), aes(x = year, ymin = min_rate, ymax = max_rate),
                inherit.aes = FALSE, fill = colors_for_hiv_graph[2])+
    geom_line(data = hiv_prev_all_data_vis,
              aes(x = year, 
                  y = hiv_prevalence,
                  group = categories_data_vis,
                  color = categories_data_vis,
                  linetype = categories_data_vis),
              size = .5)+
    scale_linetype_manual(values=c("dashed", "dashed", "solid", "solid"))+
    scale_color_manual(values=c('darkorchid1','forestgreen', 'darkorchid4', 'darkgreen'))+
    ggtitle(paste0('Projected HIV prevalence: HPV-HIV model projections\ncompared to TB-HIV model projections (average over best calibration parameters)'))+
    ylab('HIV prevalence')+
    scale_x_continuous(name = 'year', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
    scale_y_continuous(name = 'HIV prevalence', limits = c(0, 1))+
    theme(text = element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=10),
          plot.title = element_text(hjust = .5, size=12))+
    guides(linetype = guide_legend(nrow = 2, byrow = TRUE))+
    geom_vline(xintercept = 2017, linetype="dashed", 
               color = "darkgrey", size=1.5)+
    annotate("text", x=2013, y=.5, label= "calibration period", size = 4)+
    annotate("text", x=2021, y=.5, label= "evaluation period", size = 4)
  
  
  
  #mort overlay
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
    #filter(year < 2018)%>%
    mutate(G_compartment = if_else(grepl('female', calibration_group), '2', '1'))%>%
    left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()
  
  mort_calibration_df$HIV_status = if_else(grepl('pos', mort_calibration_df$calibration_group),
                                           'positive', 'negative')
  
  setwd(outdir_best_analysis_TB_MORT_GBD_OVERLAY)
  
  mort_calibration_df_vis<-mort_calibration_df%>%
    group_by(calibration_group, G_compartment, year, HIV_status)%>%
    summarise(max_rate = max(mort_rate_per_100K),
              min_rate = min(mort_rate_per_100K),
              expected_rate = mean(mort_rate_per_100K))
  
  for (g in c(1,2)){
    
    gender_temp<-if_else(g == 1, 'Males', 'Females')
    
    file_name <- paste0('TB_mortality_projections_for_', gender_temp, '.png')
    
    df_temp_model <- mort_calibration_df_vis%>%
      filter(G_compartment == g, year < 2028)
    
    df_neg_model<- df_temp_model%>%filter(HIV_status == 'negative')
    df_pos_model<- df_temp_model%>%filter(HIV_status == 'positive')
    
    df_temp_GBD<-GBD_calibration_df%>%
      filter(G_compartment == g)%>%
      mutate(HIV_status = if_else(grepl('neg', calibration_group),
                                  'negative', 'positive'))
    
    df_neg_GBD <-df_temp_GBD%>%filter(HIV_status == 'negative')
    df_pos_GBD <-df_temp_GBD%>%filter(HIV_status == 'positive')
    
    df_temp_expected_rate_model<-df_temp_model%>%
      select(c('year', 'HIV_status', 'expected_rate'))%>%
      mutate(calibration_group = paste0('Model Projections: HIV ', HIV_status))
    df_temp_expected_rate_GBD<-df_temp_GBD%>%
      select(c('year', 'HIV_status', 'expected_rate'))%>%
      mutate(calibration_group = paste0('GBD Projections: HIV ', HIV_status))
    
    df_temp_expected_rate_all_df<-rbind(df_temp_expected_rate_model, df_temp_expected_rate_GBD)
    
    #https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
    colors_for_graph = palette.colors(palette = "Paired")
    
    graph_temp <- ggplot()+
      geom_ribbon(data = df_neg_GBD, aes(x = year, ymin = min_rate, ymax = max_rate),
                  inherit.aes = FALSE, fill = colors_for_graph[1])+
      geom_ribbon(data = df_pos_GBD, aes(x = year, ymin = min_rate, ymax = max_rate),
                  inherit.aes = FALSE, fill = colors_for_graph[3])+
      geom_ribbon(data=df_neg_model,aes(x = year, ymin = min_rate, ymax = max_rate), 
                  inherit.aes = FALSE,fill = colors_for_graph[5])+
      geom_ribbon(data=df_pos_model,aes(x = year, ymin = min_rate, ymax = max_rate), 
                  inherit.aes = FALSE,fill = colors_for_graph[7])+
      geom_line(data = df_temp_expected_rate_all_df, aes(x = year, y = expected_rate, 
                                    color = calibration_group))+
      scale_color_manual(values=c(colors_for_graph[c(2,4,6,8)]))+
      labs(title = paste0('Tuberculosis death rate per 100K ', gender_temp))+
      scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
      scale_y_continuous(name = 'mortality rate', limits = c(0, 1000), breaks=(seq(0, 1000, 200)))+
      theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="top", legend.title = element_blank(), legend.text=element_text(size=14),
            plot.title = element_text(hjust = .5, size=22))+
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      geom_vline(xintercept = 2017, linetype="dashed", 
                 color = "darkgrey", size=1.5)+
      annotate("text", x=2013, y=450, label= "calibration period", size = 6)+
      annotate("text", x=2021, y=450, label= "evaluation period", size = 6)
      
    
    
    png(file_name, width = 800, height = 600)
    print(graph_temp)
    dev.off()
    
  }
  
}

vis_all_top_graphs_func(mse_df_top)

