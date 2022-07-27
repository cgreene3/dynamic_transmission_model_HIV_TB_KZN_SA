#Policy Analysis Dec 7 to 2028

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops', 'data.table', 'extrafont'), require, character.only=T)

#############Set in directory and out directory###########
###########Make sure KZN_south_africa.Rproj is open, otherwise will need to change wd manually########
policy_eval_date<- '2021_12_07/'
indir_ref_data<-paste0(here(), '/param_files/policy_eval_params/')
indir_state_prog<-paste0(here(), '/policy_runs/results/', policy_eval_date)
outdir_analysis <- paste0(here(), '/policy_runs/results/', policy_eval_date, 'analysis')
dir.create(file.path(outdir_analysis))

#######Read in reference df for sim IDs tested#########
setwd(indir_ref_data)
mse_df_top<-read.csv('mse_df_top.csv')

#create 
outdir_analysis_TB_MORT<-paste0(outdir_analysis, '/TB_MORT')
outdir_analysis_HIV_MORT<-paste0(outdir_analysis, '/HIV_MORT')
outdir_analysis_TB_INCIDENCE<-paste0(outdir_analysis, '/TB_INCIDENCE')

dir.create(file.path(outdir_analysis_TB_MORT))
dir.create(file.path(outdir_analysis_HIV_MORT))
dir.create(file.path(outdir_analysis_TB_INCIDENCE))

#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
colors_for_graph = palette.colors(palette = "Paired")
color_lookup_range<-c(1,3,5)
color_lookup_expected<-c(2,4,6)

mort_start_indicators_df_4_cats<-data.frame()

#####Evaluating Simulations#####
####merge and visualize outputs
vis_TB_mort<-function(mse_df_top){
  
  policy_analysis_df_all_top<-data.frame()
  state_prog_df<-data.frame()
  total_pop_in_gender_df<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in mse_df_top$sim_id){
    
    file_name<-paste0('out_df_sim_id_', sim_id_itr, '.csv')
    print(sim_id_itr)
    temp_data <- fread(file_name) 
    temp_data$sim_id <- rep(sim_id_itr, times = nrow(temp_data))
    #for each iteration, bind the new data to the building dataset
    policy_analysis_df_all_top <- rbind(policy_analysis_df_all_top, temp_data) 
  }
  
  
  
  state_prog_df <-reshape2::melt(policy_analysis_df_all_top, 
                                 id.vars = c("time", "year",
                                             'sim_id', 'policy_id'))
  
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
    group_by(time, year, G_compartment, sim_id, policy_id)%>%
    summarise(total_gender_pop = sum(value))%>%
    ungroup()%>%
    group_by(year, sim_id, G_compartment)%>%
    summarise(total_gender_pop = median(total_gender_pop))
  
  total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  
  TB_mort_df <-policy_analysis_df_all_top%>%
    select(c("time", "year",
             'sim_id', 'policy_id',
             'TB_mort_HIV_neg_male', 'TB_mort_HIV_neg_female',
             'TB_mort_HIV_pos_male', 'TB_mort_HIV_pos_female'))
  
  TB_mort_df <-reshape2::melt(TB_mort_df, 
                                 id.vars = c("time", "year",
                                             'sim_id', 'policy_id'))
  
  TB_mort_df<-TB_mort_df%>%
    group_by(year, sim_id, variable, policy_id)%>%
    summarise(cum_mort = max(value))%>% #for 2017
    as.data.frame()
  
  #mort_calibration_df<-reshape2::melt(mort_calibration_df, id = c('year', 'sim_id'))
  TB_mort_df<-TB_mort_df%>%
    group_by(variable, sim_id, policy_id)%>%
    mutate(mort_in_year = cum_mort - lag(cum_mort))%>%
    ungroup()%>%
    mutate(mort_in_year = if_else(year == 1990, cum_mort, mort_in_year))%>%
    mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
    left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()%>%
    select(-c('cum_mort'))
  
  #pull indicators start of intervention
  mort_start_indicators_df_4_cats<<-TB_mort_df%>%
    filter(year == 2017,
           policy_id == 1)%>%
    group_by(variable)%>%
    summarise(avg = mean(mort_in_year),
              max = max(mort_in_year),
              min = min(mort_in_year))
  
  TB_mort_df$HIV_status = if_else(grepl('pos', TB_mort_df$variable),
                                           'positive', 'negative')
  
  
  
  setwd(outdir_analysis_TB_MORT)
  
  TB_mort_df_vis<-TB_mort_df%>%
    group_by(variable, G_compartment, year, HIV_status, policy_id)%>%
    summarise(max_rate = max(mort_rate_per_100K),
              min_rate = min(mort_rate_per_100K),
              expected_rate = mean(mort_rate_per_100K))%>%
    mutate(policy_name = paste0('Policy: ', policy_id))%>%
    filter(year <= 2028) #because does not really evaluate 2029
  
  for (g in c(1,2)){
    gender_temp<-if_else(g == 1, 'Males', 'Females')
    for (hs in c('negative', 'positive')){
      
      file_name <- paste0('TB_mortality_projections_by_policy_for_', 
                          gender_temp, hs, '.png')
      
      #first graph trends (calibration period)
      trend_df<-TB_mort_df_vis%>%
        filter(G_compartment == g, 
               year <= 2017,
               policy_id == 1, #all are the same before 2017
               HIV_status == hs)
      
      graph_temp <- ggplot()+
        geom_ribbon(data=trend_df, aes(x = year, ymin = min_rate, ymax = max_rate), 
                    inherit.aes = FALSE, fill = 'grey75')+
        geom_line(data=trend_df, aes(x = year, y = expected_rate))
      
      #then graph for each policy
      
      #first range
      for (p in c(1,2,3)){
        df_temp<-TB_mort_df_vis%>%
          filter(G_compartment == g, 
                 year >= 2017,
                 policy_id == p,
                 HIV_status == hs)
        
        graph_temp<-graph_temp+
          geom_ribbon(data=df_temp,aes(x = year, ymin = min_rate, ymax = max_rate), 
                      inherit.aes = FALSE,
                      fill = colors_for_graph[color_lookup_range[p]])
      }
      
      #then graph lines for each policy
      df_temp2<-TB_mort_df_vis%>%
        filter(G_compartment == g, 
                 year >= 2017,
                 HIV_status == hs)
      
      graph_temp<-graph_temp+
        geom_line(data = df_temp2, aes(x = year, y = expected_rate, 
                                      color = policy_name), size = 1)+
        scale_color_manual(values=c(colors_for_graph[color_lookup_expected]))
      
      
      max_graph = if_else(hs == 'negative', 150, 400)
      breaks_graph=if_else(hs == 'negative', 50, 90)
      
      graph_temp<-graph_temp+
        labs(title = paste0('Tuberculosis death rate per 100K ', gender_temp, ' and HIV ', hs))+
        scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
        scale_y_continuous(name = 'mortality rate', limits = c(0, max_graph), breaks=(seq(0, max_graph, breaks_graph)))+
        theme(text = element_text(size=20, family="Times New Roman"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.position="top", legend.title = element_blank(), legend.text=element_text(size=20),
              plot.title = element_text(hjust = .5, size=22))+
        guides(color=guide_legend(nrow=1,byrow=TRUE))+
        geom_vline(xintercept = 2017, linetype="dashed", 
                   color = "darkgrey", size=1.5)+
        annotate("text", x=2013, y=10, label= "calibration period", size = 6)+
        annotate("text", x=2021, y=10, label= "evaluation period", size = 6)
      
      png(file_name, width = 800, height = 600)
      print(graph_temp)
      dev.off()
    }
  }
}

vis_TB_mort(mse_df_top)
write.csv(mort_start_indicators_df_4_cats, 'mort_start_indicators_df_4_cats.csv')

incidence_start_indicators_df_4_cats<-data.frame()

vis_TB_incidence<-function(mse_df_top){
  
  policy_analysis_df_all_top<-data.frame()
  state_prog_df<-data.frame()
  total_pop_in_gender_df<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in mse_df_top$sim_id){
    
    file_name<-paste0('out_df_sim_id_', sim_id_itr, '.csv')
    print(sim_id_itr)
    temp_data <- fread(file_name) 
    temp_data$sim_id <- rep(sim_id_itr, times = nrow(temp_data))
    #for each iteration, bind the new data to the building dataset
    policy_analysis_df_all_top <- rbind(policy_analysis_df_all_top, temp_data) 
  }
  
  
  
  state_prog_df <-reshape2::melt(policy_analysis_df_all_top, 
                                 id.vars = c("time", "year",
                                             'sim_id', 'policy_id'))
  
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
  
  incidence_trend_df<<-state_prog_df%>%
    filter(time == 28,
           policy_id == 1)
  
  total_pop_in_gender_df<-state_prog_df%>%
    filter(G_compartment %in% c(1,2))%>%
    group_by(time, year, G_compartment, sim_id, policy_id)%>%
    summarise(total_gender_pop = sum(value))%>%
    ungroup()%>%
    group_by(year, sim_id, G_compartment)%>%
    summarise(total_gender_pop = median(total_gender_pop))
  
  total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  
  TB_incidence_df <-policy_analysis_df_all_top%>%
    select(c("time", "year",
             'sim_id', 'policy_id',
             'TB_incidence_HIV_neg_male', 'TB_incidence_HIV_neg_female',
             'TB_incidence_HIV_pos_male', 'TB_incidence_HIV_pos_female'))
  
  TB_incidence_df <-reshape2::melt(TB_incidence_df, 
                              id.vars = c("time", "year",
                                          'sim_id', 'policy_id'))
  
  TB_incidence_df<-TB_incidence_df%>%
    group_by(year, sim_id, variable, policy_id)%>%
    summarise(cum_incidence = max(value))%>% #for 2017
    as.data.frame()
  
  #mort_calibration_df<-reshape2::melt(mort_calibration_df, id = c('year', 'sim_id'))
  TB_incidence_df<-TB_incidence_df%>%
    group_by(variable, sim_id, policy_id)%>%
    mutate(incidence_in_year = cum_incidence - lag(cum_incidence))%>%
    ungroup()%>%
    mutate(incidence_in_year = if_else(year == 1990, cum_incidence, incidence_in_year))%>%
    mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
    left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(incidence_rate_per_100K = ((incidence_in_year*100000)/total_gender_pop))%>%
    as.data.frame()%>%
    select(-c('cum_incidence'))
  
  #pull indicators start of intervention
  incidence_start_indicators_df_4_cats<<-TB_incidence_df%>%
    filter(year == 2017,
           policy_id == 1)%>%
    group_by(variable)%>%
    summarise(avg = mean(incidence_in_year),
              max = max(incidence_in_year),
              min = min(incidence_in_year))
  
  TB_incidence_df$HIV_status = if_else(grepl('pos', TB_incidence_df$variable),
                                  'positive', 'negative')
  
  
  
  setwd(outdir_analysis_TB_INCIDENCE)
  
  TB_incidence_df_vis<-TB_incidence_df%>%
    group_by(variable, G_compartment, year, HIV_status, policy_id)%>%
    summarise(max_rate = max(incidence_rate_per_100K),
              min_rate = min(incidence_rate_per_100K),
              expected_rate = mean(incidence_rate_per_100K))%>%
    mutate(policy_name = paste0('Policy ', policy_id))%>%
    filter(year <= 2028) #because does not really evaluate 2029
  
  for (g in c(1,2)){
    gender_temp<-if_else(g == 1, 'Males', 'Females')
    for (hs in c('negative', 'positive')){
      
      file_name <- paste0('TB_incidence_projections_by_policy_for_', 
                          gender_temp, hs, '.png')
      
      #first graph trends (calibration period)
      trend_df<-TB_incidence_df_vis%>%
        filter(G_compartment == g, 
               year <= 2017,
               policy_id == 1, #all are the same before 2017
               HIV_status == hs)
      
      graph_temp <- ggplot()+
        geom_ribbon(data=trend_df, aes(x = year, ymin = min_rate, ymax = max_rate), 
                    inherit.aes = FALSE, fill = 'grey75')+
        geom_line(data=trend_df, aes(x = year, y = expected_rate))
      
      #then graph for each policy
      
      #first range
      for (p in c(1,2,3)){
        df_temp<-TB_incidence_df_vis%>%
          filter(G_compartment == g, 
                 year >= 2017,
                 policy_id == p,
                 HIV_status == hs)
        
        graph_temp<-graph_temp+
          geom_ribbon(data=df_temp,aes(x = year, ymin = min_rate, ymax = max_rate), 
                      inherit.aes = FALSE,
                      fill = colors_for_graph[color_lookup_range[p]])
      }
      
      #then graph lines for each policy
      df_temp2<-TB_incidence_df_vis%>%
        filter(G_compartment == g, 
               year >= 2017,
               HIV_status == hs)
      
      graph_temp<-graph_temp+
        geom_line(data = df_temp2, aes(x = year, y = expected_rate, 
                                       color = policy_name),#, linetype = policy_name), 
                  size = 1.2)+
        scale_color_manual(values=c(colors_for_graph[color_lookup_expected]))
      
      
      max_graph = if_else(hs == 'negative', 1500, 3500)
      breaks_graph=if_else(hs == 'negative', 500, 500)
      
      graph_temp<-graph_temp+
        labs(title = paste0('Tuberculosis incidence rate per 100K ', gender_temp, ' and HIV ', hs))+
        scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
        scale_y_continuous(name = 'TB incidence rate', limits = c(0, max_graph), breaks=(seq(0, max_graph, breaks_graph)))+
        theme(text = element_text(size=22, family="Times New Roman"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.position="top", legend.title = element_blank(), legend.text=element_text(size=24),
              plot.title = element_text(hjust = .5, size=24),
              axis.text.x = element_text(size = 22, family="Times New Roman"),
              axis.text.y = element_text(size = 22, family="Times New Roman"),
              axis.title.x = element_text(size = 22, family="Times New Roman"),
              axis.title.y = element_text(size = 22, family="Times New Roman"))+
        guides(color=guide_legend(nrow=1,byrow=TRUE))+
        geom_vline(xintercept = 2017, linetype="dashed", 
                   color = "darkgrey", size=1.5)+
        annotate("text", x=2012, y=10, label= "calibration period", size = 8, family="Times New Roman")+
        annotate("text", x=2022, y=10, label= "evaluation period", size = 8, family="Times New Roman")
      
      png(file_name, width = 800, height = 600)
      print(graph_temp)
      dev.off()
    }
  }
}

vis_TB_incidence(mse_df_top)
write.csv(incidence_start_indicators_df_4_cats, 'incidence_start_indicators_df_4_cats.csv')
