#Calibration Analysis May 21

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
current_date <- gsub('-', '_', Sys.Date())
warmup_date<-'2021_05_24'
indir<-paste0(here(),'/model_outputs/calibration/', current_date)
indir_sim_ref<-paste0(here(),'/model_outputs/calibration_data_sets/', warmup_date)
outdir_analysis <- paste0(here(),'/model_outputs/calibration/', current_date, '/analysis')

#will give warning if wd already exists
#ignore warning!
dir.create(file.path(outdir_analysis))

file_list <- list.files(path=indir)
file_list<-file_list[file_list != "analysis"]

#######Create calibration ref df#########
setwd(indir_sim_ref)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')

# #####Graphs that describe model states overtime####

sim_ids<-1:length(file_list)

#subset list so the analysis will run
subset_list<-seq(from = 244, to = length(file_list)+243, by = 243)

setwd(paste0(here(), '/param_files'))
GBD_calibration_df<-read.csv('calibration_rates_df.csv')

mse_df<-data.frame()

lapply(subset_list, function(s){
  
  #print(s)
  #print(Sys.time())
  upper_sim <-s
  lower_sim <-upper_sim-243
  
  calibration_df_all_temp<-data.frame()
  
  setwd(indir)
  for (i in lower_sim:upper_sim){
    print(i)
    temp_data <- fread(file_list[i]) 
    temp_data$sim_id <- rep(i, times = nrow(temp_data))
    
    calibration_df_all_temp <- rbind(calibration_df_all_temp, temp_data) #for each iteration, bind the new data to the building dataset
  }
  #print(Sys.time() - time)
  
  calibration_df_all_temp<-subset(calibration_df_all_temp, select = -V1)
  
  state_prog_df <-reshape2::melt(calibration_df_all_temp, 
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
    group_by(time, year, G_compartment, sim_id)%>%
    summarise(total_gender_pop = sum(value))%>%
    ungroup()%>%
    group_by(year, sim_id, G_compartment)%>%
    summarise(total_gender_pop = median(total_gender_pop))
  
  total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  mort_calibration_df<-calibration_df_all_temp%>%
    select(c('year', 'sim_id', 'time', 'HIV_neg_male', 
             'HIV_neg_female', 'HIV_pos_male', 'HIV_pos_female'))%>%
    group_by(year, sim_id)%>%
    summarise(cum_mort_hiv_neg_male = max(HIV_neg_male),
              cum_mort_hiv_neg_female= max(HIV_neg_female),
              cum_mort_hiv_pos_male= max(HIV_pos_male),
              cum_mort_hiv_pos_female= max(HIV_pos_female))%>%
    as.data.frame()
  
  mort_calibration_df<-reshape2::melt(mort_calibration_df, id = c('year', 'sim_id'))
  mort_calibration_df<-mort_calibration_df%>%
    group_by(variable, sim_id)%>%
    mutate(mort_in_year = value - lag(value))%>%
    ungroup()%>%
    mutate(mort_in_year = if_else(year == 1990, value, mort_in_year))%>%
    filter(year < 2018)%>%
    mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
    left_join(total_pop_in_gender_df%>%filter(sim_id < s)%>%
                filter(sim_id >= lower_sim), by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()
  
  mort_calibration_df_all<-mort_calibration_df[,c('year', 'variable', 'sim_id', 'mort_rate_per_100K')]
  
  calibration_group_name<-c()
  print(paste0('total mort rows', nrow(mort_calibration_df_all)))
  for (n in 1:nrow(mort_calibration_df_all)){
    print(n)
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
  
  mort_calibration_df_all$calibration_group<-calibration_group_name
  mort_calibration_df_all<-mort_calibration_df_all%>%
    rename(model_output = mort_rate_per_100K)%>%
    select(c('year', 'calibration_group', 'model_output', 'sim_id'))%>%
    left_join(GBD_calibration_df, by = c('year','calibration_group'))%>%
    filter(year < 2017)
  
  mse_df_temp<-mort_calibration_df_all%>%
    mutate(diff = (model_output - expected_rate),
           mse = diff^2)%>%
    group_by(sim_id, calibration_group)%>%
    summarise(calibration_group_mse = sum(mse))%>%
    ungroup()
  
  mse_df_temp<-reshape2::dcast(mse_df_temp, sim_id ~ calibration_group)
  mse_df_temp$total_mse <-rowSums(mse_df_temp[,2:5])
  
  setwd(outdir_analysis)
  write.csv(mse_df_temp, paste0('mse_df_upper_', s, '.csv'))
  mse_df<<-rbind(mse_df, mse_df_temp)
  
})

setwd(outdir_analysis)
mse_df<-mse_df%>%
  mutate(rank_total = rank(total_mse),
         rank_1 = rank(`HIV/TB_coinfection_Male`),
         rank_2 = rank(`HIV/TB_coinfection_Female`),
         rank_3 = rank(TB_only_Male),
         rank_4 = rank(TB_only_Female))%>%
  left_join(sim_calibration_ref_df%>%
              rename(sim_id = sim_calib_id), 
            by = c('sim_id'))%>%
  as.data.frame()

write.csv(mse_df, 'mse_df.csv')

  
top_graphs_func<-function(s_top){
  
  calibration_df_all_top<-data.frame()
  
  setwd(indir)
  for (i in s_top){
    file_name_temp<-paste0('out_df_sim_id_', i,'.csv')
    temp_data <- fread(file_name_temp) 
    calibration_df_all_top <- rbind(calibration_df_all_top, temp_data) #for each iteration, bind the new data to the building dataset
  }
  calibration_df_all_top<-subset(calibration_df_all_top, select = -V1)
  
  setwd(outdir_analysis)
  
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
  
  state_prog_df<-state_prog_df%>%
    filter(G_compartment %in% c(1,2))
  
  total_pop_in_gender_df<-state_prog_df%>%
    group_by(time, year, G_compartment, sim_id)%>%
    summarise(total_gender_pop = sum(value))%>%
    ungroup()%>%
    group_by(year, sim_id, G_compartment)%>%
    summarise(total_gender_pop = median(total_gender_pop))%>%
    filter(G_compartment %in% c(1,2))
  
  total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  mort_calibration_df<-calibration_df_all_top%>%
    select(c('year', 'sim_id', 'time', 'HIV_neg_male', 
             'HIV_neg_female', 'HIV_pos_male', 'HIV_pos_female'))%>%
    group_by(year, sim_id)%>%
    summarise(cum_mort_hiv_neg_male = max(HIV_neg_male),
              cum_mort_hiv_neg_female= max(HIV_neg_female),
              cum_mort_hiv_pos_male= max(HIV_pos_male),
              cum_mort_hiv_pos_female= max(HIV_pos_female))%>%
    as.data.frame()
  
  mort_calibration_df<-reshape2::melt(mort_calibration_df, id = c('year', 'sim_id'))
  mort_calibration_df<-mort_calibration_df%>%
    group_by(variable, sim_id)%>%
    mutate(mort_in_year = value - lag(value))%>%
    ungroup()%>%
    mutate(mort_in_year = if_else(year == 1990, value, mort_in_year))%>%
    filter(year < 2018)%>%
    mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
    left_join(total_pop_in_gender_df%>%filter(sim_id %in% s_top),
              by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()
  
  mort_calibration_df_all<-mort_calibration_df[,c('year', 'variable', 'sim_id', 'mort_rate_per_100K')]
  
  calibration_group_name<-c()
  #print(paste0('total mort rows', nrow(mort_calibration_df_all)))
  for (n in 1:nrow(mort_calibration_df_all)){
    #print(n)
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
  
  mort_calibration_df_all$calibration_group<-calibration_group_name
  mort_calibration_df_all<-mort_calibration_df_all%>%
    rename(model_output = mort_rate_per_100K)%>%
    select(c('year', 'calibration_group', 'model_output', 'sim_id'))%>%
    left_join(GBD_calibration_df, by = c('year','calibration_group'))%>%
    filter(year < 2017)
  
  for (s in s_top){
    
    mort_calibration_df_all_temp<-mort_calibration_df_all%>%
      filter(sim_id == s)
    
    print(mort_calibration_df_all_temp)
    
    #graph most promising policies
    for (g in unique(mort_calibration_df_all_temp$sex_name)){
      
      file_name <- paste0('TB_moratality_calibration', g, '_simid_', s, '.png')
      
      df_temp <- mort_calibration_df_all_temp%>%filter(sex_name == g)
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
        labs(title = paste0('Deaths, rate per 100K, ', g, '\n sim id', s))+
        scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
        scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
        theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      
      print(getwd())
      png(file_name, width = 480, height = 480)
      print(graph_temp)
      dev.off()
    }
  }
}

top_graphs_func(unique(mse_df_top$sim_id))

write.csv(mse_df_top, 'best_calabs_df.csv')

