#Calibration Analysis Jun8

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
start_eval_date <- '2021_06_08'
indir_ref_data<-paste0(here(),'/param_files/calib_ref_data/', start_eval_date)

#set calib outputs to local file on computer (otherwise files will overwhelm github)
indir_state_prog<-"~/Documents/academic_posttt_2020/HIV_TB/calib_outdir/jun9"
outdir_analysis <- paste0(here(),'/param_files/calib_ref_data/', start_eval_date, '/analysis')

dir.create(file.path(outdir_analysis))

#######Create calibration ref df#########
setwd(indir_ref_data)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')

# #####Graphs that describe model states overtime####

#subset list so the analysis will run
subset_list<-seq(from = 288, to = nrow(sim_calibration_ref_df)+1, by = 288)

setwd(paste0(here(), '/param_files'))
HIV_calibration_df<-read.csv('hiv_transmission_df.csv')
HIV_calibration_df<-HIV_calibration_df%>%
  mutate(G_compartment = if_else(gender == 'Females', 2, 1))%>%
  dplyr::select('year', 'G_compartment', 'hiv_prevalence')



GBD_calibration_df<-read.csv('calibration_rates_df.csv')
#match calibration group with model output calibration groups
GBD_calibration_df$calibration_group<-if_else(GBD_calibration_df$calibration_group == "TB_only_Female",
                                              "HIV_neg_female",
                                              if_else(GBD_calibration_df$calibration_group == "TB_only_Male",
                                                      "HIV_neg_male",
                                                      if_else(GBD_calibration_df$calibration_group == "HIV/TB_coinfection_Male",
                                                              "HIV_pos_male", "HIV_pos_female")))


GBD_calibration_df$G_compartment <-if_else(GBD_calibration_df$sex_name == 'male',
                                           1, 2)

mse_df<-data.frame()


lapply(subset_list, function(s){
  
  upper_sim <-s
  lower_sim <-upper_sim-288
  
  calibration_df_all_temp<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in lower_sim:upper_sim){
    if(sim_id_itr %in% sim_calibration_ref_df$sim_calib_id){
      file_name<-paste0('out_df_sim_id_', sim_id_itr, '.csv')
      print(sim_id_itr)
      temp_data <- fread(file_name) 
      temp_data$sim_id <- rep(sim_id_itr, times = nrow(temp_data))
      #for each iteration, bind the new data to the building dataset
      calibration_df_all_temp <- rbind(calibration_df_all_temp, temp_data) 
    }
  }
  
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
    filter(G_compartment %in% c(1,2))%>%
    group_by(time, year, G_compartment, sim_id)%>%
    summarise(total_gender_pop = sum(value))%>%
    ungroup()%>%
    group_by(year, sim_id, G_compartment)%>%
    summarise(total_gender_pop = median(total_gender_pop))
  
  total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
  ###hiv prev calcs###
  HIV_prev_df<-state_prog_df%>%
    filter(HIV_compartment %in% c(2,3,4))%>%
    group_by(year, time, G_compartment, sim_id)%>%
    summarise(total_prev = sum(value))%>%
    left_join(total_pop_in_gender_df, by = c('year', 'G_compartment', 'sim_id'))%>%
    mutate(HIV_prev_per_100K = total_prev*(100000/total_gender_pop))%>%
    group_by(year, G_compartment, sim_id)%>%
    summarise(HIV_prev_per_100K = mean(HIV_prev_per_100K))
  
  HIV_prev_df<-HIV_prev_df%>%
    mutate(G_compartment = as.integer(G_compartment))%>%
    left_join(HIV_calibration_df, by = c('year', 'G_compartment'))%>%
    mutate(hiv_proj_per_100K = hiv_prevalence*100000)%>%
    mutate(hiv_prev_perc_diff = (HIV_prev_per_100K - hiv_proj_per_100K)^2,
           G_compartment = as.integer(G_compartment))%>%
  filter(year < 2018)%>%
    group_by(sim_id)%>%
    summarise(hiv_prev_mse_male = sum(if_else(G_compartment == 1, 
                                              hiv_prev_perc_diff, 0)),
              hiv_prev_mse_female = sum(if_else(G_compartment == 2, 
                                              hiv_prev_perc_diff, 0)),
              hiv_prev_mse = hiv_prev_mse_male+hiv_prev_mse_female)
  
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
    filter(year < 2018)%>%
    mutate(G_compartment = if_else(grepl('female', calibration_group), '2', '1'))%>%
    left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()
  
  mort_calibration_df_all<-mort_calibration_df[,c('sim_id',
                                                  'year', 
                                                  'calibration_group', 
                                                  'mort_rate_per_100K')]
  
  mort_calibration_df_all<-mort_calibration_df_all%>%
    left_join(GBD_calibration_df, 
              by = c('year','calibration_group'))
  
  mse_df_temp<-mort_calibration_df_all%>%
    mutate(diff = (mort_rate_per_100K - expected_rate),
           mse = diff^2)%>%
    group_by(sim_id, calibration_group)%>%
    summarise(calibration_group_mse = sum(mse))%>%
    ungroup()
    
  
  mse_df_temp<-reshape2::dcast(mse_df_temp, sim_id ~ calibration_group)
  mse_df_temp$tb_mort_mse_total <-rowSums(mse_df_temp[,2:5])
  
  mse_df_temp<-mse_df_temp%>%
    left_join(HIV_prev_df, by = c('sim_id'))
  
  mse_df<<-rbind(mse_df, mse_df_temp)
  
})

mse_df<-mse_df[!duplicated(mse_df), ]


setwd(outdir_analysis)
mse_df2<-mse_df%>%
  left_join(sim_calibration_ref_df%>%
              rename(sim_id = sim_calib_id), 
            by = c('sim_id'))%>%
  mutate(total_mse = tb_mort_mse_total+hiv_prev_mse)%>%
  mutate(tbmort_rank_total = rank(tb_mort_mse_total),
         tbmort_HIV_neg_male = rank(HIV_neg_male),
         tbmort_HIV_neg_female = rank(HIV_neg_female),
         tbmort_HIV_pos_male = rank(HIV_pos_male),
         tbmort_HIV_pos_female = rank(HIV_pos_female),
         hivprev_rank_total = rank(hiv_prev_mse),
         hivprev_male = rank(hiv_prev_mse_male),
         hivprev_female = rank(hiv_prev_mse_female),
         prev_mort_mse_rank = rank(total_mse))%>%
  as.data.frame()

write.csv(mse_df2, 'mse_df.csv')

mse_df_top<-mse_df2%>%
  filter(tbmort_rank_total <=3)
    #tbmort_rank_total <= 10|
          # hivprev_rank_total <= 10)#|
           #rank_1 <= 3|
           #rank_2 <= 3|
           #rank_3 <= 3|
           #rank_4 <= 3)

write.csv(mse_df_top, 'mse_df_top.csv')
  
top_graphs_func<-function(s_top){
  
  calibration_df_all_top<-data.frame()
  state_prog_df<-data.frame()
  total_pop_in_gender_df<-data.frame()
  
  setwd(indir_state_prog)
  
  for (sim_id_itr in s_top){
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
  
  ###mort calcs###
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
    filter(year < 2018)%>%
    mutate(G_compartment = if_else(grepl('female', calibration_group), '2', '1'))%>%
    left_join(total_pop_in_gender_df, by = c('G_compartment', 'year', 'sim_id'))%>%
    mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))%>%
    as.data.frame()
  
  mort_calibration_df_all<-mort_calibration_df[,c('sim_id',
                                                  'year', 
                                                  'calibration_group',
                                                  'G_compartment',
                                                  'mort_rate_per_100K')]
  
  setwd(outdir_analysis)
  
  for (s in s_top){
    
    mort_calibration_df_all_temp<-mort_calibration_df_all%>%
      filter(sim_id == s)
    
    #graph most promising policies
    for (g in unique(mort_calibration_df_all_temp$G_compartment)){
      
      file_name <- paste0('TB_moratality_calibration_g_', g, '_simid_', s, '.png')
      
      df_temp <- mort_calibration_df_all_temp%>%filter(G_compartment == g)
      df_temp_model<-df_temp
      df_temp_model$calibration_group<-paste0('Model--', 
                                              df_temp_model$calibration_group)
      df_temp_GBD<-GBD_calibration_df%>%filter(G_compartment == g)
      df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--', 
                                            df_temp_GBD$calibration_group)
      df1<-df_temp_GBD%>%filter(calibration_group == unique(df_temp_GBD$calibration_group)[2])
      df2<-df_temp_GBD%>%filter(calibration_group == unique(df_temp_GBD$calibration_group)[1])
    
      gender_temp<-if_else(g == 1, 'Males', 'Females')
      
      graph_temp <- ggplot()+
        geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), 
                    inherit.aes = FALSE,fill = "plum2")+
        geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), 
                    inherit.aes = FALSE,fill = "darkseagreen1")+
        geom_line(data = df_temp_model, aes(x = year, y = mort_rate_per_100K, 
                                            color = calibration_group),
                  linetype="dashed", size = 1.2)+
        geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate, 
                                          color = calibration_group))+
        scale_color_manual(values=c('purple', 'green', 'darkorchid4', 'darkgreen'))+
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

top_graphs_func(unique(mse_df_top$sim_id))

write.csv(mse_df_top, 'mse_df_top.csv', row.names = FALSE)

#test prevalence
setwd(indir_state_prog)
sim_id_test<-mse_df_top$sim_id[3]
file_name<-paste0('out_df_sim_id_', sim_id_test, '.csv')
best_calib_df <- fread(file_name) 

state_prog_df <-reshape2::melt(best_calib_df, 
                               id.vars = c("time", "year", 'sim_id'))

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
  group_by(time, year, G_compartment)%>%
  summarise(total_gender_pop = sum(value))%>%
  ungroup()%>%
  group_by(year, G_compartment)%>%
  summarise(total_gender_pop = median(total_gender_pop))

total_pop_in_gender_df<-as.data.frame(total_pop_in_gender_df)
  
TB_prev_df<-state_prog_df%>%
  filter(TB_compartment == 6)%>%
  #mutate(HIV_group = if_else(HIV_compartment == 1, 'HIV_negative', 'HIV_positive'))%>%
  group_by(time, G_compartment)%>%
  summarise(TB_prev_per_100K = sum(value)*2)

graph_temp<-TB_prev_df%>%
  ggplot(aes(x = time+1990, y = TB_prev_per_100K,
             group = G_compartment, 
             color = G_compartment))+
  geom_line()


HIV_prev_df<-state_prog_df%>%
  filter(HIV_compartment %in% c(2,3,4))%>%
  group_by(year, time, G_compartment)%>%
  summarise(total_prev = sum(value))%>%
  left_join(total_pop_in_gender_df, by = c('year', 'G_compartment'))%>%
  mutate(HIV_prev_per_100K = total_prev*(100000/total_gender_pop))

hiv_prev_graph<-HIV_prev_df%>%
  ggplot(aes(x = time+1990, y = HIV_prev_per_100K,
             group = G_compartment, 
             color = G_compartment))+
  geom_line()+
  ggtitle(paste0('sim id ', sim_id_test))

