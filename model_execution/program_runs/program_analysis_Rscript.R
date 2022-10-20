#Program Analysis
#generate graphs
#calib graphs per 100K males females
#program graphs per 100K males females

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'extrafont', 'RColorBrewer'), require, character.only=T)

##GITHUB/LOCAL###
#location where we sent outputs
library(here)
indir_outputs<-paste0(here(), '/test/program_outputs')
indir_target_calibration_estimates<-paste0(here(), '/param_files/target_calibration_estimates')
outdir_graphs<-paste0(here(), '/test/program_graphs')

setwd(indir_target_calibration_estimates)
HIV_prev_df<-read.csv('KZN_GBD_HIV_prev_rate_calibration_df.csv')%>%
  mutate(sex = as.character(sex))
TB_inc_df<-read.csv('KZN_GBD_TB_inc_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))
TB_mort_df<-read.csv('KZN_GBD_TB_mort_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))

setwd(outdir_graphs)
graph_ref<-read.csv('graph_ref.csv')
graph_ref$variable<-as.character(graph_ref$variable)
graph_ref$measure<-as.character(graph_ref$measure)
graph_ref$sex<-as.character(graph_ref$sex)
graph_ref$hs<-as.character(graph_ref$hs)
graph_ref$TB_HIV_coinfection<-as.character(graph_ref$TB_HIV_coinfection)

##combine simulation estimates##
setwd(indir_outputs)
all_files_in_outdir <- list.files(pattern="*.csv")

outputs_combined_df<-read.csv(all_files_in_outdir[1])
for (i in 2:length(all_files_in_outdir)){
  temp<-read.csv(all_files_in_outdir[i])
  outputs_combined_df<-rbind(outputs_combined_df, temp)
}

#display.brewer.pal(n = 10, name = "Paired")
colors_for_graph <- brewer.pal(n = 10, name = "Paired")
color_lookup_range<-c(1,3,5)
color_lookup_expected<-c(2,4,6)

###calibration graphs####
for(current_var_eval in graph_ref$variable){
  
  outputs_df_select<-outputs_combined_df%>%
    filter(year <= 2017)%>%
    select(c('year',
             'sim_id',
             current_var_eval))
  
  graph_ref_temp<-graph_ref%>%filter(variable == current_var_eval)
  TB_HIV_coinfection_temp<-graph_ref_temp$TB_HIV_coinfection[1]
  sex_temp<-graph_ref_temp$sex[1]
  
  if (grepl('inc', current_var_eval)){
    GBD_df<-TB_inc_df%>%
      filter(TB_HIV_coinfection == TB_HIV_coinfection_temp,
             sex == sex_temp)%>%
      filter(year <= 2017)
  } else if (grepl('mort', current_var_eval)){
    GBD_df<-TB_mort_df%>%
      filter(TB_HIV_coinfection == TB_HIV_coinfection_temp,
             sex == sex_temp)%>%
      filter(year <= 2017)
  } else {
    GBD_df<-HIV_prev_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2017)
  }
  
  
}


###graphs with program differences###
for(current_var_eval in graph_ref$variable){
  outputs_df_select<-outputs_combined_df%>%
    select(c('year',
             'program_id',
             'sim_id',
             current_var_eval))
  
  outputs_df_select_melt<-melt(outputs_df_select, 
                               id = c('year',
                                      'program_id',
                                      'sim_id'))
  
  program_eval_summarised<-outputs_df_select_melt%>%
    group_by(year, program_id, variable)%>%
    summarise(max_rate = max(value),
              min_rate = min(value),
              expected_rate = mean(value))
  
  graph_temp <- ggplot()+
    geom_ribbon(data=program_eval_summarised%>%
                  filter(year <= 2018),
                aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE, 
                fill = 'grey75')+
    geom_line(data=program_eval_summarised%>%
                filter(year <= 2018), 
              aes(x = year, y = expected_rate))
  
  for (p in c(1,2,3)){
    df_temp<-program_eval_summarised%>%
      filter(year >= 2018,
             program_id == p)
    
    graph_temp<-graph_temp+
      geom_ribbon(data=df_temp, aes(x = year, ymin = min_rate, ymax = max_rate), 
                  inherit.aes = FALSE,
                  fill = colors_for_graph[color_lookup_range[p]])
    
  }
  
  graph_temp<-graph_temp+
    geom_line(data = program_eval_summarised%>%
                filter(year >= 2018), 
              aes(x = year, y = expected_rate, 
                  color = as.factor(program_id)), size = 1)+
    scale_color_manual(values=c(colors_for_graph[color_lookup_expected]))
  
  measure = if_else(grepl('inc', current_var_eval), 
                    'Tuberculosis incidence rate per 100K ', 
                    if_else(grepl('mort', current_var_eval),
                            'Tuberculosis mortality rate per 100K ',
                            'HIV prevalence rate per 100K '))
  
  measure2 = if_else(grepl('inc', current_var_eval), 
                    'TB incidence rate', 
                    if_else(grepl('mort', current_var_eval),
                            'TB mortality rate',
                            'HIV prevalence rate'))
  
  gender_temp = if_else(grepl('female', current_var_eval),
                        'Females', 'Males')
  
  hs = if_else(grepl('pos', current_var_eval), ' and HIV positive', 
               if_else(grepl('neg', current_var_eval), ' and HIV negative', ''))
  
  
  
  graph_temp<-graph_temp+
    labs(title = paste0(measure, gender_temp, hs))+
    scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
    scale_y_continuous(name = measure2)+#, limits = c(0, max_graph), breaks=(seq(0, max_graph, breaks_graph)))+
    theme(text = element_text(size=20, family="Times New Roman"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=20),
          plot.title = element_text(hjust = .5, size=22))+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+
    geom_vline(xintercept = 2018, linetype="dashed", 
               color = "darkgrey", size=1.5)+
    annotate("text", x=2013, y=10, label= "calibration period", size = 6)+
    annotate("text", x=2023, y=10, label= "evaluation period", size = 6)
  
  setwd(outdir_graphs)
  file_name = paste0(measure2, gender_temp, hs, '.png')
  png(file_name, width = 800, height = 600)
  print(graph_temp)
  dev.off()
}
