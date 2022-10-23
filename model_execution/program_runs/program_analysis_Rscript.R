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
#indir_outputs<-paste0(here(), '/test/program_outputs')
indir_outputs<-paste0(here(), '/program_outputs')
indir_target_calibration_estimates<-paste0(here(), '/param_files/target_calibration_estimates')
outdir_graphs_prog<-paste0(here(), '/graphs_results/program')
outdir_graphs_calib<-paste0(here(), '/graphs_results/calibration')

setwd(indir_target_calibration_estimates)
HIV_prev_df<-read.csv('KZN_GBD_HIV_prev_rate_calibration_df.csv')%>%
  mutate(sex = as.character(sex))
TB_inc_df<-read.csv('KZN_GBD_TB_inc_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))
TB_mort_df<-read.csv('KZN_GBD_TB_mort_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))

setwd(paste0(here(), '/graphs_results'))
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

display.brewer.pal(n = 10, name = "Paired")
colors_for_graphs <- brewer.pal(n = 10, name = "Paired")

colors_for_calib_graph_fill_tb<-colors_for_graphs[c(1, 3, 5, 7)]
colors_for_calib_graph_line_tb<-colors_for_graphs[c(2, 4, 6, 8)]

colors_for_calib_graph_fill_hiv<-colors_for_graphs[c(7, 9)]
colors_for_calib_graph_line_hiv<-colors_for_graphs[c(8, 10)]

###calibration graphs####
outputs_combined_df_calibration<-outputs_combined_df%>%
  filter(program_id == 1)%>%
  select(-c('program_id'))
outputs_combined_df_calibration<-melt(outputs_combined_df_calibration, id = c("year", "sim_id"))
outputs_combined_df_calibration<-outputs_combined_df_calibration%>%
  left_join(graph_ref, by = c('variable'))

outputs_combined_df_calibration$title<-paste0(outputs_combined_df_calibration$measure, 
                                              outputs_combined_df_calibration$sex, "s")

for(current_measure_eval in unique(outputs_combined_df_calibration$title)){
  
  outputs_df_select<-outputs_combined_df_calibration%>%
    filter(title == current_measure_eval)
  
  if(!grepl('hiv', as.character(current_measure_eval), ignore.case = TRUE)){
    outputs_df_summarised<-outputs_df_select%>%
      group_by(year, hs)%>%
      summarise(upper_rate = max(value),
                lower_rate = min(value),
                val_rate = mean(value))%>%
      mutate(group = paste0('Model Projections: ', hs))%>%
      select(-c('hs'))
  } else {
    outputs_df_summarised<-outputs_df_select%>%
      group_by(year, hs)%>%
      summarise(upper_rate = max(value),
                lower_rate = min(value),
                val_rate = mean(value))%>%
      select(-c('hs'))
    outputs_df_summarised$group = "Model Projections"
  }
  
  
  sex_temp<-unique(outputs_df_select$sex)
  
  if (grepl('inc', current_measure_eval)){
    GBD_df<-TB_inc_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2018)
    measure_temp<-"TB incidence rate"
    max_min_graph_y_axis<-c(0, 2400, 400)
  } else if (grepl('mort', current_measure_eval)){
    GBD_df<-TB_mort_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2018)
    measure_temp<-"TB mortality rate"
    max_min_graph_y_axis<-c(0, 1000, 100)
  } else {
    GBD_df<-HIV_prev_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2018)
    measure_temp<-"HIV prevalence rate"
    max_min_graph_y_axis<-c(0, 50000, 10000)
  }
  
  if(!grepl('hiv', as.character(current_measure_eval), ignore.case = TRUE)) {
    GBD_df$group = paste0("GBD Projections: ", 
                          if_else(GBD_df$TB_HIV_coinfection == "no", 
                                  "HIV negative", "HIV positive"))
    } else {
      GBD_df$group <- "GBD Projections"
  }
  
  calibration_graph_df<-rbind(data.frame(GBD_df%>%select(year, val_rate, upper_rate, lower_rate, group)),
                              data.frame(outputs_df_summarised))
  
  if(!grepl('hiv', as.character(current_measure_eval), ignore.case = TRUE)){
    calib_graph_temp<-ggplot(calibration_graph_df)+
      geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
      geom_line(aes(x = year, y = val_rate, colour = group, linetype = group), size = .7)+
      theme(text = element_text(size=18, family="Times New Roman"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="top", legend.title = element_blank(), legend.text=element_text(size=15),
            plot.title = element_text(hjust = .5, size=20))+
      geom_vline(xintercept = 2018, linetype="dashed", 
                 color = "darkgrey", size=1.5)+
      annotate("text", x=2013, y=max_min_graph_y_axis[2]*.9, label= "calibration period", size = 6)+
      annotate("text", x=2023, y=max_min_graph_y_axis[2]*.9, label= "evaluation period", size = 6)+
      ggtitle(current_measure_eval)+
      scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
      #scale_y_continuous(name = measure_temp, breaks = c(seq(from = max_min_graph_y_axis[1],
      #                                                       to = max_min_graph_y_axis[2],
      #                                                       by = max_min_graph_y_axis[3])))+
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      scale_color_manual(values=colors_for_calib_graph_line_tb)+
      scale_fill_manual(values=colors_for_calib_graph_fill_tb)+
      ylim(max_min_graph_y_axis[1:2])+
      ylab(measure_temp)
  } else {
    calib_graph_temp<-ggplot(calibration_graph_df)+
      geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
      geom_line(aes(x = year, y = val_rate, colour = group, linetype = group), size = .7)+
      theme(text = element_text(size=20, family="Times New Roman"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="top", legend.title = element_blank(), legend.text=element_text(size=15),
            plot.title = element_text(hjust = .5, size=18))+
      geom_vline(xintercept = 2018, linetype="dashed", 
                 color = "darkgrey", size=1.5)+
      annotate("text", x=2013, y=max_min_graph_y_axis[2]*.9, label= "calibration period", size = 6)+
      annotate("text", x=2023, y=max_min_graph_y_axis[2]*.9, label= "evaluation period", size = 6)+
      ggtitle(current_measure_eval)+
      scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
      guides(color=guide_legend(nrow=1,byrow=TRUE))+
      scale_color_manual(values=colors_for_calib_graph_line_hiv)+
      scale_fill_manual(values=colors_for_calib_graph_fill_hiv)+
      ylim(max_min_graph_y_axis[1:2])+
      ylab(measure_temp)
  }
  
  file_name<-paste0("calibration_",
                    str_replace_all(current_measure_eval, " ", "_"),
                    ".png")
                    
  
  setwd(outdir_graphs_calib)
  png(file_name, width = 800, height = 600)
  print(calib_graph_temp)
  dev.off()
}
  

###graphs with program differences###
colors_for_program_graph_fill<-colors_for_graphs[c(1, 3, 5)]
colors_for_program_graph_line<-colors_for_graphs[c(2, 4, 6)]


outputs_combined_df_program<-melt(outputs_combined_df, 
                                      id = c("year", "sim_id", "program_id"))
outputs_combined_df_program<-outputs_combined_df_program%>%
  left_join(graph_ref, by = c('variable'))%>%
  filter(grepl('Tuberculosis', measure))

outputs_combined_df_program$title<-paste0(outputs_combined_df_program$measure, 
                                          outputs_combined_df_program$sex, 
                                          "s and ",
                                          outputs_combined_df_program$hs)


for(current_var_eval in unique(outputs_combined_df_program$title)){
  
  outputs_df_select<-outputs_combined_df_program%>%
    filter(title == current_var_eval)%>%
    mutate(group = paste0('Program ', program_id))
  
  program_eval_summarised<-outputs_df_select%>%
    group_by(year, group)%>%
    summarise(upper_rate = max(value),
              lower_rate = min(value),
              val_rate = mean(value))
  
  program_eval_summarised_connect<-program_eval_summarised%>%
    filter(year == 2018)
  
  program_eval_summarised$group<-if_else(program_eval_summarised$year <= 2018,
                                         "Calibration Period", program_eval_summarised$group)
  
  program_eval_summarised<-rbind(program_eval_summarised, program_eval_summarised_connect)
  
  if (grepl('inc', current_var_eval)){
    measure_temp<-"TB incidence rate"
    if(grepl('negative', current_var_eval)){
      max_min_graph_y_axis<-c(0, 600, 100)
    } else {
      max_min_graph_y_axis<-c(0, 2400, 400)
    }
  } else if (grepl('mort', current_var_eval)){
    measure_temp<-"TB mortality rate"
    if(grepl('negative', current_var_eval)){
      max_min_graph_y_axis<-c(0, 150, 20)
      } else{
      max_min_graph_y_axis<-c(0, 600, 150)
    }
  }
  
  program_graph_temp<-ggplot(program_eval_summarised)+
    geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
    geom_line(aes(x = year, y = val_rate, colour = group, linetype = group), size = .7)+
    theme(text = element_text(size=18, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=15),
          plot.title = element_text(hjust = .5, size=20))+
    geom_vline(xintercept = 2018, linetype="dashed", 
               color = "darkgrey", size=1.5)+
    annotate("text", x=2013, y=1, label= "calibration period", size = 6)+
    annotate("text", x=2023, y=1, label= "evaluation period", size = 6)+
    ggtitle(current_var_eval)+
    scale_x_continuous(name = 'time', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
    scale_color_manual(values=c("black", colors_for_program_graph_line))+
    scale_fill_manual(values=c('grey75', colors_for_program_graph_fill))+
    ylim(max_min_graph_y_axis[1:2])+
    ylab(measure_temp)
  
  file_name<-paste0(str_replace_all(current_var_eval, " ", "_"),".png")
  
  setwd(outdir_graphs_prog)
  png(file_name, width = 800, height = 600)
  print(program_graph_temp)
  dev.off()
}
