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
indir_outputs<-paste0(here(), '/results/program_outputs')
indir_target_calibration_estimates<-paste0(here(), '/param_files/target_calibration_estimates')
outdir_stats<-paste0(here(), '/results/stats')
outdir_graphs_prog<-paste0(here(), '/results/graphs/program')
outdir_graphs_prog2<-paste0(here(), '/results/graphs/program2')
outdir_graphs_prog3<-paste0(here(), '/results/graphs/program3')
outdir_graphs_prog4<-paste0(here(), '/results/graphs/program4')
outdir_graphs_calib<-paste0(here(), '/results/graphs/calibration')
outdir_graphs_calib2<-paste0(here(), '/results/graphs/calibration2')

setwd(indir_target_calibration_estimates)


HIV_prev_df<-read.csv('KZN_GBD_HIV_prev_rate_calibration_df.csv')%>%
  mutate(sex = as.character(sex))
TB_inc_df<-read.csv('KZN_GBD_TB_inc_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))
TB_mort_df<-read.csv('KZN_GBD_TB_mort_rate_calibration_df.csv')%>%
  mutate(TB_HIV_coinfection = as.character(TB_HIV_coinfection),
         sex = as.character(sex))

setwd(paste0(here(), '/results/graphs'))
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
  print(i)
}

#missing_sims_prog_df<-as.data.frame(unique(outputs_combined_df$sim_id))

#program stats

#by gender in 2017
stats_df_TB_inc_mort_by_g_2017<-outputs_combined_df%>%
  filter(year == 2017,
         program_id == 1)%>%
  select(c("sim_id",
           'Tb_inc_neg_female_Y_per_100k_female', 
           'Tb_inc_pos_female_Y_per_100k_female',
           "Tb_inc_neg_male_Y_per_100k_male",
           "Tb_inc_pos_male_Y_per_100k_male",
           'Tb_mort_neg_female_Y_per_100k_female', 
           'Tb_mort_pos_female_Y_per_100k_female',
           "Tb_mort_neg_male_Y_per_100k_male",
           "Tb_mort_pos_male_Y_per_100k_male"))

stats_df_TB_inc_mort_by_g_2017<-melt(stats_df_TB_inc_mort_by_g_2017, id = "sim_id")
stats_df_TB_inc_mort_by_g_2017<-stats_df_TB_inc_mort_by_g_2017%>%
  mutate(variable = str_replace(variable, "pos", ""),
         variable = str_replace(variable, "neg", ""))%>%
  group_by(variable, sim_id)%>% #sum over hiv pos and hiv neg
  summarise(total = sum(value))%>%
  group_by(variable)%>%
  summarise(mean = round(mean(total)),
            min = round(min(total)),
            max = round(max(total)))

#by gender in 2027

#program 1
stats_df_TB_inc_mort_by_g_2027_P1<-outputs_combined_df%>%
  filter(year == 2027,
         program_id == 1)%>%
  select(c("sim_id",
           'Tb_inc_neg_female_Y_per_100k_female', 
           'Tb_inc_pos_female_Y_per_100k_female',
           "Tb_inc_neg_male_Y_per_100k_male",
           "Tb_inc_pos_male_Y_per_100k_male",
           'Tb_mort_neg_female_Y_per_100k_female', 
           'Tb_mort_pos_female_Y_per_100k_female',
           "Tb_mort_neg_male_Y_per_100k_male",
           "Tb_mort_pos_male_Y_per_100k_male"))

stats_df_TB_inc_mort_by_g_2027_P1<-melt(stats_df_TB_inc_mort_by_g_2027_P1, id = "sim_id")
stats_df_TB_inc_mort_by_g_2027_P1<-stats_df_TB_inc_mort_by_g_2027_P1%>%
  mutate(variable = str_replace(variable, "pos", ""),
         variable = str_replace(variable, "neg", ""))%>%
  group_by(variable, sim_id)%>% #sum over hiv pos and hiv neg
  summarise(total = sum(value))%>%
  group_by(variable)%>%
  summarise(mean = round(mean(total)),
            min = round(min(total)),
            max = round(max(total)))

###Program 3
stats_df_TB_inc_mort_by_g_2027_P3<-outputs_combined_df%>%
  filter(year == 2027,
         program_id == 3)%>%
  select(c("sim_id",
           'Tb_inc_neg_female_Y_per_100k_female', 
           'Tb_inc_pos_female_Y_per_100k_female',
           "Tb_inc_neg_male_Y_per_100k_male",
           "Tb_inc_pos_male_Y_per_100k_male",
           'Tb_mort_neg_female_Y_per_100k_female', 
           'Tb_mort_pos_female_Y_per_100k_female',
           "Tb_mort_neg_male_Y_per_100k_male",
           "Tb_mort_pos_male_Y_per_100k_male"))

stats_df_TB_inc_mort_by_g_2027_P3<-melt(stats_df_TB_inc_mort_by_g_2027_P3, id = "sim_id")
stats_df_TB_inc_mort_by_g_2027_P3<-stats_df_TB_inc_mort_by_g_2027_P3%>%
  mutate(variable = str_replace(variable, "pos", ""),
         variable = str_replace(variable, "neg", ""))%>%
  group_by(variable, sim_id)%>% #sum over hiv pos and hiv neg
  summarise(total = sum(value))%>%
  group_by(variable)%>%
  summarise(mean = round(mean(total)),
            min = round(min(total)),
            max = round(max(total)))

#by gender and hiv status in 2017 (and hiv prev)
stats_df_by_g_hiv_status_2017<-outputs_combined_df%>%
  filter(year == 2017,
         program_id == 1)%>%
  select(c("sim_id",
           'Tb_inc_neg_female_Y_per_100k_female', 
           'Tb_inc_pos_female_Y_per_100k_female',
           "Tb_inc_neg_male_Y_per_100k_male",
           "Tb_inc_pos_male_Y_per_100k_male",
           'Tb_mort_neg_female_Y_per_100k_female', 
           'Tb_mort_pos_female_Y_per_100k_female',
           "Tb_mort_neg_male_Y_per_100k_male",
           "Tb_mort_pos_male_Y_per_100k_male",
           "hiv_prev_per_100k_female",
           "hiv_prev_per_100k_male"))

stats_df_by_g_hiv_status_2017<-melt(stats_df_by_g_hiv_status_2017, 
                                    id = "sim_id")
stats_df_by_g_hiv_status_2017<-stats_df_by_g_hiv_status_2017%>%
  group_by(variable)%>%
  summarise(mean = round(mean(value)),
            min = round(min(value)),
            max = round(max(value)))

##2027 stats
stats_df_by_program_2027<-outputs_combined_df%>%
  filter(year == 2027)%>%
  select(c('sim_id',
           "program_id",
           "Tb_inc_Y_per_100k_ppl",
           "Tb_mort_Y_per_100k_ppl"))

stats_df_by_program_2027<-melt(stats_df_by_program_2027, 
                                    id = c("sim_id",
                                    "program_id"))

stats_df_by_program_2027_2<-stats_df_by_program_2027%>%
  group_by(variable, program_id)%>%
  summarise(mean = round(mean(value)),
            min = round(min(value)),
            max = round(max(value)))

##percent decrease
stats_df_by_program_2027_3<-stats_df_by_program_2027%>%
  mutate(program_id = paste0("p", program_id))

stats_df_by_program_2027_3<-dcast(stats_df_by_program_2027_3, 
                                 sim_id+variable~program_id)
stats_df_by_program_2027_3<-stats_df_by_program_2027_3%>% 
  mutate(P2_v_P1 = (1- (p2/p1)),
         P3_v_P1 = (1- (p3/p1)),
         P3_v_P2 = (1- (p3/p2)))%>%
  select(-c('sim_id', 'p1', 'p2', 'p3'))

stats_df_by_program_2027_3<-melt(stats_df_by_program_2027_3,
                                 id = "variable")

colnames(stats_df_by_program_2027_3)<-c("variable",
                                        "comparison",
                                        "value")

stats_df_by_program_2027_3<-stats_df_by_program_2027_3%>%
  group_by(variable, comparison)%>%
  summarise(mean = round(mean(value)*100, 1),
            min = round(min(value)*100, 1),
            max = round(max(value)*100, 1))


##
stats_df_by_g_hiv_status_2027<-outputs_combined_df%>%
  filter(year == 2027)%>%
  select(c("program_id",
    "sim_id",
           'Tb_inc_neg_female_Y_per_100k_female', 
           'Tb_inc_pos_female_Y_per_100k_female',
           "Tb_inc_neg_male_Y_per_100k_male",
           "Tb_inc_pos_male_Y_per_100k_male",
           'Tb_mort_neg_female_Y_per_100k_female', 
           'Tb_mort_pos_female_Y_per_100k_female',
           "Tb_mort_neg_male_Y_per_100k_male",
           "Tb_mort_pos_male_Y_per_100k_male"))

stats_df_by_g_hiv_status_2027<-melt(stats_df_by_g_hiv_status_2027,
                                    id = c("program_id", 
                                           "sim_id"))

stats_df_by_g_hiv_status_2027_summarised_total<-stats_df_by_g_hiv_status_2027%>%
  group_by(program_id, variable)%>%
  summarise(mean = round(mean(value)),
            min = round(min(value)),
            max = round(max(value)))%>%
  mutate(ref = paste0(mean, " (range ", min, "-", max, ")"))

stats_df_by_g_hiv_status_2027_perc_reduction_summarised<-stats_df_by_g_hiv_status_2027
stats_df_by_g_hiv_status_2027_perc_reduction_summarised$program_id<-paste0('p', stats_df_by_g_hiv_status_2027_perc_reduction_summarised$program_id)
stats_df_by_g_hiv_status_2027_perc_reduction_summarised<-dcast(stats_df_by_g_hiv_status_2027_perc_reduction_summarised,
                            sim_id+variable ~ program_id)
stats_df_by_g_hiv_status_2027_perc_reduction_summarised<-stats_df_by_g_hiv_status_2027_perc_reduction_summarised%>%
  mutate(p3_p1_perc_reduction = 1-p3/p1)%>%
  group_by(variable)%>%
  summarise(mean = mean(p3_p1_perc_reduction),
            min = min(p3_p1_perc_reduction),
            max = max(p3_p1_perc_reduction))%>%
  mutate(text = paste0(round(mean*100, 1), "% (range ", round(min*100, 1), "% - ", round(max*100, 1), "%)"))

stats_df_by_g_2027_1<-stats_df_by_g_hiv_status_2027%>%
  mutate(variable = str_replace(variable, "pos", ""),
         variable = str_replace(variable, "neg", ""))%>%
  group_by(sim_id, program_id, variable)%>%
  summarise(value = round(sum(value)))%>%
  ungroup()

stats_df_by_g_2027_1<-stats_df_by_g_2027_1%>%
  mutate(program_id = paste0("p", program_id))

stats_df_by_g_2027_1<-dcast(stats_df_by_g_2027_1,
                                     sim_id+variable ~ program_id)

stats_df_by_g_2027_1<-stats_df_by_g_2027_1%>% 
  mutate(P2_v_P1 = (1- (p2/p1)),
         P3_v_P1 = (1- (p3/p1)),
         P3_v_P2 = (1- (p3/p2)))%>%
  select(-c('sim_id', 'p1', 'p2', 'p3'))

stats_df_by_g_2027_1<-melt(stats_df_by_g_2027_1,
                                 id = "variable")

colnames(stats_df_by_g_2027_1)<-c("variable",
                                        "comparison",
                                        "value")

stats_df_by_g_2027_1<-stats_df_by_g_2027_1%>%
  group_by(variable, comparison)%>%
  summarise(mean = round(mean(value)*100, 1),
            min = round(min(value)*100, 1),
            max = round(max(value)*100, 1))


###negative inc and mort
stats_df_by_hiv_status_2027_2<-stats_df_by_g_hiv_status_2027%>%
  filter(grepl('neg', variable))%>%
  group_by(program_id, variable)%>%
  summarise(mean = round(mean(value)),
            min = round(min(value)),
            max = round(max(value)))%>%
  filter(program_id != 2)

#####COLORS FOR PROGRAM GRAPHS#######
colors_for_graphs <- brewer.pal(n = 10, name = "Paired")
colors_for_program_graph_fill<-colors_for_graphs[c(1, 3, 5)]
colors_for_program_graph_line<-colors_for_graphs[c(2, 4, 6)]


#add transparecy
for (i in 1:length(colors_for_program_graph_fill)){
  rbg_temp<-as.vector(col2rgb(colors_for_program_graph_fill[i]))
  mycol <- rgb(rbg_temp[1], rbg_temp[2], rbg_temp[3], 
               max = 255, alpha = 185, 
               names = colors_for_program_graph_fill[i])
  colors_for_program_graph_fill[i]<-mycol
}

#####graphs with program1 differences by hiv status gender and metric######
outputs_combined_df_program_hiv_status_gender_metric<-melt(outputs_combined_df, 
                                  id = c("year", "sim_id", "program_id"))

outputs_combined_df_program_hiv_status_gender_metric<-outputs_combined_df_program_hiv_status_gender_metric%>%
  left_join(graph_ref, by = c('variable'))%>%
  filter(grepl('Tuberculosis', measure))

outputs_combined_df_program_hiv_status_gender_metric$title<-
  paste0(outputs_combined_df_program_hiv_status_gender_metric$measure, 
         outputs_combined_df_program_hiv_status_gender_metric$sex, 
         "s and ",
         outputs_combined_df_program_hiv_status_gender_metric$hs)


stats_program_hiv_status_gender_metric_2027<-outputs_combined_df_program_hiv_status_gender_metric%>%
  filter(year == 2027)%>%
  mutate(program_id = paste0('program_', program_id))%>%
  group_by(sim_id, program_id, title)%>% #summarise over hs
  summarise(value = sum(value))

cast_stats_program_hiv_status_gender_metric_2027<-dcast(stats_program_hiv_status_gender_metric_2027,
                                                  title+sim_id ~ program_id)

cast_stats_program_hiv_status_gender_metric_2027<-cast_stats_program_hiv_status_gender_metric_2027%>%
  mutate(program_1_2_diff = (program_1 - program_2)/program_1,
         program_2_3_diff = (program_2 - program_3)/program_2,
         program_1_3_diff = (program_1 - program_3)/program_1)

cast_stats_program_hiv_status_gender_metric_2027<-cast_stats_program_hiv_status_gender_metric_2027%>%
  select(-c("program_1", "program_2", "program_3"))

melt_stats_program_hiv_status_gender_metric_2027<-melt(cast_stats_program_hiv_status_gender_metric_2027,
                                                 id = c("title", "sim_id"))

stats_program_hiv_status_gender_metric_2027_2<-melt_stats_program_hiv_status_gender_metric_2027%>%
  group_by(title, variable)%>%
  summarise(max=max(value),
            min=min(value),
            mean=mean(value))

setwd(outdir_stats)
write.csv(stats_program_hiv_status_gender_metric_2027_2, 
          'stats_program_hiv_status_gender_metric_2027.csv',
          row.names = FALSE)


for(current_var_eval in unique(outputs_combined_df_program_hiv_status_gender_metric$title)){
  
  outputs_df_select<-outputs_combined_df_program_hiv_status_gender_metric%>%
    filter(title == current_var_eval)%>%
    mutate(group = paste0('Program ', program_id, "\n(p = ", program_id, ")"))
  
  program_eval_summarised<-outputs_df_select%>%
    group_by(year, group)%>%
    summarise(upper_rate = max(value),
              lower_rate = min(value),
              val_rate = mean(value))
  
  program_eval_summarised_connect<-program_eval_summarised%>%
    filter(year == 2018)
  
  program_eval_summarised_connect2<-program_eval_summarised%>%
    filter(year == 2017)
  
  last_yr_calib_est<-program_eval_summarised_connect2$val_rate
  
  program_eval_summarised_connect2<-rbind(program_eval_summarised_connect2,
                                          program_eval_summarised_connect2,
                                          program_eval_summarised_connect2)
  
  program_eval_summarised_connect2$group<-c('Program 1\n(p = 1)',
                                            'Program 2\n(p = 2)',
                                            'Program 3\n(p = 3)')
  
  program_eval_summarised$group<-if_else(program_eval_summarised$year <= 2017,
                                         "Calibration Period", 
                                         program_eval_summarised$group)
  
  program_eval_summarised<-rbind(program_eval_summarised,
                                 program_eval_summarised_connect2,
                                 program_eval_summarised_connect)
  
  
  if (grepl('inc', current_var_eval)){
    measure_temp<-"TB incidence rate"
    if(grepl('negative', current_var_eval)){
      max_min_graph_y_axis<-c(0, 800, 200)
    } else {
      max_min_graph_y_axis<-c(0, 2400, 400)
    }
  } else if (grepl('mort', current_var_eval)){
    measure_temp<-"TB mortality rate"
    if(grepl('negative', current_var_eval)){
      max_min_graph_y_axis<-c(0, 150, 20)
    } else{
      max_min_graph_y_axis<-c(0, 800, 150)
    }
  }
  
  program_graph_temp<-ggplot(program_eval_summarised)+
    geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
    geom_point(aes(x = year, y = val_rate, colour = group, shape = group), size = 2)+
    geom_point(aes(x = 2017, y = last_yr_calib_est), colour = "black", size = 2)+
    theme(text = element_text(size=18, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=22),
          plot.title = element_text(hjust = .5, size=20),
          legend.background = element_rect(colour = "lightgrey"),
          legend.box.background = element_rect(colour = "black"))+
    geom_vline(xintercept = 2017, linetype="dashed", 
               color = "darkgrey", size=1.5)+
    annotate("text", x=2012, y=max_min_graph_y_axis[2]*.95, label= "calibration period", size = 7)+
    annotate("text", x=2022, y=max_min_graph_y_axis[2]*.95, label= "intervention period", size = 7)+
    #ggtitle(current_var_eval)+
    scale_x_continuous(name = 'year Y', breaks=c(seq(from = 1990, to = 2028, by = 3)))+
    scale_color_manual(values=c("black", colors_for_program_graph_line))+
    scale_fill_manual(values=c('grey75', colors_for_program_graph_fill))+
    ylim(max_min_graph_y_axis[1:2])+
    ylab(measure_temp)+
    scale_shape_manual(values=c(16, 16, 8, 5))
  
  file_name<-paste0(str_replace_all(current_var_eval, " ", "_"),".png")
  
  setwd(outdir_graphs_prog)
  png(file_name, width = 800, height = 600)
  print(program_graph_temp)
  dev.off()
}


#####graphs with program2 differences by gender and metric (not HIV status)######
#for manuscript

outputs_combined_df_program_gender_metric<-melt(outputs_combined_df, 
                                  id = c("year", "sim_id", "program_id"))
outputs_combined_df_program_gender_metric<-outputs_combined_df_program_gender_metric%>%
  left_join(graph_ref, by = c('variable'))%>%
  filter(grepl('Tuberculosis', measure))%>%
  filter(!grepl('ppl', variable))


outputs_combined_df_program_gender_metric$title<-paste0(outputs_combined_df_program_gender_metric$measure, 
                                          outputs_combined_df_program_gender_metric$sex, "s")


program_gender_metric_summarised_2027<-outputs_combined_df_program_gender_metric%>%
  filter(year == 2027)%>%
  mutate(program_id = paste0('program_', program_id))%>%
  group_by(sim_id, program_id, title)%>% #summarise over hs
  summarise(value = sum(value))

cast_program_gender_metric_summarised_2027<-dcast(program_gender_metric_summarised_2027,
        title+sim_id ~ program_id)

cast_program_gender_metric_summarised_2027<-cast_program_gender_metric_summarised_2027%>%
  mutate(program_1_2_diff = (program_1 - program_2)/program_1,
         program_2_3_diff = (program_2 - program_3)/program_2,
         program_1_3_diff = (program_1 - program_3)/program_1)

cast_program_gender_metric_summarised_2027<-cast_program_gender_metric_summarised_2027%>%
  select(-c("program_1", "program_2", "program_3"))

melt_program_gender_metric_summarised_2027<-melt(cast_program_gender_metric_summarised_2027,
                                                 id = c("title", "sim_id"))

program_gender_metric_summarised_2027_2<-melt_program_gender_metric_summarised_2027%>%
  group_by(title, variable)%>%
  summarise(max=max(value),
            min=min(value),
            mean=mean(value))



setwd(outdir_stats)
write.csv(program_gender_metric_summarised_2027_2,
          'program_gender_metric_summarised_2027.csv',
          row.names = FALSE)

# gender_metric_prog_diff_stats_2027<-outputs_combined_df_program_gender_metric%>%
#   filter(year == 2027)%>%
#   mutate(program_id = paste0("program_", program_id))%>%
#   group_by(program_id, title)%>%
#   summarise(val = mean(value),
#             max = max(value),
#             min = min(value))
# 
# gender_metric_prog_diff_stats_2027_cast<-dcast(gender_metric_prog_diff_stats_2027,
#                                                title ~ program_id)
# 
# gender_metric_prog_diff_stats_2027_cast<-gender_metric_prog_diff_stats_2027_cast%>%
#   mutate(program_1_2_diff = 1 - program_2/program_1,
#          program_2_3_diff = 1 - program_3/program_2)
# 
# setwd(outdir_stats)
# write.csv(gender_metric_prog_diff_stats_2027_cast, 
#           'gender_metric_prog_diff_stats_2027.csv',
#           row.names = FALSE)


for(current_var_eval in unique(outputs_combined_df_program_gender_metric$title)){
  
  outputs_df_select<-outputs_combined_df_program_gender_metric%>%
    filter(title == current_var_eval)%>%
    mutate(group = paste0('Program ', program_id))%>%
    group_by(year, group, sim_id)%>%
    summarise(value = sum(value))
  
  program_eval_summarised<-outputs_df_select%>%
    group_by(year, group)%>%
    summarise(upper_rate = max(value),
              lower_rate = min(value),
              val_rate = mean(value))
  
  ##to connect 2017 ending at calibration to 2018 starting at evaluation
  program_eval_summarised_connect<-program_eval_summarised%>%
    filter(year == 2018)
  
  program_eval_summarised_connect2<-program_eval_summarised%>%
    filter(year == 2017)
  
  last_yr_calib_est<-program_eval_summarised_connect2$val_rate
  
  program_eval_summarised_connect2<-rbind(program_eval_summarised_connect2,
                                          program_eval_summarised_connect2,
                                          program_eval_summarised_connect2)
  program_eval_summarised_connect2$group<-c('Program 1',
                                            'Program 2',
                                            'Program 3')
  
  program_eval_summarised$group<-if_else(program_eval_summarised$year <= 2017,
                                         "Calibration Period", program_eval_summarised$group)
  
  program_eval_summarised<-rbind(program_eval_summarised,
                                 program_eval_summarised_connect2,
                                 program_eval_summarised_connect)
  
  if (grepl('inc', current_var_eval)){
    measure_temp<-"TB incidence rate"
    max_min_graph_y_axis<-c(0, 2500, 500)
  } else if (grepl('mort', current_var_eval)){
    measure_temp<-"TB mortality rate"
    max_min_graph_y_axis<-c(0, 700, 100)
  }
  
  yaxis_label<-gsub(" per", ", \nper", current_var_eval)
  yaxis_label<-paste0(yaxis_label, ", per year")
  
  title<-gsub(" per", ", per", current_var_eval)
  title<-paste0(title, ", per year")
  
  program_graph_temp<-ggplot(program_eval_summarised)+
    geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
    geom_point(aes(x = year, y = val_rate, colour = group, shape = group), 
               size = 2)+
    geom_point(aes(x = 2017, y = last_yr_calib_est), colour="black")+
    theme(text = element_text(size=18, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=18),
          plot.title = element_text(hjust = .5, size=20))+
    geom_vline(xintercept = 2017, linetype="dashed", 
               color = "darkgrey", size=1.5)+
    annotate("text", x=2012, y=max_min_graph_y_axis[2]*.95, label= "calibration period", size = 6.5)+
    annotate("text", x=2022, y=max_min_graph_y_axis[2]*.95, label= "intervention period", size = 6.5)+
    ggtitle(title)+
    scale_x_continuous(name = 'Year', breaks=c(seq(from = 1990, to = 2028, by = 4)))+
    scale_color_manual(values=c("black", colors_for_program_graph_line))+
    scale_fill_manual(values=c('grey75', colors_for_program_graph_fill))+
    ylim(max_min_graph_y_axis[1:2])+
    ylab(yaxis_label)+
    scale_shape_manual(values=c(16, 16, 8, 5))
  
  file_name<-paste0(str_replace_all(current_var_eval, " ", "_"),".png")
  
  setwd(outdir_graphs_prog2)
  png(file_name, width = 800, height = 600)
  print(program_graph_temp)
  dev.off()
}


######program3 mort prev and incidence differences over eval period######

metric_diff_stats_2027<-outputs_combined_df%>%
  filter(year == 2027)%>%
  mutate(program_id = paste0("program_", program_id))%>%
  select(c('program_id', "sim_id", 
           "tb_prev_Y_per_100K_ppl", 
           "Tb_inc_Y_per_100k_ppl",
           "Tb_mort_Y_per_100k_ppl"))

melt_metric_diff_stats_2027<-melt(metric_diff_stats_2027,
                                  id = c("program_id", "sim_id"))

cast_program_metric_summarised_2027<-dcast(melt_metric_diff_stats_2027,
                                                  sim_id + variable ~ program_id)

cast_program_metric_summarised_2027<-cast_program_metric_summarised_2027%>%
  mutate(program_1_2_diff = (program_1 - program_2)/program_1,
         program_2_3_diff = (program_2 - program_3)/program_2,
         program_1_3_diff = (program_1 - program_3)/program_1)

cast_program_metric_summarised_2027<-cast_program_metric_summarised_2027%>%
  mutate(measure = variable)%>%
  select(-c("program_1", "program_2", "program_3", "variable"))

melt_2_program_metric_summarised_2027<-melt(cast_program_metric_summarised_2027,
                                                 id = c("measure", "sim_id"))

program_metric_summarised_2027_2<-melt_2_program_metric_summarised_2027%>%
  group_by(variable, measure)%>%
  summarise(max=max(value),
            min=min(value),
            mean=mean(value))

setwd(outdir_stats)
write.csv(program_metric_summarised_2027_2,
          'program_metric_summarised_2027.csv',
          row.names = FALSE)

####program 3 for appendix (with p = ...)####
for(measure in c('mort', 'prev', 'inc')){
  
  measure2<-paste0('tb_', measure, '_y_per_100k_ppl')
  measure_y_axis<-if_else(grepl('mort', measure), 'TB mortality rate,\nper 100,000 individuals, per year',
                          if_else(grepl('prev', measure), 
                                  'TB prevalence rate,\nper 100,000 individuals, per year',
                                  'TB incidence rate,\nper 100,000 individuals, per year'))
  measure_title<-if_else(grepl('mort', measure), 'TB mortality rate, per 100,000 individuals, per year',
                         if_else(grepl('prev', measure), 
                                 'TB prevalence rate, per 100,000 individuals, per year',
                                 'TB incidence rate, per 100,000 individuals, per year'))
  
  
  outputs_current_measure_df<-outputs_combined_df
  
  names(outputs_current_measure_df)<-tolower(names(outputs_current_measure_df))
  
  outputs_current_measure_df<-outputs_current_measure_df%>%
    select(c('year', 'program_id', measure2))%>%
    rename(value = measure2)%>%
    group_by(year, program_id)%>%
    summarise(lower_rate = min(value),
              upper_rate = max(value),
              val_rate = mean(value))%>%
    mutate(Program = paste0('Program ', program_id, "\n(p = ", program_id, ")"))%>%
    filter(year >= 2017)%>%
    ungroup()
  
  
  ###so that all programs start at the same place before eval period
  program_eval_summarised_connect2<-outputs_current_measure_df%>%
    filter(year == 2017)
  
  last_yr_calib_est<-program_eval_summarised_connect2$val_rate
  
  program_eval_summarised_connect2<-rbind(program_eval_summarised_connect2,
                                          program_eval_summarised_connect2,
                                          program_eval_summarised_connect2)
  program_eval_summarised_connect2$Program<-c('Program 1\n(p = 1)',
                                              'Program 2\n(p = 2)',
                                              'Program 3\n(p = 3)')
  
  outputs_current_measure_df<-rbind(outputs_current_measure_df,
                                    program_eval_summarised_connect2)
  
  prog1_val<-outputs_current_measure_df%>%
    filter(year == 2027)%>%
    filter(program_id == 1)%>%
    select('val_rate')
  
  prog2_val<-outputs_current_measure_df%>%
    filter(year == 2027)%>%
    filter(program_id == 2)%>%
    select('val_rate')
  
  prog3_val<-outputs_current_measure_df%>%
    filter(year == 2027)%>%
    filter(program_id == 3)%>%
    select('val_rate')
  
  prog1_v_2<-program_metric_summarised_2027_2%>%
    mutate(measure = tolower(measure))%>%
    filter(measure == measure2,
           variable == "program_1_2_diff")
  
  prog2_v_3<-program_metric_summarised_2027_2%>%
    mutate(measure = tolower(measure))%>%
    filter(measure == measure2,
           variable == "program_2_3_diff")
  
  prog1_v_2<-prog1_v_2$mean[[1]]
  prog1_v_2<-paste(round(100*prog1_v_2, 1), "%", sep="")
  
  prog2_v_3<-prog2_v_3$mean[[1]]
  prog2_v_3<-paste(round(100*prog2_v_3, 1), "%", sep="")
  
  program_graph_temp<-ggplot(outputs_current_measure_df)+
    geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = Program))+
    geom_point(aes(x = year, y = val_rate, colour = Program, shape = Program), size = 3)+
    theme(text = element_text(size=18, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=18),
          plot.title = element_text(hjust = .5, size=20))+
    scale_x_continuous(name = 'year Y',
                       breaks=c(seq(from = 2018, to = 2027, by = 3)),
                       expand = expansion(mult=c(0,.1)))+
    scale_color_manual(values=c(colors_for_program_graph_line))+
    scale_fill_manual(values=c(colors_for_program_graph_fill))+
    ylab(measure_y_axis)+
    geom_segment(aes(x = 2027.5, y = prog1_val[[1]], xend = 2027.2, 
                     yend = prog1_val[[1]]),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    geom_segment(aes(x = 2027.5, y = prog1_val[[1]], xend = 2027.5, 
                     yend = prog2_val[[1]]*1.02))+
    geom_segment(aes(x = 2027.5, y = prog2_val[[1]]*1.02, xend = 2027.2, 
                     yend = prog2_val[[1]]*1.02),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    geom_segment(aes(x = 2027.5, y = prog2_val[[1]]*.98, xend = 2027.2, 
                     yend = prog2_val[[1]]*.98),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    geom_segment(aes(x = 2027.5, y = prog2_val[[1]]*.98, xend = 2027.5, 
                     yend = prog3_val[[1]]))+
    geom_segment(aes(x = 2027.5, y = prog3_val[[1]], xend = 2027.2, 
                     yend = prog3_val[[1]]),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    annotate("text", x=2028.8, y=((prog1_val[[1]]+prog2_val[[1]])/2)*1.1, 
             label= paste0('Program 2\n', prog1_v_2, ' reduction\nvs Program 1\nin 2027'), 
             #label= paste0('Community-ART\n', prog1_v_2, ' reduction\nin 2027'), 
             size = 6.5)+
    annotate("text", x=2028.8, y=(prog2_val[[1]]*.9+prog3_val[[1]])/2, 
             label= paste0('Program 3\n', prog2_v_3, ' reduction\nvs Program 2\nin 2027'), 
             #label= paste0('Community-IPT\n', prog2_v_3, ' reduction\nin 2027'), 
             size = 6.5)+
    scale_shape_manual(values=c(16, 8, 5))
  
  
  
  
  
  file_name<-paste0(measure,".png")
  
  setwd(outdir_graphs_prog3)
  png(file_name, width = 800, height = 600)
  print(program_graph_temp)
  dev.off()
}

####program 3 for manuscript####
for(measure in c('mort', 'prev', 'inc')){
  
  measure2<-paste0('tb_', measure, '_y_per_100k_ppl')
  measure_y_axis<-if_else(grepl('mort', measure), 'TB mortality rate,\nper 100,000 individuals, per year',
                          if_else(grepl('prev', measure), 
                                  'TB prevalence rate,\nper 100,000 individuals, per year',
                                  'TB incidence rate,\nper 100,000 individuals, per year'))
  measure_title<-if_else(grepl('mort', measure), 'TB mortality rate for each program over the intervention period',
                         if_else(grepl('prev', measure), 
                                 'TB prevalence rate for each program over the intervention period',
                                 'TB incidence rate for each program over the intervention period'))
  
  
  outputs_current_measure_df<-outputs_combined_df
  
  names(outputs_current_measure_df)<-tolower(names(outputs_current_measure_df))
  
  outputs_current_measure_df<-outputs_current_measure_df%>%
    select(c('year', 'program_id', measure2))%>%
    rename(value = measure2)%>%
    group_by(year, program_id)%>%
    summarise(lower_rate = min(value),
              upper_rate = max(value),
              val_rate = mean(value))%>%
    mutate(Program = paste0('Program ', program_id))%>%
    filter(year >= 2017)%>%
    ungroup()
  
  
  ###so that all programs start at the same place before eval period
  program_eval_summarised_connect2<-outputs_current_measure_df%>%
    filter(year == 2017)
  
  last_yr_calib_est<-program_eval_summarised_connect2$val_rate
  
  program_eval_summarised_connect2<-rbind(program_eval_summarised_connect2,
                                          program_eval_summarised_connect2,
                                          program_eval_summarised_connect2)
  program_eval_summarised_connect2$Program<-c('Program 1',
                                              'Program 2',
                                              'Program 3')
  
  outputs_current_measure_df<-rbind(outputs_current_measure_df,
                                    program_eval_summarised_connect2)
  
  prog1_val<-outputs_current_measure_df%>%
    filter(year == 2027)%>%
    filter(program_id == 1)%>%
    select('val_rate')
  
  prog2_val<-outputs_current_measure_df%>%
    filter(year == 2027)%>%
    filter(program_id == 2)%>%
    select('val_rate')
  
  prog3_val<-outputs_current_measure_df%>%
    filter(year == 2027)%>%
    filter(program_id == 3)%>%
    select('val_rate')
  
  prog1_v_2<-program_metric_summarised_2027_2%>%
    mutate(measure = tolower(measure))%>%
    filter(measure == measure2,
           variable == "program_1_2_diff")
  
  prog2_v_3<-program_metric_summarised_2027_2%>%
    mutate(measure = tolower(measure))%>%
    filter(measure == measure2,
           variable == "program_2_3_diff")
  
  prog1_v_2<-prog1_v_2$mean[[1]]
  prog1_v_2<-paste(round(100*prog1_v_2, 1), "%", sep="")
  
  prog2_v_3<-prog2_v_3$mean[[1]]
  prog2_v_3<-paste(round(100*prog2_v_3, 1), "%", sep="")
  
  program_graph_temp<-ggplot(outputs_current_measure_df)+
    geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = Program))+
    geom_point(aes(x = year, y = val_rate, colour = Program, shape = Program), size = 3)+
    theme(text = element_text(size=18, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=18),
          plot.title = element_text(hjust = .5, size=20))+
    scale_x_continuous(name = 'Year',
                       breaks=c(seq(from = 2018, to = 2027, by = 3)),
                       expand = expansion(mult=c(0,.1)))+
    scale_color_manual(values=c(colors_for_program_graph_line))+
    scale_fill_manual(values=c(colors_for_program_graph_fill))+
    ylab(measure_y_axis)+
    geom_segment(aes(x = 2027.5, y = prog1_val[[1]], xend = 2027.2, 
                     yend = prog1_val[[1]]),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    geom_segment(aes(x = 2027.5, y = prog1_val[[1]], xend = 2027.5, 
                     yend = prog2_val[[1]]*1.02))+
    geom_segment(aes(x = 2027.5, y = prog2_val[[1]]*1.02, xend = 2027.2, 
                     yend = prog2_val[[1]]*1.02),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    geom_segment(aes(x = 2027.5, y = prog2_val[[1]]*.98, xend = 2027.2, 
                     yend = prog2_val[[1]]*.98),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    geom_segment(aes(x = 2027.5, y = prog2_val[[1]]*.98, xend = 2027.5, 
                     yend = prog3_val[[1]]))+
    geom_segment(aes(x = 2027.5, y = prog3_val[[1]], xend = 2027.2, 
                     yend = prog3_val[[1]]),
                 arrow = arrow(length = unit(0.3, "cm"), type = 'closed'))+
    annotate("text", x=2028.8, y=((prog1_val[[1]]+prog2_val[[1]])/2)*1.1, 
             label= paste0('Program 2\n', prog1_v_2, ' reduction\nvs Program 1\nin 2027'), 
             #label= paste0('Community-ART\n', prog1_v_2, ' reduction\nin 2027'), 
             size = 6.5)+
    annotate("text", x=2028.8, y=(prog2_val[[1]]*.9+prog3_val[[1]])/2, 
             label= paste0('Program 3\n', prog2_v_3, ' reduction\nvs Program 2\nin 2027'), 
             #label= paste0('Community-IPT\n', prog2_v_3, ' reduction\nin 2027'), 
             size = 6.5)+
    scale_shape_manual(values=c(16, 8, 5))+
    ggtitle(measure_title)
  
  
  
  
  
  file_name<-paste0(measure,".png")
  
  setwd(outdir_graphs_prog4)
  png(file_name, width = 800, height = 600)
  print(program_graph_temp)
  dev.off()
}

###COLORS FOR CALIBRATION GRAPHS####

#display.brewer.pal(n = 10, name = "Paired")
colors_model_projections_fill<-colors_for_graphs[c(7, 9, 3)] #orange, purple, green
colors_model_projections_fill_transparent<-colors_for_graphs[c(7, 9, 3)] #orange, purple, green
for (i in 1:length(colors_model_projections_fill)){
  rbg_temp<-as.vector(col2rgb(colors_model_projections_fill[i]))
  mycol <- rgb(rbg_temp[1], rbg_temp[2], rbg_temp[3], 
               max = 255, alpha = 50, 
               names = colors_model_projections_fill[i])
  colors_model_projections_fill_transparent[i]<-mycol
}

colors_model_projections_line<-c("#d95f02", "#7570b3", "#1b9e77")

colors_GBD_estimates_line<-colors_for_graphs[c(2)]

#####calibration graphs1 by HIV status, metric and gender#######

outputs_combined_df_calibration<-outputs_combined_df%>%
  filter(program_id == 1)%>%
  select(-c('program_id'))
outputs_combined_df_calibration<-melt(outputs_combined_df_calibration, 
                                      id = c("year", "sim_id"))

outputs_combined_df_calibration<-outputs_combined_df_calibration%>%
  filter(!grepl('ppl', variable))

outputs_combined_df_calibration<-outputs_combined_df_calibration%>%
  left_join(graph_ref, by = c('variable'))

outputs_combined_df_calibration$title<-paste0(outputs_combined_df_calibration$measure, 
                                              outputs_combined_df_calibration$sex, 
                                              if_else(grepl("Tuberculosis", 
                                                            outputs_combined_df_calibration$measure),
                                                      paste0("s and ",
                                                      outputs_combined_df_calibration$hs), "s"))


                                          
for(current_var_eval in unique(outputs_combined_df_calibration$title)){
  
  outputs_df_select<-outputs_combined_df_calibration%>%
    filter(title == current_var_eval)
  
  outputs_df_summarised<-outputs_df_select%>%
    group_by(year)%>%
    summarise(upper_rate = max(value),
              lower_rate = min(value),
              val_rate = mean(value))%>%
    filter(year <= 2017)
  outputs_df_summarised$group = "Model Projections"
  
  sex_temp<-unique(outputs_df_select$sex)
  TB_HIV_coinfection_temp<-unique(outputs_df_select$TB_HIV_coinfection)
  
  if (grepl('inc', current_var_eval)){
    measure_temp<-"TB incidence rate"
    GBD_df<-TB_inc_df%>%
      filter(TB_HIV_coinfection == TB_HIV_coinfection_temp,
             sex == sex_temp)%>%
      filter(year == 2017|year == 2005)
    if(grepl('negative', current_var_eval)){
      max_min_graph_y_axis<-c(0, 1200, 200)
    } else {
      max_min_graph_y_axis<-c(0, 2400, 400)
    }
  } else if (grepl('mort', current_var_eval)){
    measure_temp<-"TB mortality rate"
    GBD_df<-TB_mort_df%>%
      filter(TB_HIV_coinfection == TB_HIV_coinfection_temp,
             sex == sex_temp)%>%
      filter(year == 2017|year == 2005)
    if(grepl('negative', current_var_eval)){
      max_min_graph_y_axis<-c(0, 260, 40)
    } else{
      max_min_graph_y_axis<-c(0, 800, 200)
    }
  } else {
    GBD_df<-HIV_prev_df%>%
      filter(sex == sex_temp)%>%
      filter(year == 2017|year == 2005)
    measure_temp<-"HIV prevalence rate"
    max_min_graph_y_axis<-c(0, 50000, 10000)
  }
  
  calib_graph_temp<-ggplot(outputs_df_summarised)+
    geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year), 
                fill = colors_model_projections_fill[1])+
    geom_line(aes(x = year, y = val_rate), size = .7,
              color = colors_model_projections_line[1])+
    theme(text = element_text(size=20, family="Times New Roman"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="top", legend.title = element_blank(), legend.text=element_text(size=18),
          plot.title = element_text(hjust = .5, size=20))+
    ggtitle(current_var_eval)+
    scale_x_continuous(name = 'year Y', breaks=c(seq(from = 1990, to = 2017, by = 3)))+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+
    ylim(max_min_graph_y_axis[1:2])+
    ylab(measure_temp)+
    geom_segment(x = 2005, y = unlist(GBD_df%>%
                       filter(year == 2005)%>%
                       select(c("lower_rate"))), 
                     xend = 2005, yend = min(unlist(GBD_df%>%
                       filter(year == 2005)%>%
                       select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                     col = colors_GBD_estimates_line, size = 1)+
    geom_segment(x = 2017, y = unlist(GBD_df%>%
                                        filter(year == 2017)%>%
                                        select(c("lower_rate"))), 
                 xend = 2017, yend = min(unlist(GBD_df%>%
                                                  filter(year == 2017)%>%
                                                  select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                 col = colors_GBD_estimates_line, size = 1)+
    geom_point(data = GBD_df, aes(x = year, y = val_rate), 
               col = colors_GBD_estimates_line, size = 6)
  
  file_name<-paste0("calibration_",
                    str_replace_all(current_var_eval, " ", "_"),
                    ".png")
                    
  
  setwd(outdir_graphs_calib)
  png(file_name, width = 800, height = 600)
  print(calib_graph_temp)
  dev.off()
}
  


#display.brewer.pal(n = 10, name = "Paired")
colors_for_graphs <- brewer.pal(n = 10, name = "Paired")

colors_for_calib_graph_fill_tb<-colors_for_graphs[c(1, 3, 5, 7)]


for (i in 1:length(colors_for_calib_graph_fill_tb)){
  rbg_temp<-as.vector(col2rgb(colors_for_calib_graph_fill_tb[i]))
  mycol <- rgb(rbg_temp[1], rbg_temp[2], rbg_temp[3],
               max = 255, alpha = 125,
               names = colors_for_calib_graph_fill_tb[i])
  colors_for_calib_graph_fill_tb[i]<-mycol
}

###COLORS FOR CALIBRATION GRAPHS both HIV NEG AND POSITIVE on the same graph#####
colors_for_calib_graph_line_tb<-colors_for_graphs[c(2, 4, 6, 8)]

colors_for_calib_graph_fill_hiv<-colors_for_graphs[c(3, 7)]

for (i in 1:length(colors_for_calib_graph_fill_hiv)){
  rbg_temp<-as.vector(col2rgb(colors_for_calib_graph_fill_hiv[i]))
  mycol <- rgb(rbg_temp[1], rbg_temp[2], rbg_temp[3], 
               max = 255, alpha = 125, 
               names = colors_for_calib_graph_fill_hiv[i])
  colors_for_calib_graph_fill_hiv[i]<-mycol
}

colors_for_calib_graph_line_hiv<-colors_for_graphs[c(4, 8)]


###calibration graphs2 by gender and metric (includes both HIV status) ####
outputs_combined_df_calibration<-outputs_combined_df%>%
  filter(program_id == 1)%>%
  select(-c('program_id'))
outputs_combined_df_calibration<-melt(outputs_combined_df_calibration, id = c("year", "sim_id"))
outputs_combined_df_calibration<-outputs_combined_df_calibration%>%
  left_join(graph_ref, by = c('variable'))%>%
  filter(!grepl('ppl', variable))

outputs_combined_df_calibration$title<-paste0(outputs_combined_df_calibration$measure, 
                                              outputs_combined_df_calibration$sex, "s")

outputs_combined_df_calibration$yaxis<-paste0(outputs_combined_df_calibration$measure, 
                                              outputs_combined_df_calibration$sex, "s", 
                                              " per year")

for(current_measure_eval in unique(outputs_combined_df_calibration$title)){
  
  outputs_df_select<-outputs_combined_df_calibration%>%
    filter(title == current_measure_eval)
  
  yaxis_label<-gsub(" per", ", \nper", current_measure_eval)
  yaxis_label<-paste0(yaxis_label, ", per year")
  
  
  if(!grepl('hiv', as.character(current_measure_eval), ignore.case = TRUE)){
    outputs_df_summarised<-outputs_df_select%>%
      group_by(year, hs)%>%
      summarise(upper_rate = max(value),
                lower_rate = min(value),
                val_rate = mean(value))%>%
      mutate(group = paste0('Model Projections: ', hs))%>%
      select(-c('hs'))%>%
      filter(year <= 2017)
  } else {
    outputs_df_summarised<-outputs_df_select%>%
      group_by(year, hs)%>%
      summarise(upper_rate = max(value),
                lower_rate = min(value),
                val_rate = mean(value))%>%
      select(-c('hs'))%>%
      filter(year <= 2017)
    outputs_df_summarised$group = "Model Projections"
  }
  
  
  sex_temp<-unique(outputs_df_select$sex)
  
  if (grepl('inc', current_measure_eval)){
    GBD_df<-TB_inc_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2017)
    measure_temp<-"TB incidence rate"
    max_min_graph_y_axis<-c(0, 2400, 400)
  } else if (grepl('mort', current_measure_eval)){
    GBD_df<-TB_mort_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2017)
    measure_temp<-"TB mortality rate"
    max_min_graph_y_axis<-c(0, 1000, 100)
  } else {
    GBD_df<-HIV_prev_df%>%
      filter(sex == sex_temp)%>%
      filter(year <= 2017)
    measure_temp<-"HIV prevalence rate"
    max_min_graph_y_axis<-c(0, 50000, 10000)
  }
  
  if(!grepl('hiv', as.character(current_measure_eval), ignore.case = TRUE)) {
    GBD_df$group = paste0("GBD Projections: ", 
                          if_else(GBD_df$TB_HIV_coinfection == "no", 
                                  "HIV-", "HIV+"))
  } else {
    GBD_df$group <- "GBD Projections"
  }
  
  calibration_graph_df<-rbind(data.frame(GBD_df%>%select(year, val_rate, upper_rate, lower_rate, group)),
                              data.frame(outputs_df_summarised))
  
  if(!grepl('hiv', as.character(current_measure_eval), ignore.case = TRUE)){
    calib_graph_temp<-ggplot(calibration_graph_df)+
      geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
      geom_point(data = calibration_graph_df%>%filter(grepl("Model Projections", group),
                                                     grepl("-", group)), 
                aes(x = year, y = val_rate), size = 2, color = "darkred")+
      geom_point(data = calibration_graph_df%>%filter(grepl("Model Projections", group),
                                                     !grepl("-", group)), 
                aes(x = year, y = val_rate), size = 2, color = "darkorange")+
      theme(text = element_text(size=22, family="Times New Roman"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="top", legend.title = element_blank(), legend.text=element_text(size=22),
            plot.title = element_text(hjust = .5, size=20),
            legend.background = element_rect(colour = "lightgrey"),
            legend.box.background = element_rect(colour = "black"))+
      #ggtitle(current_measure_eval)+
      scale_x_continuous(name = 'year Y', breaks=c(seq(from = 1990, to = 2018, by = 3)))+
      
      
      #scale_y_continuous(name = yaxis_label, breaks = c(seq(from = max_min_graph_y_axis[1],
      #                                                       to = max_min_graph_y_axis[2],
      #                                                       by = max_min_graph_y_axis[3])))+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))+
      #scale_color_manual(values=colors_for_calib_graph_line_tb[1,3])+
      scale_fill_manual(values=colors_for_calib_graph_fill_tb)+
      ylim(max_min_graph_y_axis[1:2])+
      ylab(yaxis_label)+
      geom_segment(x = 2017-0.2, y = unlist(GBD_df%>%
                                          filter(year == 2017, grepl("-", group))%>%
                                          select(c("lower_rate"))), 
                   xend = 2017-0.2, yend = min(unlist(GBD_df%>%
                                                    filter(year == 2017, 
                                                           grepl("-", group))%>%
                                                    select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                   col = "darkblue", size = 1.2, 
                   arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
      geom_segment(x = 2005-0.2, y = unlist(GBD_df%>%
                                          filter(year == 2005, grepl("-", group))%>%
                                          select(c("lower_rate"))), 
                   xend = 2005-0.2, yend = min(unlist(GBD_df%>%
                                                    filter(year == 2005, 
                                                           grepl("-", group))%>%
                                                    select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                   col = "darkblue", size = 1.2, 
                   arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
      geom_segment(x = 2017+0.2, y = unlist(GBD_df%>%
                                          filter(year == 2017, !grepl("-", group))%>%
                                          select(c("lower_rate"))), 
                   xend = 2017+0.2, yend = min(unlist(GBD_df%>%
                                                    filter(year == 2017, 
                                                           !grepl("-", group))%>%
                                                    select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                   col = "darkgreen", size = 1.2, 
                   arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
      geom_segment(x = 2005+0.2, y = unlist(GBD_df%>%
                                          filter(year == 2005, !grepl("-", group))%>%
                                          select(c("lower_rate"))), 
                   xend = 2005+0.2, yend = min(unlist(GBD_df%>%
                                                    filter(year == 2005, 
                                                           !grepl("-", group))%>%
                                                    select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                   col = "darkgreen", size = 1.2, 
                   arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))
    
  } else {
    calib_graph_temp<-ggplot(calibration_graph_df)+
      geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, x = year, fill = group))+
      geom_point(data = calibration_graph_df%>%filter(group == "Model Projections"), 
                aes(x = year, y = val_rate), size = 2, color = "darkorange")+
      theme(text = element_text(size=22, family="Times New Roman"), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position="top", legend.title = element_blank(), legend.text=element_text(size=18),
            plot.title = element_text(hjust = .5, size=20),
            legend.background = element_rect(colour = "lightgrey"),
            legend.box.background = element_rect(colour = "black"))+
      #ggtitle(current_measure_eval)+
      scale_x_continuous(name = 'year Y', breaks=c(seq(from = 1990, to = 2017, by = 3)))+
      guides(color=guide_legend(nrow=1,byrow=TRUE))+
      scale_color_manual(values=colors_for_calib_graph_line_hiv)+
      scale_fill_manual(values=colors_for_calib_graph_fill_hiv)+
      ylim(max_min_graph_y_axis[1:2])+
      ylab(yaxis_label)+
      geom_segment(x = 2005, y = unlist(GBD_df%>%
                                          filter(year == 2005)%>%
                                          select(c("lower_rate"))), 
                   xend = 2005, yend = min(unlist(GBD_df%>%
                                                    filter(year == 2005)%>%
                                                    select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                   col = "darkgreen", size = 1.2, 
                   arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))+
      geom_segment(x = 2017, y = unlist(GBD_df%>%
                                          filter(year == 2017)%>%
                                          select(c("lower_rate"))), 
                   xend = 2017, yend = min(unlist(GBD_df%>%
                                                    filter(year == 2017)%>%
                                                    select(c("upper_rate"))), max_min_graph_y_axis[2]), 
                   col = "darkgreen", size = 1.2, 
                   arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")))
  }
  
  file_name<-paste0("calibration2_",
                    str_replace_all(current_measure_eval, " ", "_"),
                    ".png")
  
  
  setwd(outdir_graphs_calib2)
  png(file_name, width = 800, height = 500)
  print(calib_graph_temp)
  dev.off()
}




