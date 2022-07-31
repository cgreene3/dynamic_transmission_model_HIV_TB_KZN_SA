#Estimating background mortality rate from GBD
#Uses GBD 2019
#(Total Mort - HIV Mort - TB mort)/pop estimate = base mort
#KZN population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa for here to work
indir <- paste0(here(),'/param_files/calculated_param_gen/input_data/GBD')
outdir <- paste0(here(),'/param_files/input_parameters')
graph_outdir<-paste0(here(),'/param_files/dynamic_param_graphs')

setwd(indir)
#read in pop estimates
pop_df<-read.csv('pop_estimates_all_ages.csv')

#read in all cause mort num and summarise
all_cause_mort_df<-read.csv('all_cause_mort_num_rate_pulled_2_22_22.csv')%>%
  filter(metric == 'Number')%>%
  group_by(year, sex)%>%
  summarise(all_cause_val = sum(val))

disease_mort_df<-read.csv('disease_mort_num_pulled_2_22_22.csv')%>%
  group_by(year, sex)%>%
  summarise(disease_val = sum(val))

base_mort_df<-pop_df%>%
  left_join(all_cause_mort_df, by = c('year', 'sex'))%>%
  left_join(disease_mort_df, by = c('year', 'sex'))%>%
  mutate(val = ((all_cause_val-disease_val)/expected_total_pop))%>%
  select('year', 'sex', 'val')

base_mort_df$max = base_mort_df$val*1.25
base_mort_df$min = base_mort_df$val*.75

setwd(outdir)
write.csv(base_mort_df, 'base_mort_df.csv', row.names = FALSE)

#make plots for appendix

base_mort_df_graph_male<-base_mort_df%>%
  filter(sex == 'Male')%>%
  select(-c('sex'))

base_mort_df_graph_male<-melt(base_mort_df_graph_male, id = c('year'))
  
baseline_mort_plot_male<-ggplot()+
  geom_ribbon(data = base_mort_df%>%filter(sex == 'Male'),
              aes(ymin = min, ymax = max, x = year), fill = "darkseagreen1")+
  #geom_line(data = base_mort_df%>%filter(sex == 'Male'), aes(x = year, y = val), color = 'darkgreen', size = .5)+
  geom_point(data = base_mort_df_graph_male,
             aes(x=year, y=value, group = variable, shape = variable), size = .8, color = 'darkgreen')+
  labs(title=(bquote(atop(mu[{1}]^{BASELINE}~(tau)~~Baseline~mortality, "           rate calibration ranges, Males"))))+
  scale_x_continuous(name = 'Year Y', breaks=seq(from = 1990, to = 2017, by = 5))+
  scale_y_continuous(name = 'baseline mortality rate', limits = c(0, .012), breaks=(seq(0, .012, .003)))+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())+
  scale_shape_manual(values=c(8, 17, 6))

base_mort_df_graph_female<-base_mort_df%>%
  filter(sex == 'Female')%>%
  select(-c('sex'))

base_mort_df_graph_female<-melt(base_mort_df_graph_female, id = c('year'))

baseline_mort_plot_female<-ggplot()+
  geom_ribbon(data = base_mort_df%>%filter(sex == 'Female'),
              aes(ymin = min, ymax = max, x = year), fill = "plum2")+
  geom_point(data = base_mort_df_graph_female,
             aes(x=year, y=value, group = variable, shape = variable), size = .8, color = 'darkorchid4')+
  labs(title=(bquote(atop(mu[{2}]^{BASELINE}~(tau)~~Baseline~mortality, "           rate calibration ranges, Females"))))+
  scale_x_continuous(name = 'Year Y', breaks=seq(from = 1990, to = 2017, by = 5))+
  scale_y_continuous(name = 'baseline mortality rate', limits = c(0, .012), breaks=(seq(0, .012, .003)))+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_blank())+
  scale_shape_manual(values=c(8, 17, 6))

setwd(graph_outdir)

png("baseline_mort_plot_male.png")
print(baseline_mort_plot_male)
dev.off()

png("baseline_mort_plot_female.png")
print(baseline_mort_plot_female)
dev.off()
