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
pop_df<-read.csv('pop_estimates_15_59.csv')

#read in all cause mort num and summarise
all_cause_mort_df<-read.csv('all_cause_mort_num_rate.csv')%>%
  filter(metric == 'Number')%>%
  group_by(year, sex)%>%
  summarise(all_cause_val = sum(val))

disease_mort_df<-read.csv('disease_mort_num.csv')%>%
  group_by(year, sex)%>%
  summarise(disease_val = sum(val))

base_mort_df<-pop_df%>%
  left_join(all_cause_mort_df, by = c('year', 'sex'))%>%
  left_join(disease_mort_df, by = c('year', 'sex'))%>%
  mutate(val = ((all_cause_val-disease_val)/expected_total_pop))%>%
  select('year', 'sex', 'val')

base_mort_df<-base_mort_df%>%
  filter(year <= 2017)

base_mort_df_2018_mort_vals_female<-base_mort_df%>%filter(year == 2017, sex == "Female")
base_mort_df_2018_mort_vals_male<-base_mort_df%>%filter(year == 2017, sex == "Male")


base_mort_df_2018_2028<-data.frame(year = rep(2018:2028, times = 2),
                                   sex = rep(c("Female", "Male"),
                                             each = length(2018:2028)),
                                   val = rep(c(base_mort_df_2018_mort_vals_female$val,
                                               base_mort_df_2018_mort_vals_male$val),
                                             each = length(2018:2028)))

base_mort_df<-rbind(base_mort_df, base_mort_df_2018_2028)
base_mort_df$max = base_mort_df$val*1.25
base_mort_df$min = base_mort_df$val*.75

setwd(outdir)
write.csv(base_mort_df, 'base_mort_df.csv', row.names = FALSE)

#make plots for appendix

base_mort_df_graph_male<-base_mort_df%>%
  filter(sex == 'Male')%>%
  select(-c('sex'))

#without max/min
baseline_mort_plot_male<-ggplot(data=base_mort_df_graph_male, 
       aes(x=year, y=val, group=1)) +
  geom_line(color="grey", size = 1)+
  geom_point(color="blue", size = 2)+
  geom_vline(xintercept = 2017, linetype="dashed", 
             color = "darkgrey", size=1)+
  annotate("text", x=2010, y=.011, label= "calibration period", size = 5.5)+
  annotate("text", x=2024, y=.011, label= "evaluation period", size = 5.5)+
  #labs(title=(bquote(atop("Baseline mortality rate values, Males"))))+
  scale_x_continuous(name = 'Year Y', breaks=seq(from = 1990, to = 2027, by = 4))+
  scale_y_continuous(name = bquote(atop(mu[{1}]^{VAL}~(Y))), limits = c(0, .012), breaks=(seq(0, .012, .003)))+
  theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))
# with max/min
# base_mort_df_graph_male<-melt(base_mort_df_graph_male, id = c('year'))
#   
# baseline_mort_plot_male<-ggplot()+
#   geom_ribbon(data = base_mort_df%>%filter(sex == 'Male'),
#               aes(ymin = min, ymax = max, x = year), fill = "darkseagreen1")+
#   #geom_line(data = base_mort_df%>%filter(sex == 'Male'), aes(x = year, y = val), color = 'darkgreen', size = .5)+
#   geom_point(data = base_mort_df_graph_male,
#              aes(x=year, y=value, group = variable, shape = variable), size = 1.5, color = 'black')+
#   labs(title=(bquote(atop("Baseline mortality rate calibration ranges, Males",
#                           mu[{1}]^{BASELINE}~(tau)~~"for"~all~tau~"in"~Year~Y))))+
#   scale_x_continuous(name = 'Year Y', breaks=seq(from = 1990, to = 2030, by = 10))+
#   scale_y_continuous(name = 'baseline mortality rate', limits = c(0, .012), breaks=(seq(0, .012, .003)))+
#   theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())+
#   scale_shape_manual(values=c(8, 17, 6))+
#   geom_vline(xintercept = 2017, linetype="dashed", 
#              color = "darkgrey", size=1)+
#   annotate("text", x=2010, y=.011, label= "calibration period", size = 5.5)+
#   annotate("text", x=2024, y=.011, label= "evaluation period", size = 5.5)

base_mort_df_graph_female<-base_mort_df%>%
  filter(sex == 'Female')%>%
  select(-c('sex'))

baseline_mort_plot_female<-ggplot(data=base_mort_df_graph_female, 
                                aes(x=year, y=val, group=1)) +
  geom_line(color="grey", size = 1)+
  geom_point(color="blue", size = 2)+
  geom_vline(xintercept = 2017, linetype="dashed", 
             color = "darkgrey", size=1)+
  annotate("text", x=2010, y=.011, label= "calibration period", size = 5.5)+
  annotate("text", x=2024, y=.011, label= "evaluation period", size = 5.5)+
  #labs(title=(bquote(atop("Baseline mortality rate values, Females"))))+
  scale_x_continuous(name = 'Year Y', breaks=seq(from = 1990, to = 2027, by = 4))+
  scale_y_continuous(name = bquote(atop(mu[{2}]^{VAL}~(Y))), limits = c(0, .012), breaks=(seq(0, .012, .003)))+
  theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))

# base_mort_df_graph_female<-melt(base_mort_df_graph_female, id = c('year'))
# 
# baseline_mort_plot_female<-ggplot()+
#   geom_ribbon(data = base_mort_df%>%filter(sex == 'Female'),
#               aes(ymin = min, ymax = max, x = year), fill = "plum2")+
#   geom_point(data = base_mort_df_graph_female,
#              aes(x=year, y=value, group = variable, shape = variable), size = 1.5, color = 'black')+
#   labs(title=(bquote(atop("Baseline mortality rate calibration ranges, Females",
#                           mu[{2}]^{BASELINE}~(tau)~~"for"~all~tau~"in"~Year~Y))))+
#   scale_x_continuous(name = 'Year Y', breaks=seq(from = 1990, to = 2030, by = 10))+
#   scale_y_continuous(name = 'baseline mortality rate', limits = c(0, .012), breaks=(seq(0, .012, .003)))+
#   theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.title=element_blank())+
#   scale_shape_manual(values=c(8, 17, 6))+
#   geom_vline(xintercept = 2017, linetype="dashed", 
#              color = "darkgrey", size=1)+
#   annotate("text", x=2010, y=.011, label= "calibration period", size = 5.5)+
#   annotate("text", x=2024, y=.011, label= "evaluation period", size = 5.5)

setwd(graph_outdir)

png("baseline_mort_plot_male.png", width = 540, height = 440)
print(baseline_mort_plot_male)
dev.off()

png("baseline_mort_plot_female.png", width = 540, height = 440)
print(baseline_mort_plot_female)
dev.off()

