#Estimating HIV incidence rate from GBD
#Uses GBD 2019
#(num HIV incidence)/pop estimate = incidence rate by gender
#KZN population estimates

#scales up linearly between 1980-1990 
#(when hiv incidence estimates don't exist from GBD)
#GBD from 1990-2017
#scales down according to reduced incidence projections from DO ART HPV-HIV model

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa for here to work
indir_GBD <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data/GBD')
indir_DOART<- paste0(here(), '/param_files/calculated_param_gen/raw_input_data/DO_ART')
outdir <- paste0(here(),'/param_files/input_parameters')
graph_outdir<-paste0(here(),'/param_files/dynamic_param_graphs')

setwd(indir_GBD)

#read in pop estimates
pop_df<-read.csv('pop_estimates_15_59.csv')

#read in incidence rate df
hiv_inc_df<-read.csv('hiv_inc_num.csv')%>%
  group_by(year, sex)%>%
  summarise(val = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex'))%>%
  mutate(val = (val/expected_total_pop))%>%
  select('year', 'sex', 'val')%>%
  ungroup()


hiv_inc_df_2007_less_male<-hiv_inc_df%>%
  filter(year <= 2004,
         sex == "Male")

hiv_inc_df_2007_less_female<-hiv_inc_df%>%
  filter(year <= 2004,
         sex == "Female")

male_1990_inc<-hiv_inc_df%>%
  filter(year == 1990,
         sex == 'Male')

female_1990_inc<-hiv_inc_df%>%
  filter(year == 1990,
         sex == 'Female')

mult_lm_male<-lm(c(hiv_inc_df_2007_less_male$val) ~ 
                   (c(male_1990_inc$val[1]*hiv_inc_df_2007_less_male$year)))

mult_lm_female<-lm(c(hiv_inc_df_2007_less_female$val) ~ 
                   (c(female_1990_inc$val[1]*hiv_inc_df_2007_less_female$year)))

male_rate<-mult_lm_male[["coefficients"]][2]
female_rate<-mult_lm_female[["coefficients"]][2]

est_yrs<-1989:1980

male_vals<-rep(0, times = length(est_yrs))
female_vals<-rep(0, times = length(est_yrs))

for (yr in 1:length(est_yrs)){
  male_vals[yr]<-male_1990_inc$val[1]/((1+male_rate)^yr)
  female_vals[yr]<-female_1990_inc$val[1]/((1+female_rate)^yr)
}

male_1980_1989_df<-data.frame(year = est_yrs,
                              sex = rep('Male', times = length(1980:1989)))
male_1980_1989_df$val = male_vals

female_1980_1989_df<-data.frame(year = est_yrs,
                              sex = rep('Female', times = length(1980:1989)))
female_1980_1989_df$val = female_vals

hiv_inc_df2<-rbind.data.frame(as.data.frame(hiv_inc_df),
                   as.data.frame(female_1980_1989_df),
                              as.data.frame(male_1980_1989_df))

hiv_inc_df2$program_id <- 1

hiv_inc_df2<-hiv_inc_df2%>%
  filter(year <= 2017)

#projections from 2018 to 2028 (evaluation years)
GBD_2017_female_val<-hiv_inc_df2%>%
  filter(sex == 'Female',
         year == 2017)

GBD_2017_male_val<-hiv_inc_df2%>%
  filter(sex == 'Male',
         year == 2017)

setwd(indir_DOART)
incidence_do_art_df_SOC<-readxl::read_excel('incidence_est.xlsx', sheet = 'option_1')%>%
  filter(year > 2017)%>%
  select(c('year', 'sex', 'Percent_Decrease_SOC'))%>%
  arrange(sex, year)%>%
  group_by(sex)%>%
  mutate(cum_perc_decrease = if_else(year == 2018, 
                                     1+Percent_Decrease_SOC,
                                     cumprod(1+Percent_Decrease_SOC)))%>%
  mutate(val = if_else(sex == "Female", 
                       (cum_perc_decrease*GBD_2017_female_val$val),
                       (cum_perc_decrease*GBD_2017_male_val$val)))%>%
  select(c('year', 'sex', 'val'))%>%
  ungroup()

incidence_do_art_df_SOC$program_id = 1
incidence_do_art_df_SOC<-as.data.frame(incidence_do_art_df_SOC)

incidence_do_art_df_DO_ART<-readxl::read_excel('incidence_est.xlsx', sheet = 'option_1')%>%
  filter(year > 2017)%>%
  select(c('year', 'sex', 'Percent_Decrease_DO_ART'))%>%
  arrange(sex, year)%>%
  group_by(sex)%>%
  mutate(cum_perc_decrease = if_else(year == 2018, 
                                     1+Percent_Decrease_DO_ART,
                                     cumprod(1+Percent_Decrease_DO_ART)))%>%
  mutate(val = if_else(sex == "Female", 
                       (cum_perc_decrease*GBD_2017_female_val$val),
                       (cum_perc_decrease*GBD_2017_male_val$val)))%>%
  select(c('year', 'sex', 'val'))%>%
  ungroup()

incidence_do_art_df_DO_ART<-as.data.frame(incidence_do_art_df_DO_ART)

incidence_do_art_df_DO_ART_2<-incidence_do_art_df_DO_ART
incidence_do_art_df_DO_ART_2$program_id = 2

incidence_do_art_df_DO_ART_3<-incidence_do_art_df_DO_ART
incidence_do_art_df_DO_ART_3$program_id = 3

hiv_inc_df3<-rbind(hiv_inc_df2, incidence_do_art_df_SOC, incidence_do_art_df_DO_ART_2,
                   incidence_do_art_df_DO_ART_3)


hiv_inc_df3$max = hiv_inc_df3$val*1.25
hiv_inc_df3$min = hiv_inc_df3$val*.75

hiv_inc_df3$program_id<-as.factor(hiv_inc_df3$program_id)


#########graphing just val
hiv_inc_df_graph_male<-hiv_inc_df3%>%
  filter(sex == 'Male')%>%
  select(-c('sex'))%>%
  filter(program_id != 3)%>%
  mutate(`Policy ID` = if_else(program_id == 1, "Program 1", "Program 2 and 3"))

hiv_inc_plot_male<-ggplot(data=hiv_inc_df_graph_male) +
  geom_line(aes(x=year, y=val, group = factor(program_id)), color="grey")+
  geom_point(aes(x=year, y=val, colour=factor(`Policy ID`), 
                 shape = factor(`Policy ID`)), size = 1.5)+
  scale_color_manual(values = c("blue", "red"))+
  geom_vline(xintercept = 2017, linetype="dashed", 
             color = "darkgrey", size=1)+
  annotate("text", x=2009.5, y=.045, label= "calibration period", size = 5)+
  annotate("text", x=2023.5, y=.045, label= "evaluation period", size = 5)+
  #labs(title="HIV incidence rate, Males")+
  scale_x_continuous(name = 'Year Y', breaks=seq(from = 1980, to = 2028, by = 4))+
  scale_y_continuous(name = bquote(atop(eta[{"1,2,1"}]^{VAL}~(Y))), limits = c(0, .05), breaks=(seq(0, .05, .01)))+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5), legend.position = "top",
        legend.key = element_rect(fill = "grey94"),
        legend.background = element_rect(colour = "lightgrey"),
        legend.box.background = element_rect(colour = "black"))+
  scale_shape_manual(values=c(8, 2))

hiv_inc_df_graph_female<-hiv_inc_df3%>%
  filter(sex == 'Female')%>%
  select(-c('sex'))%>%
  filter(program_id != 3)%>%
  mutate(`Policy ID` = if_else(program_id == 1, "Program 1", "Program 2 and 3"))

hiv_inc_plot_female<-ggplot(data=hiv_inc_df_graph_female) +
  geom_line(aes(x=year, y=val, group = factor(program_id)), color="grey")+
  geom_point(aes(x=year, y=val, colour=factor(`Policy ID`), 
                 shape = factor(`Policy ID`)), size = 1.3)+
  scale_color_manual(values = c("blue", "red"))+
  geom_vline(xintercept = 2017, linetype="dashed", 
             color = "darkgrey", size=1)+
  annotate("text", x=2009.5, y=.045, label= "calibration period", size = 5)+
  annotate("text", x=2023.5, y=.045, label= "evaluation period", size = 5)+
  #labs(title="HIV incidence rate, Females")+
  scale_x_continuous(name = 'Year Y', breaks=seq(from = 1980, to = 2028, by = 4))+
  scale_y_continuous(name = bquote(atop(eta[{"1,2,2"}]^{VAL}~(Y))), limits = c(0, .05), breaks=(seq(0, .05, .01)))+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5), legend.position = "top",
        legend.key = element_rect(fill = "grey94"),
        legend.background = element_rect(colour = "lightgrey"),
        legend.box.background = element_rect(colour = "black"))+
  scale_shape_manual(values=c(8, 2))



#reshape data for ribbon #if max and min
# hiv_inc_ribbon_graph<-hiv_inc_df3%>%
#   filter(program_id != 3)%>%
#   mutate(`policy id` = if_else(program_id == 1, 'Policy 1', 'Policy 2 and 3'))
# 
# hiv_inc_ribbon_graph_female<-hiv_inc_ribbon_graph%>%filter(sex == 'Female')
# hiv_inc_ribbon_graph_male<-hiv_inc_ribbon_graph%>%filter(sex == 'Male')
# 
# hiv_inc_df_graph_female<-hiv_inc_df3%>%
#   filter(sex == 'Female')%>%
#   select(-c('sex'))
# 
# hiv_inc_df_graph_male<-hiv_inc_df3%>%
#   filter(sex == 'Male')%>%
#   select(-c('sex'))
# 
# hiv_inc_df_point_graph_female<-melt(hiv_inc_df_graph_female, id = c('year', 'program_id'))
# hiv_inc_df_point_graph_female<-hiv_inc_df_point_graph_female%>%
#   filter(program_id != 3)%>%
#   mutate(`policy id` = if_else(program_id == 1, 'Policy 1', 'Policies 2 and 3'))
# 
# hiv_inc_df_point_graph_male<-melt(hiv_inc_df_graph_male, id = c('year', 'program_id'))
# hiv_inc_df_point_graph_male<-hiv_inc_df_point_graph_male%>%
#   filter(program_id != 3)%>%
#   mutate(`policy id` = if_else(program_id == 1, 'Policy 1', 'Policies 2 and 3'))
# 
# 
# hiv_inc_plot_female<-ggplot()+
#   geom_ribbon(data = hiv_inc_ribbon_graph_female,
#               aes(ymin = min, ymax = max, 
#                   x = year, fill = `policy id`))+
#   scale_fill_manual(values = c("plum2", "purple"), name = "policy id")+
#   geom_point(data = hiv_inc_df_point_graph_female,
#              aes(x = year, y=value, group = variable, shape = variable), 
#              size = 1.5, color = "black")+
#   labs(title=(bquote(atop("HIV incidence rate calibration ranges, Females", 
#                           eta[{"1,2,2"}]~(tau)~~"for"~all~tau~"in"~Year~Y))))+
#   scale_x_continuous(name = 'Year Y', breaks=seq(from = 1980, to = 2030, by = 10))+
#   scale_y_continuous(name = 'HIV incidence rate', limits = c(0, .05), breaks=(seq(0, .05, .01)))+
#   theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.title=element_blank())+
#   scale_shape_manual(values=c(8, 17, 6))+
#   geom_vline(xintercept = 2017, linetype="dashed", 
#              color = "darkgrey", size=1)+
#   annotate("text", x=2009, y=.045, label= "calibration period", size = 5)+
#   annotate("text", x=2024.5, y=.045, label= "evaluation period", size = 5)
# 
# hiv_inc_plot_male<-ggplot()+
#   geom_ribbon(data = hiv_inc_ribbon_graph_male,
#               aes(ymin = min, ymax = max, 
#                   x = year, fill = `policy id`))+
#   scale_fill_manual(values = c("darkseagreen1", "darkgreen"), name = "policy id")+
#   geom_point(data = hiv_inc_df_point_graph_male,
#              aes(x = year, y=value, group = variable, shape = variable), 
#              size = 1.5, color = "black")+
#   labs(title=(bquote(atop("HIV incidence rate calibration ranges, Males", 
#                           eta[{"1,2,1"}]~(tau)~~"for"~all~tau~"in"~Year~Y))))+
#   scale_x_continuous(name = 'Year Y', breaks=seq(from = 1980, to = 2030, by = 10))+
#   scale_y_continuous(name = 'HIV incidence rate', limits = c(0, .05), breaks=(seq(0, .05, .01)))+
#   theme(text = element_text(size=16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.title=element_blank())+
#   scale_shape_manual(values=c(8, 17, 6))+
#   geom_vline(xintercept = 2017, linetype="dashed", 
#              color = "darkgrey", size=1)+
#   annotate("text", x=2009, y=.045, label= "calibration period", size = 5)+
#   annotate("text", x=2024.5, y=.045, label= "evaluation period", size = 5)

setwd(outdir)
write.csv(hiv_inc_df3, 'hiv_inc_df.csv', row.names = FALSE)

setwd(graph_outdir)
png("hiv_inc_plot_female.png", width = 600, height = 440)
print(hiv_inc_plot_female)
dev.off()


png("hiv_inc_plot_male.png", width =600, height = 440)
print(hiv_inc_plot_male)
dev.off()
