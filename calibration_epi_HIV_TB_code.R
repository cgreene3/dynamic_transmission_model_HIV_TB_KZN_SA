#This file takes gbd prevalence estimates from a csv generated using the results tool 
#http://ghdx.healthdata.org/gbd-results-tool and aggregates them to use in model calibration.
#Data are already
#July 28, 2020

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 'readxl', 'stringr', 'reshape2', 'ggplot2', 'varhandle', 'here', 'paletteer'), require, character.only=T)

indir<-'param_files/'
outdir<-'test_outputs/calibration'

df<-read.csv(here(paste0(indir, "GBD_prev_1990_2017_jul28.csv")))

#Section 1 - Assembling calibration dataset. For each year (1990 to 2017) generate TB prevalence rate (per 100,000 population) by HIV+, HIV-, and male/female

#Group 1 - Active TB among HIV-negative males. One value per year (1990 to 2017). Use the number in "val"
    #For each year, sum across these three estimates: sex_id==1 & (cause_id==934|cause_id==946|cause_id==947)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 2  - Active TB prevalence among HIV-negative females
    #For each year, sum across these three estimates: sex_id==2 & (cause_id==934|cause_id==946|cause_id==947)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 3 - Active TB HIV+ males
    #For each year, sum across these three estimates: sex_id==1 & (cause_id==948|cause_id==949|cause_id==950)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 4 - Active TB HIV+ females
    #For each year, sum across these three estimates: sex_id==2 & (cause_id==948|cause_id==949|cause_id==950)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 5 - LTBI prevalence over time for males
    #Use sex_id==1 and cause_id==954
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 6 - LTBI prevalence over time for females
    #Use sex_id==2 and cause_id==954
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

df<-df%>%
  mutate(group_id = if_else(sex_id == 1 & cause_id == 934|
         sex_id == 1 & cause_id == 946|
         sex_id == 1 & cause_id == 947, 1,
         if_else(sex_id == 2 & cause_id == 934|
         sex_id == 2 & cause_id == 946|
         sex_id == 2 & cause_id == 947, 2,
         if_else(sex_id == 1 & cause_id == 948|
         sex_id == 1 & cause_id == 949|
         sex_id == 1 & cause_id == 950, 3,
         if_else(sex_id == 2 & cause_id == 948|
         sex_id == 2 & cause_id == 949|
         sex_id == 2 & cause_id == 950, 4,
         if_else(sex_id == 1 & cause_id == 954, 5,
                 if_else(sex_id == 2 & cause_id == 954, 6,
                         100)))))))

#filter to only include relevant values
df<-df%>%
  filter(group_id != 100)%>%
  filter(measure_id == 5 & metric_id == 3)%>%
  group_by(group_id, year)%>%
  summarise(expected = sum(val),
            upper = sum(upper),
            lower = sum(lower))

#Section 2 - Line graph of values in groups 1-4 with years on x axis

df_active<-df%>%
  filter(group_id == 1|
           group_id == 2|
           group_id == 3|
           group_id == 4)

df_active$compartment <- if_else(df_active$group_id == 1, 'HIV-, Male',
                                 if_else(df_active$group_id == 2, 'HIV-, Female',
                                         if_else(df_active$group_id == 3, 'HIV+, Male',
                                                 'HIV+, Female'))) 

df_active$group_id<-as.factor(df_active$group_id)

active_prev_overtime <- ggplot(df_active, aes(x = year, 
                                              y = expected,
                                              group = compartment,
                                              color = compartment))+
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                 position=position_dodge(0.05))+
  labs(title = "active prevalence overtime",
       y = 'prevalence') +
  scale_color_manual(values=c(paletteer_dynamic("cartography::green.pal", 4)))

print(active_prev_overtime)
 
#Section 3 - Line graph of values in groups 5&6 with years on x-axis. Values should be substantially higher than groups 1-4.

df_ltbi<-df%>%
  filter(group_id == 5|
           group_id == 6)

df_ltbi$group_name <- if_else(df_ltbi$group_id == 5, 'Male', 'Female')

df_ltbi$group_id<-as.factor(df_ltbi$group_id)

ltbi_prev_overtime <- ggplot(df_ltbi, aes(x = year, 
                                              y = expected,
                                              group = group_name,
                                              color = group_name))+
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                 position=position_dodge(0.05))+
  labs(title = "LTBI prevalence overtime",
       y = 'prevalence') +
  scale_color_manual(values=c(paletteer_dynamic("cartography::blue.pal", 2)))

print(ltbi_prev_overtime)

