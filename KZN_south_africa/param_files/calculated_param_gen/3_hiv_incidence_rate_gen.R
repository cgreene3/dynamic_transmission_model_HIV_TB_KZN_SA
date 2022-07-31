#Estimating HIV incidence rate from GBD
#Uses GBD 2019
#(num HIV incidence)/pop estimate = incidence rate by gender
#KZN population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to KZN_south_africa for here to work
indir <- paste0(here(),'/param_files/calculated_param_gen/input_data/GBD')
outdir <- paste0(here(),'/param_files/input_parameters')

setwd(indir)

#read in pop estimates
pop_df<-read.csv('pop_estimates_all_ages.csv')

#read in incidence rate df
hiv_inc_df<-read.csv('HIV_inc_num_pulled_2_22_22.csv')%>%
  group_by(year, sex)%>%
  summarise(val = sum(val))%>%
  left_join(pop_df, by = c('year', 'sex'))%>%
  mutate(val = (val/expected_total_pop))%>%
  select('year', 'sex', 'val')

#assume acceleration between 1980-1990 as seen between 1990 and 1991
male_1990_inc<-hiv_inc_df%>%
  filter(year == 1990,
         sex == 'Male')

male_1991_inc<-hiv_inc_df%>%
  filter(year == 1991,
         sex == 'Male')

male_scale <-male_1991_inc$val/male_1990_inc$val
male_vals<-rep(0, times = length(1989:1980))

for (yr in 1989:1980){
  if (yr == 1989){
    male_vals[1990-yr]=male_1990_inc$val/male_scale
  } else {
    male_vals[1990-yr]=male_vals[1990-yr-1]/male_scale
  }
}

male_1980_1989_df<-data.frame(year = 1989:1980,
                              sex = rep('Male', times = length(1980:1989)))
male_1980_1989_df$val = male_vals

##females
female_1990_inc<-hiv_inc_df%>%
  filter(year == 1990,
         sex == 'Female')

female_1991_inc<-hiv_inc_df%>%
  filter(year == 1991,
         sex == 'Female')

female_scale <-female_1991_inc$val/female_1990_inc$val
female_vals<-rep(0, times = length(1989:1980))

for (yr in 1989:1980){
  if (yr == 1989){
    female_vals[1990-yr]=female_1990_inc$val/female_scale
  } else {
    female_vals[1990-yr]=female_vals[1990-yr-1]/female_scale
  }
}

female_1980_1989_df<-data.frame(year = 1989:1980,
                              sex = rep('Female', times = length(1980:1989)))
female_1980_1989_df$val = female_vals

hiv_inc_df2<-rbind(data.frame(hiv_inc_df),
                   data.frame(female_1980_1989_df),
                              data.frame(male_1980_1989_df))

hiv_inc_df2$max = hiv_inc_df2$val*1.25
hiv_inc_df2$min = hiv_inc_df2$val*.75

setwd(outdir)
write.csv(hiv_inc_df2, 'hiv_inc_df.csv', row.names = FALSE)
