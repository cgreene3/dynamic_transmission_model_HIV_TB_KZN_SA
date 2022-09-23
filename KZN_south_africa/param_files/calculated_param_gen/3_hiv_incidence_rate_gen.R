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
graph_outdir<-paste0(here(),'/param_files/dynamic_param_graphs')

setwd(indir)

#read in pop estimates
pop_df<-read.csv('pop_estimates_all_ages.csv')

#read in incidence rate df
hiv_inc_df<-read.csv('HIV_inc_num_pulled_2_22_22.csv')%>%
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

hiv_inc_df2$max = hiv_inc_df2$val*1.25
hiv_inc_df2$min = hiv_inc_df2$val*.75

setwd(outdir)
write.csv(hiv_inc_df2, 'hiv_inc_df.csv', row.names = FALSE)

##assume sigmodal scale up between 1980-1990
#sigmoid = function(x) {
#  1 / (1 + exp(-x))
#}

#yrs_sigmodal_scale_up<-seq(0, (1990-1980), 1)
#sogmoid_estimates_perc<-sigmoid(yrs_sigmodal_scale_up)

