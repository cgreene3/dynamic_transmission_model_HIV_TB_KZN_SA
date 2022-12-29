#births that CHANGE OVERTIME
#based on GBD
#clean an workspace
rm(list = ls())
gc()

#1940 pop init based on 1990 proportions (except for HIV)

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#'normalr'

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
#percent in each gender based on GBD pop estimates
indir <- paste0(here(),'/param_files/calculated_param_gen/raw_input_data/GBD')

#update parameter file
outdir <- paste0(here(),'/param_files/input_parameters')

#read data
setwd(indir)
pop_estimates_df<-read.csv('pop_estimates_15_19.csv')
all_prev_num_estimates_df<-read.csv('all_prev_num_15_19.csv')
HIV_prev<-read.csv('hiv_prev_perc_15_19.csv')

#aging into unifected (TB 1), latent recent (TB 3 and TB 4), active (TB 6)

#calculate # uninfected
uninfected_df_no_hiv<-all_prev_num_estimates_df%>%
  group_by(year, sex)%>%
  summarise(total_disease_pop = sum(val))%>%
  ungroup()%>%
  left_join(pop_estimates_df, by = c('year', 'sex'))%>%
  mutate(n_uninfected = expected_total_pop - total_disease_pop,
         TB_compartment = 1,
         DR_compartment = 1)%>%
  left_join(HIV_prev, by = c('year', 'sex'))%>%
  mutate(HIV_compartment = 1,
         HIV_adj = 1-val,
         G_compartment = if_else(sex == "Male", 1 , 2))%>%
  filter(!is.na(year))%>%
  select(c('year', 'TB_compartment', 'DR_compartment', 
           'HIV_compartment', 'G_compartment',
           'HIV_adj', 'n_uninfected', 'val'))


unif_hiv_pos<-uninfected_df_no_hiv%>%
  mutate(HIV_compartment = 2,
         HIV_adj = 1-HIV_adj)

uninfected_df<-rbind(uninfected_df_no_hiv,
                    unif_hiv_pos)

uninfected_df<-uninfected_df%>%
         mutate(pop_estimate = n_uninfected * HIV_adj)%>%
  select(c('year', 
           'TB_compartment', 
           'DR_compartment', 
           'HIV_compartment', 
           'G_compartment',
           'pop_estimate'))


recent_adj_calc<-2/15 #equally likely to have been infected with TB in first 15 years


#latent TB in GBD not seperated by DR or HIV status
latent_TB_df_recent<-all_prev_num_estimates_df%>%
  filter(cause == "Latent tuberculosis infection")%>%
  mutate(TB_compartment = 3,
         DR_compartment = 1, 
         DR_adj = 1-.036,
         recent_adj = recent_adj_calc)

latent_TB_df_recent<-rbind(latent_TB_df_recent,
                           latent_TB_df_recent%>%
          mutate(DR_compartment = 2,
                 DR_adj = .036)) #from WHO: TB country profile for SA (var epsilon)

latent_TB_df_remote<-all_prev_num_estimates_df%>%
  filter(cause == "Latent tuberculosis infection")%>%
  mutate(TB_compartment = 4,
         DR_compartment = 1, 
         DR_adj = 1-.036,
         recent_adj = 1-recent_adj_calc)

latent_TB_df_remote<-rbind(latent_TB_df_remote,
                           latent_TB_df_remote%>%
                             mutate(DR_compartment = 2,
                                    DR_adj = .036))

latent_TB_df<-rbind(latent_TB_df_recent, latent_TB_df_remote)

#add in HIV prev estimates
latent_TB_df<-latent_TB_df%>%
  left_join(HIV_prev, by = c('year', 'sex'))%>%
  mutate(HIV_compartment = 1,
         HIV_adj = 1-val.y)


latent_TB_df<-rbind(latent_TB_df,
                    latent_TB_df%>%
                      mutate(HIV_compartment = 2,
                             HIV_adj = val.y))

latent_TB_df<-latent_TB_df%>%
  mutate(G_compartment = if_else(sex == "Male", 1 , 2),
         pop_estimate = val.x * DR_adj * recent_adj* HIV_adj)%>%
  select(c('year', 
           'TB_compartment', 
           'DR_compartment', 
           'HIV_compartment', 
           'G_compartment',
         'pop_estimate'))

#calculate active
active_TB_df<-all_prev_num_estimates_df%>%
  filter(cause != "Latent tuberculosis infection")%>%
  mutate(HIV_compartment = if_else(grepl('HIV', cause),
                                   2, 1),
         DR_compartment = if_else(grepl('resistant', cause),
                                  2, 1),
         G_compartment = if_else(sex == "Male", 1 , 2),
         TB_compartment = 6)%>%
  group_by(year, TB_compartment, HIV_compartment, DR_compartment, G_compartment)%>%
  summarise(pop_estimate = sum(val))%>%
  filter(!is.na(year))


#combine
birth_perc_df_overtime_temp<-rbind(data.frame(uninfected_df),
                              data.frame(latent_TB_df),
                              data.frame(active_TB_df))


birth_perc_df<-data.frame(year = as.integer(),
                                   TB_compartment = as.integer(),
                                   DR_compartment = as.integer(),
                                   HIV_compartment = as.integer(),
                                   G_compartment = as.integer())


for (year in unique(pop_estimates_df$year)){
  for (t in 1:8){
    for (r in 1:2){
      for (h in 1:4){
        for (g in 1:2){
          birth_perc_df<-rbind(birth_perc_df,
                               c(year, t, r, h, g))
        }
      }
    }
  }
}

colnames(birth_perc_df)<-c('year', 'TB_compartment',
                           'DR_compartment', 'HIV_compartment',
                           'G_compartment')

birth_perc_df2<-birth_perc_df%>%
  left_join(birth_perc_df_overtime_temp, by = c('year', 'TB_compartment',
                                                'DR_compartment', 'HIV_compartment',
                                                'G_compartment'))%>%
  mutate(pop_estimate = if_else(is.na(pop_estimate), 0, pop_estimate))%>%
  group_by(year)%>%
  mutate(total_pop = sum(pop_estimate))%>%
  ungroup()%>%
  mutate(prop_of_pop = pop_estimate/total_pop)%>%
  select(c('year', 'TB_compartment', 'DR_compartment',
           'HIV_compartment', 'G_compartment', 'prop_of_pop'))

#used before 1990
prop_of_pop_no_hiv<-birth_perc_df2%>%
  filter(year == 1990)%>%
  group_by(TB_compartment, DR_compartment, G_compartment)%>%
  mutate(total_prop_no_HIV = sum(prop_of_pop))%>%
  ungroup()%>%
  mutate(prop_of_pop = if_else(HIV_compartment == 1, total_prop_no_HIV, 0))%>%
  mutate(year = 1989)%>%
  select(-c('total_prop_no_HIV'))

birth_perc_df3<-rbind(birth_perc_df2, prop_of_pop_no_hiv)

setwd(outdir)
write.csv(birth_perc_df3, 'birth_perc_df_overtime.csv')

#check everything adds up to 1
#test<-birth_perc_df3%>%group_by(year)%>%summarise(sum(prop_of_pop))

####pop init df###
pop_init_df_1940<-prop_of_pop_no_hiv%>%
  mutate(compartment_id = paste('N_', TB_compartment, "_", 
                                DR_compartment, "_", HIV_compartment, "_", 
                                G_compartment),
         dcompartment_id = paste0('d', compartment_id),
         total_pop = prop_of_pop*100000)%>%
  select(c('TB_compartment', 'DR_compartment', 
           'HIV_compartment', 'G_compartment',
           'compartment_id', 'dcompartment_id', 'total_pop'))


write.csv(pop_init_df_1940, 'pop_init_df_1940.csv')

