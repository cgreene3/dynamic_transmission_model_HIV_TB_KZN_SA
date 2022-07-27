#births that CHANGE OVERTIME
#based on GBD
#clean an workspace
rm(list = ls())
gc()

#1940 pop init based on 1990 proportions (except for HIV)

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr', 'normalr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
#percent in each gender based on GBD pop estimates
indir <- paste0(here(),'/param_files/calculated_param_gen/input_data/GBD')

#update parameter file
outdir <- paste0(here(),'/param_files/')

#read data
setwd(indir)
pop_estimates_df<-read.csv('pop_estimates_15_19.csv')
all_prev_num_estimates_df<-read.csv('all_prev_num_15_19.csv')
HIV_prev<-read.csv('hiv_prev_perc_15_19.csv')

#group prev estimates according to model compartments
#aging into unifected (TB 1), latent recent (TB 3)...since young!, active (TB 6)

#calculate # uninfected
unifected_df<-all_prev_num_estimates_df%>%
  group_by(year, sex_name)%>%
  summarise(total_disease_pop = sum(val))%>%
  left_join(pop_estimates_df, by = c('year', 'sex_name' = 'sex'))%>%
  mutate(n_uninfected = expected_total_pop - total_disease_pop,
         TB_compartment = 1)%>%
  left_join(HIV_prev, by = c('year', 'sex_name' = 'sex'))%>%
  mutate(HIV_compartment = 1,
         HIV_adj = 1-val)%>%
  filter(!is.na(year))

unifected_df<-rbind(unifected_df,
                    unifected_df%>%
                      mutate(HIV_compartment = 2,
                             HIV_adj = val))

unifected_df<-unifected_df%>%
  mutate(G_compartment = if_else(sex_name == "Male", 1 , 2),
         pop_estimate = n_uninfected * HIV_adj,
         DR_compartment = 1)%>%
  select(c('year', 
           'TB_compartment', 
           'DR_compartment', 
           'HIV_compartment', 
           'G_compartment',
           'pop_estimate'))

#latent TB in GBD not seperated by DR or HIV status
latent_TB_df<-all_prev_num_estimates_df%>%
  filter(cause_name == "Latent tuberculosis infection")%>%
  mutate(TB_compartment = 3,
         DR_compartment = 1, 
         DR_adj = 1-.036)

latent_TB_df<-rbind(latent_TB_df,
                    latent_TB_df%>%
          mutate(DR_compartment = 2,
                 DR_adj = .036)) #from WHO: TB country profile for SA (var epsilon)

#add in HIV prev estimates
latent_TB_df<-latent_TB_df%>%
  left_join(HIV_prev, by = c('year', 'sex_name' = 'sex'))%>%
  mutate(HIV_compartment = 1,
         HIV_adj = 1-val.y)

latent_TB_df<-rbind(latent_TB_df,
                    latent_TB_df%>%
                      mutate(HIV_compartment = 2,
                             HIV_adj = val.y))

latent_TB_df<-latent_TB_df%>%
  mutate(G_compartment = sex_id,
         pop_estimate = val.x * DR_adj * HIV_adj)%>%
  select(c('year', 
           'TB_compartment', 
           'DR_compartment', 
           'HIV_compartment', 
           'G_compartment',
         'pop_estimate'))

#calculate active
active_TB_df<-all_prev_num_estimates_df%>%
  filter(cause_name != "Latent tuberculosis infection")%>%
  mutate(HIV_compartment = if_else(grepl('HIV', cause_name),
                                   2, 1),
         DR_compartment = if_else(grepl('resistant', cause_name),
                                  2, 1),
         G_compartment = sex_id,
         TB_compartment = 6)%>%
  group_by(year, HIV_compartment, DR_compartment, G_compartment)%>%
  summarise(pop_estimate = sum(val))%>%
  filter(!is.na(year))

#combine
birth_perc_df_overtime_temp<-rbind(unifected_df,
                              latent_TB_df,
                              active_TB_df)

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

birth_perc_df<-birth_perc_df%>%
  left_join(birth_perc_df_overtime_temp, by = c('year', 'TB_compartment',
                                                'DR_compartment', 'HIV_compartment',
                                                'G_compartment'))%>%
  mutate(pop_estimate = if_else(is.na(pop_estimate), 0, pop_estimate))%>%
  group_by(year)%>%
  mutate(total_pop = sum(pop_estimate))%>%
  ungroup()%>%
  mutate(prop_of_pop = pop_estimate/total_pop)

setwd(outdir)
write.csv(birth_perc_df, 'birth_perc_df_overtime.csv')

####pop init df###
pop_init_df_1940<-birth_perc_df%>%
  filter(year == 1990)%>%
  group_by(TB_compartment, DR_compartment, G_compartment)%>%
  mutate(total_prop_no_HIV = sum(prop_of_pop))%>%
  ungroup()%>%
  mutate(compartment_id = paste0('N_', t, "_", r, "_", h, "_", g),
         dcompartment_id = paste0('d', compartment_id),
         total_pop = if_else(HIV_compartment == 1, 
                             total_prop_no_HIV,
                             0)*100000)%>%
  select(c('TB_compartment', 'DR_compartment', 
           'HIV_compartment', 'G_compartment',
           'compartment_id', 'dcompartment_id', 'total_pop'))


write.csv(pop_init_df_1940, 'pop_init_df_1940.csv')
