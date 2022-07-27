#births that DO NOT CHANGE OVERTIME
#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr', 'normalr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
#percent in each gender based on GBD pop estimates
pop_estimates_indir<-paste0(here(), '/param_files/calculated_param_gen/input_data/GBD_pop_estimates')

#update parameter file
outdir <- paste0(here(),'/param_files/')

#TB factors (don't change overtime)
prop_recent<-0.049  #proportion of LTBI that is recent
prop_remote<-1-prop_recent #proportion remote
LTBI_recent<-.49*prop_recent
LTBI_remote<-.49*prop_remote
TB_active_male<-.01#.006
TB_active_female<-.01#.012
TB_uninfected_male<-1-.49-TB_active_male
TB_uninfected_female<-1-.49-TB_active_female

#MDR factors (don't change overtime)
prop_DS<-.963
prop_MDR<-.037

#HIV+ based on the HSRC survey from 2017
HIV_infected<-c(.03, .038)

setwd(pop_estimates_indir)
#gender births overtime
pop_df<-read.csv('pop_df.csv')

pop_df<-pop_df%>%
  group_by(year)%>%
  mutate(total_pop = sum(expected_total_pop))%>%
  ungroup()%>%
  mutate(percent_gender_pop = expected_total_pop/total_pop)%>%
  group_by(sex_id)%>%
  summarise(gender_percent_births = mean(percent_gender_pop))%>%
  arrange(sex_id)

gender_adjustments <- pop_df$gender_percent_births

TB_compartments <- 1:8
HIV_compartments<-1:4
DR_compartments<-1:2
G_compartments<-1:2

TB_compartment<-c()
DR_compartment<-c()
HIV_compartment<-c()
G_compartment<-c()

tb_adj<-c()
dr_adj<-c()
g_adj<-c()
hiv_adj<-c()

for (tb in TB_compartments){
  for (dr in DR_compartments){
    for (hiv in HIV_compartments){
      for (g in G_compartments){
        TB_compartment<-c(TB_compartment, tb)
        DR_compartment<-c(DR_compartment, dr)
        HIV_compartment<-c(HIV_compartment, hiv)
        G_compartment<-c(G_compartment, g)
        
        #tb adjustments
        if (tb == 1){
          if (g == 1){
            tb_adj <- c(tb_adj, TB_uninfected_male)
            } else {
            tb_adj<- c(tb_adj, TB_uninfected_female)
            }
        } else if (tb == 3) {
          tb_adj <- c(tb_adj, LTBI_recent)
        } else if (tb == 4){
          tb_adj <- c(tb_adj, LTBI_remote)
        } else if (tb == 6){
          if (g == 1){
            tb_adj <- c(tb_adj, TB_active_male)
           } else {
             tb_adj<- c(tb_adj, TB_active_female)
            }
        } else {
          tb_adj<-c(tb_adj, 0) #no one is birthed to tb compartments 2, 5, 7, or 8
        }
        
        #dr adjustments
        if (tb == 1){ #all tb unifected is assigned to dr compartment 1
          if (dr == 1){
            dr_adj<-c(dr_adj, 1)
          } else {
            dr_adj<-c(dr_adj, 0)
          }
        } else {
          if (dr == 1){
            dr_adj<-c(dr_adj, prop_DS)
          } else{
            dr_adj<-c(dr_adj, prop_MDR)
          }
        }
        
        #hiv adjustments
        hiv_adj_temp<-if_else(hiv == 1, 1-HIV_infected[g],
                if_else(hiv == 2, HIV_infected[g],
                        0))
        
        hiv_adj<-c(hiv_adj, hiv_adj_temp)
        
        #gender adjustments
        g_adj<-c(g_adj, gender_adjustments[g])
        
      }
    }
  }
}

birth_perc_df<-data.frame(TB_compartment, DR_compartment,
                          HIV_compartment, G_compartment, 
                          tb_adj, dr_adj, hiv_adj, g_adj)

birth_perc_df<-birth_perc_df%>%
  mutate(total_adj = tb_adj*dr_adj*hiv_adj*g_adj)%>%
  mutate(prop_of_pop = total_adj/sum(total_adj))

######### #to make sure all proportions add to 1 (fixing rounding error) currently not running#########
# test<-birth_perc_df%>%
#   group_by(year, G_compartment)%>%
#   summarise(total_birth_perc = sum(prop_of_pop))
# 
# birth_perc_df<-birth_perc_df%>%
#   left_join(test, by = c('year'))%>%
#   mutate(prop_of_pop = as.double(round(prop_of_pop*(1/total_birth_perc),2)))%>%
#   select(c('year', 'TB_compartment', 'DR_compartment',
#            'HIV_compartment', 'G_compartment', 'prop_of_pop'))%>%
#   mutate(prop_of_pop2 = if_else((TB_compartment == 1) & 
#                                  (HIV_compartment == 1) & 
#                                  (DR_compartment == 1) &
#                                  (G_compartment == 1), prop_of_pop + 1-(sum(prop_of_pop)), 
#                                prop_of_pop))%>%
#   mutate(prop_of_pop3 = as.double(if_else((TB_compartment == 1) & 
#                                   (HIV_compartment == 1) & 
#                                   (DR_compartment == 2) &
#                                   (G_compartment == 1), prop_of_pop2 + 1-(sum(prop_of_pop2)), 
#                                 prop_of_pop2)))%>%
#   mutate(prop_of_pop4 = if_else(prop_of_pop3 < 0, 0, prop_of_pop3))%>%
#   select(c('year', 'TB_compartment', 'DR_compartment',
#            'HIV_compartment', 'G_compartment', 'prop_of_pop4'))%>%
#   rename(prop_of_pop = prop_of_pop4)

#######write file to param####
setwd(outdir)
write.csv(birth_perc_df, 'birth_perc_df.csv', row.names = FALSE)
