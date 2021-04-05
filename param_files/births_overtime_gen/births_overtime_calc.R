#births overtime 1990-2017

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr', 'normalr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
#note: only hiv prev changing over time according to caras model
#percent in each gender based on GBD pop estimates

hiv_indir<-paste0(here(),'/param_files/hiv_param_gen')
pop_estimates_indir<-paste0(here(), '/param_files/GBD_pop_estimates')
outdir <- paste0(here(),'/param_files')

#to name plots according to date generated
data_gen_date<-Sys.Date()
data_gen_date<-str_replace_all(data_gen_date, '-', '_')

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

setwd(hiv_indir)
#hiv adjustments overtime
hiv_prop_df<-read_excel('hiv_input_gen_data.xlsx')

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
        
        #gender adjustments
        g_adj<-c(g_adj, gender_adjustments[g])
        
      }
    }
  }
}

birth_rate_df<-data.frame(TB_compartment, DR_compartment,
                          HIV_compartment, G_compartment, 
                          tb_adj, dr_adj, g_adj)

birth_rate_df<-do.call('rbind', 
                       replicate(2017-1990+1, 
                                 birth_rate_df, simplify = FALSE))

birth_rate_df$year <- rep(1990:2017, each = 128)

hiv_adj<-c()

for (row in 1:nrow(birth_rate_df)){
  yr_temp<- birth_rate_df[row, 'year']
  g_temp <- if_else(birth_rate_df[row, 'G_compartment'] == 1,
                    'Males', 
                    'Females')
  if (birth_rate_df[row, 'HIV_compartment'] == 1){
    hiv_adj <- c(hiv_adj, (unique(hiv_prop_df$n1_prop[(hiv_prop_df$year == yr_temp) &
                                                      (hiv_prop_df$gender== g_temp)])))
  } else if (birth_rate_df[row, 'HIV_compartment'] == 2){
    hiv_adj <- c(hiv_adj, (unique(hiv_prop_df$hiv_prevalence[(hiv_prop_df$year == yr_temp) &
                                                      (hiv_prop_df$gender== g_temp)])))
  } else{
    hiv_adj<-c(hiv_adj, 0)
  }
}

birth_rate_df$hiv_adj <- hiv_adj

rm(tb_adj, hiv_adj, dr_adj, gender_adjustments, g_adj)

birth_rate_df<-birth_rate_df%>%
  mutate(total_adj = tb_adj*dr_adj*hiv_adj*g_adj)%>%
  group_by(year)%>%
  mutate(prop_of_pop = total_adj/sum(total_adj))

# #to make sure all proportions add to 1 (fixing rounding error)
# test<-birth_rate_df%>%
#   group_by(year, G_compartment)%>%
#   summarise(total_birth_perc = sum(prop_of_pop))
# 
# birth_rate_df<-birth_rate_df%>%
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



setwd(outdir)
write.csv(birth_rate_df, 'birth_rate_df.csv', row.names = FALSE)
