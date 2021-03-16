#Pop init 1990

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
#note: only hiv prev changing over time according to caras model

hiv_indir<-paste0(here(),'/param_files/hiv_param_gen')
outdir <- paste0(here(),'/param_files')

#to name plots according to date generated
data_gen_date<-Sys.Date()
data_gen_date<-str_replace_all(data_gen_date, '-', '_')

#TB factors (don't change overtime)
prop_recent<-0.049  #proprotion of LTBI that is recent
prop_remote<-0.951 #proportion remote
prop_IPT<-.01
prop_LTBI<-.49
LTBI_recent<-(prop_LTBI*prop_recent)*(1-prop_IPT)
LTBI_remote<-(prop_LTBI*prop_remote)*(1-prop_IPT)
LTBI_on_IPT<-(prop_LTBI*prop_IPT)
TB_active_male<-.006
TB_active_female<-.012
TB_uninfected_male<-1-.49-TB_active_male
TB_uninfected_female<-1-.49-TB_active_female

#MDR factors (don't change overtime)
prop_DS<-.963
prop_MDR<-.037

setwd(hiv_indir)
#hiv adjustments overtime
hiv_prop_df<-read_excel('hiv_input_gen_data.xlsx')%>%
  filter(year == 1990)

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
            tb_adj <- c(tb_adj, TB_uninfected_male*(1-prop_IPT))
          } else {
            tb_adj<- c(tb_adj, TB_uninfected_female*(1-prop_IPT))
          }
        } else if(tb == 2){
          if (g == 1){
            tb_adj <- c(tb_adj, TB_uninfected_male*(prop_IPT))
          } else {
            tb_adj<- c(tb_adj, TB_uninfected_female*(prop_IPT))
          }
        }
        else if (tb == 3) {
          tb_adj <- c(tb_adj, LTBI_recent)
        } else if (tb == 4){
          tb_adj <- c(tb_adj, LTBI_remote)
        } else if (tb == 5){
          tb_adj <-c(tb_adj, LTBI_on_IPT)
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
        g_temp <- if_else(g == 1, 'Males', 'Females')
        if (hiv == 1){
          hiv_adj <- c(hiv_adj, (unique(hiv_prop_df$n1_prop[(hiv_prop_df$gender== g_temp)])))
        } else if (hiv == 2){
          hiv_adj <- c(hiv_adj, (unique(hiv_prop_df$n2_prop[(hiv_prop_df$gender== g_temp)])))
        } else if (hiv == 3){
          hiv_adj <- c(hiv_adj, (unique(hiv_prop_df$n3_prop[(hiv_prop_df$gender== g_temp)])))
        } else {
          hiv_adj <- c(hiv_adj, (unique(hiv_prop_df$n4_prop[(hiv_prop_df$gender== g_temp)])))
        }
      }
    }
  }
}

pop_init_df<-data.frame(TB_compartment, DR_compartment,
                          HIV_compartment, G_compartment, 
                          tb_adj, dr_adj, hiv_adj)

rm(tb_adj, hiv_adj, dr_adj)

pop_init_df<-pop_init_df%>%
  mutate(total_adj = round(tb_adj*dr_adj*hiv_adj,4),
         prop_of_pop = round(total_adj/sum(total_adj),5))%>%
  mutate(prop_of_pop = prop_of_pop*(1/sum(prop_of_pop)))%>% #fix rounding errors make sure add to 1
  mutate(total_pop = prop_of_pop*100000)%>%
  select(c('TB_compartment', 'DR_compartment', 
           'HIV_compartment', 'G_compartment',
           'total_pop'))


setwd(outdir)
write.csv(pop_init_df, 'pop_init_df.csv', row.names = FALSE)