#initialize populations at start of warmup period

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr'), require, character.only=T)


#Make sure you have the epi_model_HIV_TB/KNZ_south_africa.Rproj open (top right)
indir<- paste0(here(), '/param_files/calculated_param_gen/input_data/GBD')
outdir <- paste0(here(),'/param_files')

#read in pop data
setwd(indir)
pop_df<-read.csv('pop_estimates.csv')

################ DEFINE SETS ###############################

#######8 TB set description (TB)#########
#1:Uninfected, not on IPT;
#2:Uninfected, on IPT; 
#3:LTBI, infected recently (within the past two-years)
#4: LTBI, infected remotely (more than two-years ago)
#5: LTBI, on IPT
#6: Active
#7: Recovered/Treated
#8: LTBI, after IPT

TB_SET<-1:8

######4 HIV compartments description (HIV)#########
#1 : HIV Negative
#2 : HIV Positive CD4 > 200 - No ART
#3 : HIV Positive CD4 =<: 200 - No Art
#4 : HIV Positive - ART 

HIV_SET<-1:4

######2 Drug Resistance compartments description (DR)#########
#1 : Drug Susceptible
#2 : Multi Drug Resistant

DR_SET<-1:2

#######2 Gender compartments description (G)########
#1: Male
#2: Female

G_SET<-1:2

#######Parameter extraction########

#pop init 1990
#establish compartment and dcompartment ids
#TB factors 
prop_recent<-0.049  #proprotion of LTBI that is recent
prop_remote<-0.951 #proportion remote
prop_IPT<-0#.01
prop_LTBI<-.49
LTBI_recent<-(prop_LTBI*prop_recent)*(1-prop_IPT)
LTBI_remote<-(prop_LTBI*prop_remote)*(1-prop_IPT)
LTBI_on_IPT<-(prop_LTBI*prop_IPT)
TB_active_male<-.01
TB_active_female<-.01
TB_uninfected_male<-1-.49-TB_active_male
TB_uninfected_female<-1-.49-TB_active_female
prop_DS<-1-0.037
prop_MDR<-0.037

TB_compartment<-c()
DR_compartment<-c()
G_compartment<-c()

tb_adj<-c()
dr_adj<-c()
hiv_adj<-c()
g_adj<-c()

#calculate gender init adjustments from GBD data
pop_df<-pop_df%>%
  group_by(year)%>%
  mutate(total_pop = sum(expected_total_pop))%>%
  ungroup()%>%
  mutate(percent_gender_pop = expected_total_pop/total_pop)%>%
  group_by(sex)%>%
  summarise(gender_percent = mean(percent_gender_pop))%>%
  mutate(sex_id = if_else(sex == 'Male', 1, 2))%>%
  arrange(sex_id)

gender_adjustments <- pop_df$gender_percent

for (tb in TB_SET){
  for (dr in DR_SET){
    for (g in G_SET){
      TB_compartment<-c(TB_compartment, tb)
      DR_compartment<-c(DR_compartment, dr)
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
      
      #gender adjustments
      g_adj<-c(g_adj, gender_adjustments[g])
    }
  }
}


pop_init_df<-data.frame(TB_compartment, DR_compartment, G_compartment, 
                        tb_adj, dr_adj, g_adj)

rm(tb_adj, dr_adj, g_adj, TB_compartment, DR_compartment, G_compartment)

pop_init_df<-pop_init_df%>%
  mutate(total_adj = round(tb_adj*dr_adj*g_adj,4),
         prop_of_pop = round(total_adj/sum(total_adj),5))%>%
  mutate(prop_of_pop = prop_of_pop*(1/sum(prop_of_pop)))%>% #fix rounding errors make sure add to 1
  mutate(total_pop = prop_of_pop*100000)%>%
  select(c('TB_compartment', 'DR_compartment', 'G_compartment', 'total_pop'))

#set the other HIV compartments to 0
pop_init_df_tempHIV<-rbind(pop_init_df, pop_init_df, pop_init_df)
pop_init_df_tempHIV$HIV_compartment <- rep(2:4, each = nrow(pop_init_df))
pop_init_df_tempHIV$total_pop<-rep(0, times = nrow(pop_init_df_tempHIV)) 

#combine dfs
pop_init_df$HIV_compartment<-rep(1, times = nrow(pop_init_df))
pop_init_df<-rbind(pop_init_df, pop_init_df_tempHIV)

pop_init_df<- pop_init_df%>%
  mutate(compartment_id = paste0("N_", TB_compartment, 
                                 "_", DR_compartment ,
                                 "_", HIV_compartment, 
                                 "_", G_compartment),
         dcompartment_id = paste0("dN_", TB_compartment,
                                  "_", DR_compartment ,"_",
                                  HIV_compartment, "_",
                                  G_compartment))

pop_init_df<-pop_init_df%>%
  arrange(G_compartment)%>%
  arrange(HIV_compartment)%>%
  arrange(DR_compartment)%>%
  arrange(TB_compartment)


setwd(outdir)
write.csv(pop_init_df, 'pop_init_df.csv', row.names = FALSE)
