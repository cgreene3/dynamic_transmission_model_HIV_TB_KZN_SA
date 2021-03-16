#pop init 1990

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/pop_init_gen')
outdir <- paste0(here(),'/param_files')

#to name plots according to date generated
data_gen_date<-Sys.Date()
data_gen_date<-str_replace_all(data_gen_date, '-', '_')

#read in param init 
setwd(indir)
pop_init_gen_df<-read_excel('pop_init_gen_1990_data.xlsx')%>%
  filter(adj_input == 'yes')

TB_compartments <- 1:8
HIV_compartments<-1:4
DR_compartments<-1:2
G_compartments<-1:2

TB_compartment<-c()
DR_compartment<-c()
HIV_compartment<-c()
G_compartment<-c()

for (tb in TB_compartments){
  for (dr in DR_compartments){
    for (hiv in HIV_compartments){
      for (g in G_compartments){
        TB_compartment<-c(TB_compartment, tb)
        DR_compartment<-c(DR_compartment, dr)
        HIV_compartment<-c(HIV_compartment, hiv)
        G_compartment<-c(G_compartment, g)
      }
    }
  }
}

pop_init_df<-data.frame(TB_compartment,DR_compartment,HIV_compartment,G_compartment)
tb_adj<-rep(0, times = nrow(pop_init_df))
hiv_adj<-rep(0, times = nrow(pop_init_df))
dr_adj<-rep(0, times = nrow(pop_init_df))

for (n in 1:nrow(pop_init_df)){
  tb_temp <-pop_init_df$TB_compartment[n]
  hiv_temp<-pop_init_df$HIV_compartment[n]
  dr_temp<-pop_init_df$DR_compartment[n]
  g_temp<-pop_init_df$G_compartment[n]
  
  pop_init_gen_df_gender<-pop_init_gen_df%>%
    filter(g_compartment == g_temp)
  
  #adjust TB 
  #no one is initialized in compartments 7 or 8
  if (tb_temp < 7){
    tb_compartment_temp<-paste0("TB_", tb_temp)
    tb_adj_row<-which(grepl(tb_compartment_temp, pop_init_gen_df_gender$Compartment))
    tb_adj[n]<-pop_init_gen_df_gender$prop_pop[tb_adj_row]
  }
  
  #adj_hiv
  #there is a different hiv distribution 
  #for people with activeTB
  if (tb_temp == 6){
    hiv_compartment_temp<-paste0("HIV_", hiv_temp, '_activeTB')
    hiv_adj_row<-which(grepl(hiv_compartment_temp, pop_init_gen_df_gender$Compartment))
    hiv_adj[n]<-pop_init_gen_df_gender$prop_pop[hiv_adj_row]
  }else{
    hiv_compartment_temp<-paste0("HIV_", hiv_temp, '_noTB')
    hiv_adj_row<-which(grepl(hiv_compartment_temp, pop_init_gen_df_gender$Compartment))
    hiv_adj[n]<-pop_init_gen_df_gender$prop_pop[hiv_adj_row]
  }
  
  #if tb unifected than all people get assigned to DR1
  if(tb_temp <3){
    if(dr_temp==1){
      dr_adj[n]<-1
    }else{
      dr_adj[n]<-0
    }
  }else{
    if(dr_temp==1){
      dr_adj[n]<-.963
    }else{
      dr_adj[n]<-0.037
    }
  }
  
  
}

pop_init_df$tb_adj <-tb_adj
pop_init_df$hiv_adj<-hiv_adj
pop_init_df$dr_adj<-dr_adj
pop_init_df$total_adj<-pop_init_df$tb_adj*
  pop_init_df$hiv_adj*
  pop_init_df$dr_adj

pop_init_df$total_in_compartment<-pop_init_df$total_adj/sum(pop_init_df$total_adj)*
  100000

setwd(outdir)
write.csv(pop_init_df%>%select(c('TB_compartment', 'HIV_compartment',
                                 'DR_compartment', 'G_compartment', 
                                 'total_in_compartment')), 'pop_init_df.csv',
          row.names=FALSE)


