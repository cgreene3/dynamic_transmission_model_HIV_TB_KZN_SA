#Pop init

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) 
#to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/param_files')

#define compartment sets
TB_compartments <- 1:8
HIV_compartments<-1:4
DR_compartments<-1:2
G_compartments<-1:2

#to name plots according to date generated
data_gen_date<-Sys.Date()
data_gen_date<-str_replace_all(data_gen_date, '-', '_')

#TB factors 
prop_recent<-0.049  #proprotion of LTBI that is recent
prop_recovered<-.22
prop_remote<-1-prop_recent-prop_recovered #proportion remote
prop_IPT<-.01
prop_LTBI<-.49
LTBI_recent<-(prop_LTBI*prop_recent)*(1-prop_IPT)
LTBI_remote<-(prop_LTBI*prop_remote)*(1-prop_IPT)
LTBI_recovered<-(prop_LTBI*prop_recovered)*(1-prop_IPT)
LTBI_on_IPT<-(prop_LTBI*prop_IPT)
TB_active_male<-.0045
TB_active_female<-.013
TB_uninfected_male<-1-.49-TB_active_male
TB_uninfected_female<-1-.49-TB_active_female
prop_DS<-1-0.037
prop_MDR<-0.037
#HIV prev from 1990 (cara's numbers) - see hiv input gen data
prop_hiv_male<-c(0.94132, 0.05368,0.00500, 0.00000)
prop_hiv_female<-c(0.90057, 0.09300,0.00643, 0.00000)

TB_compartment<-c()
G_compartment<-c()
tb_adj<-c()

for (tb in TB_compartments){
  for (g in G_compartments){
    TB_compartment<-c(TB_compartment, tb)
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
      } else if (tb == 7){
        tb_adj<-c(tb_adj, LTBI_recovered)
      } else {
        tb_adj<-c(tb_adj, 0) #no one is birthed to tb compartment 8
      }
  }
}


pop_init_df<-data.frame(TB_compartment, G_compartment, 
                        tb_adj)
pop_init_df$tb_adj<-pop_init_df$tb_adj/2 #assuming 50-50 males/females
pop_init_df$TB_total_pop<-pop_init_df$tb_adj*100000

rm(tb_adj,TB_compartment, G_compartment)


#dr adj
pop_init_df2<-rbind(pop_init_df, pop_init_df)
pop_init_df2$DR_compartment = c(rep(1, times = nrow(pop_init_df)),
                                 rep(2, times = nrow(pop_init_df)))

dr_adj<-c()

for (n in 1:nrow(pop_init_df2)){
  tb = pop_init_df2$TB_compartment[n]
  dr = pop_init_df2$DR_compartment[n]
  
  #assign all uninfected to DR 1
  if(tb == 1 || tb == 2){
    if (dr == 1){
      dr_adj<-c(dr_adj,1)
    } else {
      dr_adj<-c(dr_adj,0)
    }
  } else {
    if (dr == 1){
      dr_adj<-c(dr_adj, prop_DS)
    } else {
      dr_adj<-c(dr_adj, prop_MDR)
    }
  }
}

pop_init_df2$dr_adj<-dr_adj

#hiv adj
pop_init_df3<-rbind(pop_init_df2, pop_init_df2, pop_init_df2, pop_init_df2)
pop_init_df3$HIV_compartment = c(rep(1, times = nrow(pop_init_df2)),
                                rep(2, times = nrow(pop_init_df2)),
                                rep(3, times = nrow(pop_init_df2)),
                                rep(4, times = nrow(pop_init_df2)))

hiv_adj<-c()
for (n in 1:nrow(pop_init_df3)){
  hiv = pop_init_df3$HIV_compartment[n]
  gender = pop_init_df3$G_compartment[n]
  
  if(gender == 1){
    hiv_adj<-c(hiv_adj, prop_hiv_male[hiv])
  } else {
    hiv_adj <-c(hiv_adj, prop_hiv_female[hiv])
  }
}

pop_init_df3$hiv_adj<-hiv_adj

#adjust all
pop_init_df<-pop_init_df3%>%
  mutate(pop_in_compartment = TB_total_pop*dr_adj*hiv_adj)%>%
  select(c('TB_compartment', 'DR_compartment', 
           'HIV_compartment', 'G_compartment', 'pop_in_compartment'))

outdir <- paste0(here(),'/param_files')
setwd(outdir)

write.csv(pop_init_df, 'pop_init_df.csv', row.names = FALSE)

test<-pop_init_df%>%
  group_by(TB_compartment, G_compartment)%>%
  summarise(total = sum(pop_in_compartment))%>%
  ungroup()%>%
  group_by(G_compartment)%>%
  mutate(perc = total/sum(total))

