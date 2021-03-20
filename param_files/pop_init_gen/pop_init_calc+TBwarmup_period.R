#Pop init
#20 year warm up TB and gender only
#then add in DR & MDR factors, and HIV to get pop init

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
prop_remote<-0.951 #proportion remote
prop_IPT<-0#.01
prop_LTBI<-.49
LTBI_recent<-(prop_LTBI*prop_recent)*(1-prop_IPT)
LTBI_remote<-(prop_LTBI*prop_remote)*(1-prop_IPT)
LTBI_on_IPT<-(prop_LTBI*prop_IPT)
TB_active_male<-.0045
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

for (tb in TB_compartments){
  for (dr in DR_compartments){
    for (g in G_compartments){
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
    }
  }
}


pop_init_df<-data.frame(TB_compartment, DR_compartment, G_compartment, 
                        tb_adj, dr_adj)

rm(tb_adj, dr_adj, TB_compartment, DR_compartment, G_compartment)

pop_init_df<-pop_init_df%>%
  mutate(total_adj = round(tb_adj*dr_adj,4),
         prop_of_pop = round(total_adj/sum(total_adj),5))%>%
  mutate(prop_of_pop = prop_of_pop*(1/sum(prop_of_pop)))%>% #fix rounding errors make sure add to 1
  mutate(total_pop = prop_of_pop*100000)%>%
  select(c('TB_compartment', 'DR_compartment', 'G_compartment', 'total_pop'))

#load packages for warmup period
sapply(c('deSolve'), require, character.only=T)

setwd(indir)

param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))

#make sure all compartments are integer type for proper indexing#
param_df$TB_compartment<-as.integer(param_df$TB_compartment)
param_df$DR_compartment<-as.integer(param_df$DR_compartment)
param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
param_df$G_compartment<-as.integer(param_df$G_compartment)
param_df$P_compartment<-as.integer(param_df$P_compartment)

pop_init_df$TB_compartment<-as.integer(pop_init_df$TB_compartment)
pop_init_df$DR_compartment<-as.integer(pop_init_df$DR_compartment)
pop_init_df$G_compartment<-as.integer(pop_init_df$G_compartment)

#group to gender and TB compartment only
pop_init_df_TB_G_temp <- pop_init_df%>%
  group_by(TB_compartment, G_compartment)%>%
  summarise(value = sum(total_pop))%>%
  mutate(compartment_id = paste0("N_", TB_compartment, '_', G_compartment),
         dcompartment_id = paste0("dN_", TB_compartment, '_', G_compartment))

#beta (effective contacts)
beta_g <- param_df%>%
  filter(notation == 'beta')

beta_g<-round(beta_g$Reference_expected_value+.002, digits = 3)

#iota (only considering DS)
iota <- param_df%>%
  filter(notation == 'iota', DR_compartment == 1)
iota <- iota$Reference_expected_value

#diminished FOI
zeta <- param_df%>%
  filter(notation == 'zeta')
zeta <-zeta$Reference_expected_value

#initiation of IPT (mean IPT initiation)
kappa_t_g <- array(data = 0, c(length(TB_compartments),
                             length(G_compartments)))

# lapply(TB_compartments, function(t){
#   lapply(G_compartments, function(g){
#     temp <- param_df%>%
#       filter(P_compartment == 1,
#              notation == 'kappa')%>%
#       group_by(TB_compartment, G_compartment)%>%
#       summarise(value = mean(Reference_expected_value))%>%
#       filter(TB_compartment == t,
#              G_compartment == g)
#     
#     if (nrow(temp) == 1){
#       kappa_t_g[t,g] <<- temp$value
#     }
#   })
# })

#adherence
varpi_g <- param_df%>%
  filter(P_compartment == 1, notation == 'varpi')

varpi_g <- varpi_g$Reference_expected_value

#IPT regimen length
omega <-param_df%>%filter(notation == 'omega')
omega <- omega$Reference_expected_value

pi_t_t <- array(data = 0, c(length(TB_compartments), length(TB_compartments)))


lapply(TB_compartments, function(t_from){
  lapply(TB_compartments, function(t_to){
    num_temp = (t_from*10) + t_to
    
    temp <- param_df%>%
      filter(notation == 'pi',
             TB_compartment == num_temp)
    
    if (nrow(temp) == 1){
      pi_t_t[t_from,t_to] <<- temp$Reference_expected_value
    }
    
  })
})

pi_t_t[3,6]<-sum(0.0866*.8,0.998*.2)
pi_t_t[4,6]<-sum(0.0005*.8,0.100*.2)
pi_t_t[3,4]<-.5
pi_t_t[7,6]<-.05

#mortality rates
mu_t_g<-array(0, dim = c(length(TB_compartments), length(G_compartments)))

mu_t_g_temp <- param_df%>%
  filter(notation == 'mu')%>%
  group_by(TB_compartment, G_compartment)%>%
  summarise(value = mean(Reference_expected_value))

lapply(TB_compartments, function(t){
  lapply(G_compartments, function(g){
    temp<-mu_t_g_temp%>%
      filter(TB_compartment == t)%>%
      filter(G_compartment == g)
    
    mu_t_g[t,g] <<- temp$value
  })
})

#aging in
alpha_in_t_g <- array(data = 0, c(length(TB_compartments), length(G_compartments)))
birth_propotions_df<-pop_init_df_TB_G_temp%>%
  filter(TB_compartment != 2)%>%
  filter(TB_compartment != 5)%>%
  filter(TB_compartment < 7)%>%
  ungroup()%>%
  mutate(value = value/sum(value))

lapply(TB_compartments, function(t){
  lapply(G_compartments, function(g){
    if (t %in% unique(birth_propotions_df$TB_compartment)){
      temp<-birth_propotions_df%>%
        filter(TB_compartment == t)%>%
        filter(G_compartment == g)
      alpha_in_t_g[t,g]<<-temp$value
    }
  })
})

alpha_out <- param_df%>%filter(notation == 'alpha^out')
alpha_out <- alpha_out$Reference_expected_value

#calculate total out for each compartment
total_out_t_g<-array(0, dim = c(length(TB_compartments), length(G_compartments)))

lapply(TB_compartments, function(t){
  lapply(G_compartments, function(g){
    total_out_t_g[t,g]<<-(mu_t_g[t,g]*(1-alpha_out))+((1-mu_t_g[t,g])*alpha_out)+(mu_t_g[t,g]*alpha_out)
  })
})

#sort for consistency
pop_init_df_TB_G_temp<-pop_init_df_TB_G_temp%>%
  arrange(TB_compartment)%>%
  arrange(G_compartment)

N_t_g_ref<-array(0, dim = c(length(TB_compartments), length(G_compartments)))

lapply(1:nrow(pop_init_df_TB_G_temp), function(n){
  t_temp <- pop_init_df_TB_G_temp$TB_compartment[n]
  g_temp <- pop_init_df_TB_G_temp$G_compartment[n]
  N_t_g_ref[t_temp,g_temp]<<-n
})


N_init <- pop_init_df_TB_G_temp$value
names(N_init) <- c(pop_init_df_TB_G_temp$compartment_id)

open_seir_model <- function(t, N_t_g, parms){
  
  dN_t_g <- array(0, dim = length(TB_compartments)*length(G_compartments))
  names(dN_t_g) <- pop_init_df_TB_G_temp$dcompartment_id
  
  B <- sum(total_out_t_g[TB_compartments, 1]*N_t_g[N_t_g_ref[TB_compartments, 1]])+
    sum(total_out_t_g[TB_compartments, 2]*N_t_g[N_t_g_ref[TB_compartments, 2]])
  
  FOI <- (sum(beta_g*N_t_g[N_t_g_ref[6,G_compartments]])/sum(N_t_g))
  
  #Uninfected, not on IPT
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[1,g]]<-((alpha_in_t_g[1,g]*B) + #entries from births
                  (omega*N_t_g[N_t_g_ref[2,g]]) - #entries from off IPT 
                  (total_out_t_g[1,g]*N_t_g[N_t_g_ref[1,g]]) - #exists from aging out and death
                  (FOI*N_t_g[N_t_g_ref[1,g]])- #exists from TB infection 
                  (kappa_t_g[1,g]*N_t_g[N_t_g_ref[1,g]])) #exists from on to IPT
  }
   
  #Uninfected, on IPT
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[2,g]]<-((kappa_t_g[1,g]*N_t_g[N_t_g_ref[1,g]])- #entries from on to IPT 
                (total_out_t_g[2,g]*N_t_g[N_t_g_ref[2,g]])- #exists from aging out and death
                (omega*N_t_g[N_t_g_ref[2,g]])- #exits from off IPT 
                (iota*FOI*N_t_g[N_t_g_ref[2,g]])) #exits from infection (diminished for IPT) 
  }
  
  #LTBI, recent
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[3,g]]<-((alpha_in_t_g[3,g]*B) + #entries from births
                (FOI*N_t_g[N_t_g_ref[1,g]])+ #infection from compartment 1 
                (iota*FOI*N_t_g[N_t_g_ref[2,g]])+ #infections from compartment 2 
                (zeta*FOI*N_t_g[N_t_g_ref[4,g]])+ #re-infection from compartment 4
                (zeta*FOI*N_t_g[N_t_g_ref[7,g]])+ #re-infection from compartment 7
                (zeta*FOI*N_t_g[N_t_g_ref[8,g]])- #re-infection from compartment 8
                (total_out_t_g[3,g]*N_t_g[N_t_g_ref[3,g]]) - #exists from aging out and death
                (pi_t_t[3,4]*N_t_g[N_t_g_ref[3,g]]) - #from recent to remote infection
                (kappa_t_g[3,g]*N_t_g[N_t_g_ref[3,g]]) - #from recent to on IPT 
                (pi_t_t[3,6]*N_t_g[N_t_g_ref[3,g]]) #from recent TB infection to active 
    )
  }
  
  #LTBI, remote
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[4,g]]<-((alpha_in_t_g[4,g]*B) + #entries from births
                    (pi_t_t[3,4]*N_t_g[N_t_g_ref[3,g]]) - #from recent to remote infection
                    (total_out_t_g[4,g]*N_t_g[N_t_g_ref[4,g]]) - #exists from aging out and death
                    (zeta*FOI*N_t_g[N_t_g_ref[4,g]])- #re-infection
                    (kappa_t_g[4,g]*N_t_g[N_t_g_ref[4,g]]) - #onto IPT
                    (pi_t_t[4,6]*N_t_g[N_t_g_ref[4,g]]) #from remote TB infection to active 
    )
  }
  
  #LTBI, on IPT
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[5,g]] <-(kappa_t_g[3,g]*N_t_g[N_t_g_ref[3,g]] + #from recent to on IPT  
                                (kappa_t_g[4,g]*N_t_g[N_t_g_ref[4,g]]) - #onto IPT
                                (total_out_t_g[5,g]*N_t_g[N_t_g_ref[5,g]]) - #exists from aging out and death
                                (pi_t_t[5,6]*N_t_g[N_t_g_ref[5,g]]) - #from TB infection on IPT to active
                                (omega*N_t_g[N_t_g_ref[5,g]])  #off IPT to after IPT
                              
    ) 
  }
  
  #active
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[6,g]] <- ((alpha_in_t_g[6,g]*B) + #entries from births
                  (pi_t_t[3,6]*N_t_g[N_t_g_ref[3,g]]) + #from recent TB infection to active 
                  (pi_t_t[4,6]*N_t_g[N_t_g_ref[4,g]]) + #from remote TB infection to active 
                  (pi_t_t[5,6]*N_t_g[N_t_g_ref[5,g]]) + #from TB infection on IPT to active
                  (pi_t_t[8,6]*N_t_g[N_t_g_ref[8,g]]) - #from TB after on IPT to active
                  (total_out_t_g[6,g]*N_t_g[N_t_g_ref[6,g]]) - #total out
                  (pi_t_t[6,7]*N_t_g[N_t_g_ref[6,g]]) + #from active to recovered
                  (pi_t_t[7,6]*N_t_g[N_t_g_ref[7,g]]) #TB relapse
    )
  }
  
  #recovered
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[7,g]] <- ((pi_t_t[6,7]*N_t_g[N_t_g_ref[6,g]]) - #from active to recovered
                               (total_out_t_g[7,g]*N_t_g[N_t_g_ref[7,g]]) - #total out
                               (zeta*FOI*N_t_g[N_t_g_ref[7,g]]) - #re-infection from compartment 7
                               (pi_t_t[7,6]*N_t_g[N_t_g_ref[7,g]]) #TB relapse
    )
  }
  
  #LTBI, after IPT
  for (g in G_compartments){
    dN_t_g[N_t_g_ref[8,g]] <- ((omega*N_t_g[N_t_g_ref[5,g]]) - #off IPT to after IPT
                  (pi_t_t[8,6]*N_t_g[N_t_g_ref[8,g]]) - #from TB after on IPT to active
                  (total_out_t_g[8,g]*N_t_g[N_t_g_ref[8,g]]) - #total out
                  (zeta*FOI*N_t_g[N_t_g_ref[8,g]])
    )
  }
  
  
  list(dN_t_g)
}

#Time Horizon 
TT<-20 #warm up for 20 years
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval)

out<-as.data.frame(ode(times = TT_SET, y = N_init, 
                       func = open_seir_model, method = 'lsoda',
                       parms = NULL))

out_melt<-melt(data = out, 
               id.vars = c("time"))

out_melt <- cbind(out_melt, 
                     data.frame(do.call('rbind', 
                                        strsplit(as.character(out_melt$variable),
                                                 '_',fixed=TRUE))))

names(out_melt)[names(out_melt) == "X2"] <- "TB_compartment"
names(out_melt)[names(out_melt) == "X3"] <- "G_compartment"

out_melt<-out_melt%>%
  select(-c('X1'))%>%
  mutate(TB_grouping = if_else(TB_compartment == 6, 'Active',
                               if_else(TB_compartment %in% c(1,2), 'Uninfected',
                                       'LTBI')))

#start from end of 20 year warm up period
out_melt_warmed_up<-out_melt%>%
  filter(time == 20)

out_melt_warmed_up_grouped_gender<-out_melt_warmed_up%>%
  group_by(TB_compartment)%>%
  summarise(value = sum(value))

out_melt_grouped<-out_melt%>%
  group_by(time, TB_grouping)%>%
  summarise(value = sum(value))

ggplot(out_melt_grouped%>%filter(TB_grouping!='Active'), aes(x = time, y = value))+
  geom_line(aes(colour = TB_grouping))

ggplot(out_melt%>%filter(TB_compartment == 6), 
       aes(x = time, y = value))+
  geom_line(aes(colour = G_compartment))


rm(list=setdiff(ls(), "out_melt_warmed_up"))


#factors considered
prop_IPT<-.01
prop_DS<-1-0.037
prop_MDR<-0.037

#HIV prev from 1990 (cara's numbers) - see hiv input gen data
prop_hiv_male<-c(0.94132, 0.05368,0.00500, 0.00000)
prop_hiv_female<-c(0.90057, 0.09300,0.00643, 0.00000)

#add in IPT percent
pop_init_df<-out_melt_warmed_up%>%
  select(-c('time'))%>%
  mutate(TB_group = if_else(TB_compartment == 6, 'active',
                            if_else(TB_compartment %in% c(3,4,5,7,8), 'LTBI',
                            'uninfected')))%>%
  mutate(IPT_group = if_else(TB_compartment %in% c(2,5), 'IPT', 
                             if_else(TB_compartment == 6, 'active', 'no IPT')))%>%
  group_by(TB_group, G_compartment)%>%
  mutate(total_in_group = sum(value))%>%
  mutate(value2 = if_else(IPT_group == 'IPT', total_in_group, value))%>%
  mutate(ipt_adj = if_else(IPT_group == 'IPT', prop_IPT, 
                           if_else(IPT_group == 'no IPT', 1-prop_IPT,
                                   1)))%>%
  ungroup()
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
  mutate(pop_in_compartment = value2*ipt_adj*dr_adj*hiv_adj)%>%
  select(c('TB_compartment', 'DR_compartment', 
           'HIV_compartment', 'G_compartment', 'pop_in_compartment'))

outdir <- paste0(here(),'/param_files')
setwd(outdir)

write.csv(pop_init_df, 'pop_init_df.csv', row.names = FALSE)
