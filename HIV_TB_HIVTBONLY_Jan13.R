#model calibration Jan 5 2021
#TB only

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 'varhandle', 'here', 'readr'), require, character.only=T)

indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/model_outputs')
setwd(indir)

param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
pop_init_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'pop_init')

names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))
names(pop_init_df)<-str_replace_all(names(pop_init_df), c(" " = "_" , "-" = "_" ))

#make sure all compartments are integer type for proper indexing#
param_df$TB_compartment<-as.integer(param_df$TB_compartment)
param_df$DR_compartment<-as.integer(param_df$DR_compartment)
param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
param_df$G_compartment<-as.integer(param_df$G_compartment)
param_df$P_compartment<-as.integer(param_df$P_compartment)

pop_init_df$TB_compartment<-as.integer(pop_init_df$TB_compartment)
pop_init_df$DR_compartment<-as.integer(pop_init_df$DR_compartment)
pop_init_df$HIV_compartment<-as.integer(pop_init_df$HIV_compartment)

pop_init_df$G_compartment<-as.integer(pop_init_df$G_compartment)

#Chelsea - can you add a note about the purpose of this data frame?
pop_init_df_TBHIV_temp <- pop_init_df%>%
  group_by(TB_compartment, HIV_compartment)%>%
  summarise(value = sum(initialized_population_in_compartment))%>%
  mutate(compartment_id = paste0("N_", TB_compartment, "_", HIV_compartment),
         dcompartment_id = paste0("dN_", TB_compartment, "_", HIV_compartment))

################ DEFINE SETS ###############################

#TB states (TB)
#1:Uninfected, not on IPT;
#2:Uninfected, on IPT; 
#3:LTBI, infected recently (within the past two-years)
#4: LTBI, infected remotely (more than two-years ago)
#5: LTBI, on IPT
#6: Active
#7: Recovered/Treated
#8: LTBI, after IPT

TB_SET<-1:8

#4 HIV compartments (HIV)#
#1 : HIV Negative
#2 : HIV Positive CD4 > 200 - No ART
#3 : HIV Positive CD4 =<: 200 - No Art
#4 : HIV Positive - ART 

HIV_SET<-1:4

#param extraction

#beta - Is this taking the median of beta_1 and beta_2 while our model doesn't have gender?
beta <- param_df%>%
  filter(notation == 'beta')%>%
  group_by(TB_compartment)%>%
  summarise(value = median(Reference_expected_value))

beta <- beta$value
#beta<-1

#phi - relative transmissibility of TB in PLWH
phi_h <- array(0, dim = length(HIV_SET))

lapply(HIV_SET, function(h){
  temp <- param_df%>%
    filter(notation == 'phi',
           HIV_compartment == h)
  
  phi_h[h] <<- temp$Reference_expected_value
})

#epsilon - MDR-TB fraction

#iota - indicator for whether infection with strain can occur while on IPT
iota <- param_df%>%
  filter(notation == 'iota', DR_compartment == 1)
iota <- iota$Reference_expected_value

#upsilon - partially protective effect of IPT after completing IPT course
upsilon <- param_df%>%
  filter(notation == 'upsilon')
upsilon <-upsilon$Reference_expected_value

#zeta <- .7 - partially protective effect of LTBI on acquiring new infection
zeta <- param_df%>%
  filter(notation == 'zeta')
zeta <-zeta$Reference_expected_value


#kappas - rates of IPT initiation by group
kappa_t_h <- array(data = 0, c(length(TB_SET), length(HIV_SET)))



lapply(TB_SET, function(t){
  lapply(HIV_SET, function(h){
    temp <- param_df%>%
      filter(P_compartment == 1,
             notation == 'kappa')%>%
      group_by(TB_compartment, HIV_compartment)%>%
      summarise(value = median(Reference_expected_value))%>%
      filter(TB_compartment == t,
             HIV_compartment == h)
    
    if (nrow(temp) == 1){
      kappa_t_h[t,h] <<- temp$value
    }
  })
})

#varpi - IPT adherence
varpi <- param_df%>%
  filter(P_compartment == 1, notation == 'varpi')%>%
  summarise(value = median(Reference_expected_value))

varpi <- varpi$value

#omega - rate of moving off IPT per year = 2
omega <-param_df%>%filter(notation == 'omega')
omega <- omega$Reference_expected_value

#pi - rates of progression through TB compartments
pi_t_t <- array(data = 0, c(length(TB_SET), length(TB_SET)))

lapply(TB_SET, function(t_from){
  lapply(TB_SET, function(t_to){
    num_temp = (t_from*10) + t_to
    
    temp <- param_df%>%
      filter(notation == 'pi',
             TB_compartment == num_temp)
    
    if (nrow(temp) == 1){
      pi_t_t[t_from,t_to] <<- temp$Reference_expected_value
    }
    
  })
})

#test impact of pi_t_t
pi_t_t[8,6] <- .02
#pi_t_t <- pi_t_t*1.5
#pi_t_t[6,7] <- .5
#pi_t_t[4,6]<-.02
#Jen question - does this show our new arrow from TB compartment 8 to 6?

#theta - HIV impact on relative risk of TB progression
theta_h <-array(0, dim = length(HIV_SET))

lapply(HIV_SET, function(h){
  temp <- param_df%>%
    filter(notation == 'theta',
           HIV_compartment == h)
  
  theta_h[h]<<- temp$Reference_expected_value
  
})

#gamma - indicator so that populations with MDR-TB can't move into LTBI after IPT

eta_i_h <- array(0, dim=c(length(HIV_SET), length(HIV_SET)))

#eta - movement between HIV compartments
lapply(HIV_SET, function(h_from){
  lapply(HIV_SET, function(h_to){
    num_temp = (h_from*10) + h_to
    
    temp <- param_df%>%
      filter(notation == 'eta',
             HIV_compartment == num_temp,
             P_compartment == 1,
             G_compartment == 1)
    
    
    
    if (nrow(temp) == 1){
      eta_i_h[h_from, h_to] <<-temp$Reference_expected_value 
    }
    
  })
})

#mu - mortality rates by HIV and TB status
mu_t_h <- array(0, dim = c(length(TB_SET), length(HIV_SET)))

lapply(TB_SET, function(t){
  lapply(HIV_SET, function(h){
    temp <- param_df%>%
      filter(notation == 'mu')%>%
      group_by(TB_compartment, HIV_compartment)%>%
      summarise(value = median(Reference_expected_value))%>%
      filter(TB_compartment == t,
             HIV_compartment == h)
    
    mu_t_h[t,h] <<- temp$value
    
  })
})

#test impacts of mu
#mu_t <- rep(.01, times = length(TB_SET))

#alpha - aging in and aging out
alpha_in_t_h <- array(data = 0, c(length(TB_SET), length(HIV_SET)))

lapply(TB_SET, function(t){
  lapply(HIV_SET, function(h){
    temp <- param_df%>%
      filter(notation == 'alpha^in')%>%
      group_by(notation, TB_compartment, HIV_compartment)%>%
      summarise(value = sum(Reference_expected_value))%>%
      filter(TB_compartment == t, HIV_compartment == h)
    
    if (nrow(temp) == 1){
      alpha_in_t_h[t,h] <<- temp$value
    }
  })
})

#alpha_in_t[1] <- 0
#alpha_in_t[3] <- .8
#alpha_in_t[4] <- 0.1923905

alpha_out <- param_df%>%filter(notation == 'alpha^out')
alpha_out <- alpha_out$Reference_expected_value

total_out_t_h <- array(0, dim = length(TB_SET)*length(HIV_SET))
count_temp <- 1

#calculate totals
lapply(TB_SET, function(t){
  lapply(HIV_SET, function(h){
    total_out_t_h[count_temp] <<- (mu_t_h[t,h]*(1-alpha_out))+
      ((1-mu_t_h[t,h])*alpha_out)+
      (mu_t_h[t,h]*alpha_out)
    count_temp <<- count_temp + 1
  })
})

N_init <- pop_init_df_TBHIV_temp$value
names(N_init) <- c(pop_init_df_TBHIV_temp$compartment_id)

N_t_h_ref <- array(0, dim = c(length(TB_SET), length(HIV_SET)))
count_temp <- 1

lapply(TB_SET, function(t){
  lapply(HIV_SET, function(h){
    N_t_h_ref[t,h] <<- count_temp
    count_temp <<- count_temp + 1
  })
})

array(0, dim = c(length(TB_SET), length(HIV_SET), 2))

open_seir_model <- function(time, N_t_h, parms){
  
  dN_t_h <- array(0, dim = length(TB_SET)*length(HIV_SET))
  names(dN_t_h) <- pop_init_df_TBHIV_temp$dcompartment_id
  
  FOI <- (beta*(sum((phi_h)*N_t_h[N_t_h_ref[6, HIV_SET]])/sum(N_t_h)))
  
  B <- sum(total_out_t_h*N_t_h)
  
  #TB compartment 1
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[1,h]]<<-((alpha_in_t_h[1,h]*B) + #entries from births
                                (omega*N_t_h[N_t_h_ref[2,h]]) - #entries from off IPT 
                                (total_out_t_h[N_t_h_ref[1,h]]*N_t_h[N_t_h_ref[1,h]]) - #exists from aging out and death
                                (FOI*N_t_h[N_t_h_ref[1,h]])- #exists from TB infection 
                                (kappa_t_h[1,h]*N_t_h[N_t_h_ref[1,h]])) + #exists from on to IPT 
      (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[1,HIV_SET]])) - #entries into HIV compartment
      (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[1,h]]) #exit from HIV compartment
  })
  
  #TB compartment 2
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[2,h]]<<-((kappa_t_h[1,h]*N_t_h[N_t_h_ref[1,h]])- #entries from on to IPT 
                                (total_out_t_h[N_t_h_ref[2,h]]*N_t_h[N_t_h_ref[2,h]])- #exists from aging out and death
                                (omega*N_t_h[N_t_h_ref[2,h]])- #exits from off IPT 
                                (iota*FOI*N_t_h[N_t_h_ref[2,h]]))# +#exits from infection (diminished for IPT)
    (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[2,HIV_SET]])) - #entries into HIV compartment
      (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[2,h]]) #exit from HIV compartment
  })
  
  #TB compartment 3
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[3, h]]<<-((alpha_in_t_h[3,h]*B) + #entries from births
                                 (FOI*N_t_h[N_t_h_ref[1,h]])+ #infection from compartment 1 
                                 (iota*FOI*N_t_h[N_t_h_ref[2,h]])+ #infections from compartment 2 
                                 (zeta*FOI*N_t_h[N_t_h_ref[4,h]])+ #re-infection from compartment 4
                                 (zeta*FOI*N_t_h[N_t_h_ref[7,h]])+ #re-infection from compartment 7
                                 (upsilon*FOI*N_t_h[N_t_h_ref[8,h]])- #re-infection from compartment 8
                                 (total_out_t_h[N_t_h_ref[3,h]]*N_t_h[N_t_h_ref[3,h]]) - #exists from aging out and death
                                 (pi_t_t[3,4]*N_t_h[N_t_h_ref[3,h]]) - #from recent to remote infection
                                 (kappa_t_h[3,h]*N_t_h[N_t_h_ref[3,h]]) - #from recent to on IPT 
                                 ((1/varpi)*theta_h[h]*pi_t_t[3,6]*N_t_h[N_t_h_ref[3,h]]) + #from recent TB infection to active
                                 (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[3,HIV_SET]])) - #entries into HIV compartment
                                 (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[3,h]]) #exit from HIV compartment
    )
  })
  
  #TB compartment 4
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[4,h]]<<-((alpha_in_t_h[4,h]*B) + #entries from births
                                (pi_t_t[3,4]*N_t_h[N_t_h_ref[3,h]]) - #from recent to remote infection
                                (total_out_t_h[N_t_h_ref[4,h]]*N_t_h[N_t_h_ref[4,h]]) - #exists from aging out and death
                                (zeta*FOI*N_t_h[N_t_h_ref[4,h]])- #re-infection
                                (kappa_t_h[4,h]*N_t_h[N_t_h_ref[4,h]]) - #onto IPT
                                ((1/varpi)*theta_h[h]*pi_t_t[4,6]*N_t_h[N_t_h_ref[4,h]]) + #from remote TB infection to active 
                                (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[4,HIV_SET]])) - #entries into HIV compartment
                                (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[4,h]]) #exit from HIV compartment
    )
  })
  
  #TB compartment 5
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[5, h]] <<-((kappa_t_h[3,h]*N_t_h[N_t_h_ref[3,h]]) + #from recent to on IPT  
                                  (kappa_t_h[4,h]*N_t_h[N_t_h_ref[4,h]]) - #onto IPT
                                  (total_out_t_h[N_t_h_ref[5, h]]*N_t_h[N_t_h_ref[5, h]]) - #exists from aging out and death
                                  ((1/varpi)*theta_h[h]*pi_t_t[5,6]*N_t_h[N_t_h_ref[5, h]]) - #from TB infection on IPT to active
                                  (omega*N_t_h[N_t_h_ref[5, h]]) + #off IPT to after IPT
                                  (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[5,HIV_SET]])) - #entries into HIV compartment
                                  (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[5,h]]) #exit from HIV compartment
    )
  })
  
  #TB compartment 6
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[6, h]] <<- ((alpha_in_t_h[6,h]*B) + #entries from births
                                   ((1/varpi)*theta_h[h]*pi_t_t[3,6]*N_t_h[N_t_h_ref[3,h]]) + #from recent TB infection to active 
                                   ((1/varpi)*theta_h[h]*pi_t_t[4,6]*N_t_h[N_t_h_ref[4,h]]) + #from remote TB infection to active 
                                   ((1/varpi)*theta_h[h]*pi_t_t[5,6]*N_t_h[N_t_h_ref[5,h]]) + #from TB infection on IPT to active
                                   ((1/varpi)*theta_h[h]*pi_t_t[8,6]*N_t_h[N_t_h_ref[8,h]]) - #from TB after on IPT to active
                                   (total_out_t_h[N_t_h_ref[6,h]]*N_t_h[N_t_h_ref[6, h]]) - #total out
                                   (pi_t_t[6,7]*N_t_h[N_t_h_ref[6, h]]) +#from active to recovered
                                   (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[6,HIV_SET]])) - #entries into HIV compartment
                                   (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[6,h]]) #exit from HIV compartment
    )
    
  })
  
  #TB compartment 7
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[7, h]] <<- ((pi_t_t[6,7]*N_t_h[N_t_h_ref[6, h]]) - #from active to recovered
                                   (total_out_t_h[N_t_h_ref[7, h]]*N_t_h[N_t_h_ref[7, h]]) - #total out
                                   (zeta*FOI*N_t_h[N_t_h_ref[7,h]]) + #re-infection from compartment 7
                                   (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[7,HIV_SET]])) - #entries into HIV compartment
                                   (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[7,h]]) #exit from HIV compartment
    )
  })
  
  #TB compartment 8
  lapply(HIV_SET, function(h){
    dN_t_h[N_t_h_ref[8,h]] <<- ((omega*N_t_h[N_t_h_ref[5,h]]) - #off IPT to after IPT
                                  ((1/varpi)*theta_h[h]*pi_t_t[8,6]*N_t_h[N_t_h_ref[8,h]]) - #from TB after on IPT to active
                                  (total_out_t_h[N_t_h_ref[8,h]]*N_t_h[N_t_h_ref[8,h]]) - #total out
                                  (upsilon*FOI*N_t_h[N_t_h_ref[8,h]])+ 
                                  (sum(eta_i_h[HIV_SET, h]*N_t_h[N_t_h_ref[8,HIV_SET]])) - #entries into HIV compartment
                                  (sum(eta_i_h[h,HIV_SET])*N_t_h[N_t_h_ref[8,h]]) #exit from HIV compartment
    )
  })
  
  list(dN_t_h)
}

#Time Horizon 
TT<-5 #2017-1990
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval)

out<-as.data.frame(ode(times = TT_SET, y = N_init, 
                       func = open_seir_model, method = 'lsoda',
                       parms = NULL))

#out_melt<-melt(data = out, 
#              id.vars = c("time"))

#colnames(out_melt)[2] <- 'TB_compartment'

#ggplot(out_melt, aes(x = time, y = value, group = TB_compartment, color = TB_compartment))+
#  geom_line()+
#  ylab('total in compartment')

out <-out%>%
  mutate(total_pop = rowSums(.[2:ncol(out)]))

#setwd(outdir)
#write_csv(out, 'out.csv')

out_melt = melt(out, id.vars = c('time', 'total_pop'))
temp <- strsplit(as.character(out_melt$variable), "_")

TB_compartment_temp <- rep(0, times = length(temp))
HIV_compartment_temp <- rep(0, times = length(temp))

for (n in 1:length(temp)){
  listt <- temp[[n]]
  TB_compartment_temp[n] <-listt[2]
  HIV_compartment_temp[n] <- listt[3]
}

out_melt$TB_compartment <- TB_compartment_temp
out_melt$HIV_compartment <- HIV_compartment_temp

out_melt_grouped <- out_melt%>%
  #filter(TB_compartment == '6')%>%
  group_by(TB_compartment, time)%>%
  summarise(total_pop = sum(value))

ggplot(out_melt_grouped, aes(x = time, y = total_pop, group = TB_compartment, color = TB_compartment))+
  geom_line()+
  ylab('total in compartment')

ggplot(out_melt_grouped%>%filter(TB_compartment=='6'), aes(x = time, y = total_pop))+
  geom_line()+
  ylab('total with Active TB')+
  ylim(0, 30000)