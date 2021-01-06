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

pop_init_df_TB_temp <- pop_init_df%>%
  group_by(TB_compartment)%>%
  summarise(value = sum(initialized_population_in_compartment))%>%
  mutate(compartment_id = paste0("N_", TB_compartment),
         dcompartment_id = paste0("dN_", TB_compartment))

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

#param extraction

#beta
beta <- param_df%>%
  filter(notation == 'beta')%>%
  group_by(TB_compartment)%>%
  summarise(value = median(Reference_expected_value))

beta <- beta$value
#beta<-1

#phi

#epsilon

#iota
iota <- param_df%>%
  filter(notation == 'iota', DR_compartment == 1)
iota <- iota$Reference_expected_value

upsilon <- param_df%>%
  filter(notation == 'upsilon')
upsilon <-upsilon$Reference_expected_value

zeta <- param_df%>%
  filter(notation == 'zeta')
zeta <-zeta$Reference_expected_value

#zeta <- .7

kappa_t <- array(data = 0, length(TB_SET))

lapply(TB_SET, function(t){
  temp <- param_df%>%
    filter(P_compartment == 1,
           notation == 'kappa')%>%
    group_by(TB_compartment)%>%
    summarise(value = median(Reference_expected_value))%>%
    filter(TB_compartment == t)
  
  if (nrow(temp) == 1){
    kappa_t[t] <<- temp$value
  }
})

varpi <- param_df%>%
  filter(P_compartment == 1, notation == 'varpi')%>%
  summarise(value = median(Reference_expected_value))

varpi <- varpi$value
#test
#varpi <- .02

omega <-param_df%>%filter(notation == 'omega')
omega <- omega$Reference_expected_value

omega <- 1

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

#phi

#gamma

#eta

mu_t <- param_df%>%
  filter(notation == 'mu')%>%
  group_by(TB_compartment)%>%
  summarise(value = median(Reference_expected_value))
mu_t <- mu_t$value

#test impacts of mu
#mu_t <- rep(.01, times = length(TB_SET))

alpha_in_t <- array(data = 0, length(TB_SET))

lapply(TB_SET, function(t){
  temp <- param_df%>%
    filter(notation == 'alpha^in')%>%
    group_by(notation, TB_compartment)%>%
    summarise(value = sum(Reference_expected_value))%>%
    filter(TB_compartment == t)
  
  if (nrow(temp) == 1){
    alpha_in_t[t] <<- temp$value
  }
})

#alpha_in_t[1] <- 0
#alpha_in_t[3] <- .8
#alpha_in_t[4] <- 0.1923905

alpha_out <- param_df%>%filter(notation == 'alpha^out')
alpha_out <- alpha_out$Reference_expected_value

total_out_t <- (mu_t*(1-alpha_out))+((1-mu_t)*alpha_out)+(mu_t*alpha_out)

N_init <- pop_init_df_TB_temp$value
names(N_init) <- c(pop_init_df_TB_temp$compartment_id)


open_seir_model <- function(t, N_t, parms){
  
  dN_t <- array(0, dim = length(TB_SET))
  names(dN_t) <- pop_init_df_TB_temp$dcompartment_id
  
  FOI <- beta*(N_t[6]/sum(N_t))
  
  B <- sum(total_out_t*N_t)
  
  dN_t[1]<-((alpha_in_t[1]*B) + #entries from births
            (omega*N_t[2]) - #entries from off IPT 
              (total_out_t[1]*N_t[1]) - #exists from aging out and death
              (FOI*N_t[1])- #exists from TB infection 
              (kappa_t[1]*N_t[1])) #exists from on to IPT 
  
  dN_t[2]<-((kappa_t[1]*N_t[1])- #entries from on to IPT 
              (total_out_t[2]*N_t[2])- #exists from aging out and death
              (omega*N_t[2])- #exits from off IPT 
              (iota*FOI*N_t[2])) #exits from infection (diminished for IPT) 
  
  dN_t[3]<-((alpha_in_t[3]*B) + #entries from births
              (FOI*N_t[1])+ #infection from compartment 1 
              (iota*FOI*N_t[2])+ #infections from compartment 2 
              (zeta*FOI*N_t[4])+ #re-infection from compartment 4
              (zeta*FOI*N_t[7])+ #re-infection from compartment 7
              (upsilon*FOI*N_t[8])- #re-infection from compartment 8
              (total_out_t[3]*N_t[3]) - #exists from aging out and death
              (pi_t_t[3,4]*N_t[3]) - #from recent to remote infection
              (kappa_t[3]*N_t[3]) - #from recent to on IPT 
              ((1/varpi)*pi_t_t[3,6]*N_t[3]) #from recent TB infection to active 
            )
  
  dN_t[4]<-((alpha_in_t[4]*B) + #entries from births
              (pi_t_t[3,4]*N_t[3]) - #from recent to remote infection
              (total_out_t[4]*N_t[4]) - #exists from aging out and death
              (zeta*FOI*N_t[4])- #re-infection
              (kappa_t[4]*N_t[4]) - #onto IPT
              ((1/varpi)*pi_t_t[4,6]*N_t[4]) #from remote TB infection to active 
  )
  
  dN_t[5] <-(kappa_t[3]*N_t[3] + #from recent to on IPT  
               (kappa_t[4]*N_t[4]) - #onto IPT
               (total_out_t[5]*N_t[5]) - #exists from aging out and death
               ((1/varpi)*pi_t_t[5,6]*N_t[5]) - #from TB infection on IPT to active
               (omega*N_t[5])  #off IPT to after IPT
    
  )
  
  dN_t[6] <- ((alpha_in_t[6]*B) + #entries from births
                ((1/varpi)*pi_t_t[3,6]*N_t[3]) + #from recent TB infection to active 
                ((1/varpi)*pi_t_t[4,6]*N_t[4]) + #from remote TB infection to active 
                ((1/varpi)*pi_t_t[5,6]*N_t[5]) + #from TB infection on IPT to active
                ((1/varpi)*pi_t_t[8,6]*N_t[8]) - #from TB after on IPT to active
                (total_out_t[6]*N_t[6]) - #total out
                (pi_t_t[6,7]*N_t[6]) #from active to recovered
              )
  
  dN_t[7] <- ((pi_t_t[6,7]*N_t[6]) - #from active to recovered
                (total_out_t[7]*N_t[7]) - #total out
                (zeta*FOI*N_t[7])) #re-infection from compartment 7
  
  dN_t[8] <- ((omega*N_t[5]) - #off IPT to after IPT
                ((1/varpi)*pi_t_t[8,6]*N_t[8]) - #from TB after on IPT to active
                (total_out_t[8]*N_t[8]) - #total out
                (upsilon*FOI*N_t[8])
              )
  
  list(dN_t)
}

#Time Horizon 
TT<-100 #2017-1990
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval)

out<-as.data.frame(ode(times = TT_SET, y = N_init, 
         func = open_seir_model, method = 'lsoda',
           parms = NULL))

out_melt<-melt(data = out, 
               id.vars = c("time"))

colnames(out_melt)[2] <- 'TB_compartment'

ggplot(out_melt, aes(x = time, y = value, group = TB_compartment, color = TB_compartment))+
  geom_line()+
  ylab('total in compartment')

out <-out%>%
  mutate(total_pop = rowSums(.[2:ncol(out)]))

setwd(outdir)
write_csv(out, 'out.csv')
