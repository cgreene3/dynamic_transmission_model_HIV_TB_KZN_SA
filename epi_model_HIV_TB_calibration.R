#model Feb 22
#TB/HIV/DR/G + mort changing over time
#+ HIV incidence and initiation changing overtime

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 'varhandle', 'here', 'readr'), require, character.only=T)

#set in directory and out directory
#Make sure you have the epi_model_HIV_TB.Rproj open, otherwise 
#you will need to change the working directory manually.
indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/model_outputs')

#read in data
setwd(indir)
param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
pop_init_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'pop_init')
mort_df <- read.csv('mort_df.csv')
hiv_transition_df<-read.csv('hiv_transmission_df.csv')

#clean dataframe column names for consistency
names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))
names(pop_init_df)<-str_replace_all(names(pop_init_df), c(" " = "_" , "-" = "_" ))

#make sure all compartments are integer type for proper indexing#
param_df$TB_compartment<-as.integer(param_df$TB_compartment)
param_df$DR_compartment<-as.integer(param_df$DR_compartment)
param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
param_df$G_compartment<-as.integer(param_df$G_compartment)
param_df$P_compartment<-as.integer(param_df$P_compartment)
mort_df$TB_compartment<-as.integer(mort_df$TB_compartment)
mort_df$HIV_compartment<-as.integer(mort_df$HIV_compartment)
mort_df$G_compartment<-as.integer(mort_df$G_compartment)

mort_df$year<-as.integer(mort_df$year)

#establish compartment and dcompartment ids
pop_init_df<- pop_init_df%>%
  mutate(compartment_id = paste0("N_", TB_compartment, "_", DR_compartment ,"_", HIV_compartment, "_", G_compartment),
         dcompartment_id = paste0("dN_", TB_compartment, "_", DR_compartment ,"_", HIV_compartment, "_", G_compartment))

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

######## PARAMETERS THAT IMPACT FORCE OF INFECTION #######

######### beta_g - Number of effective contacts for TB transmission per infectious year######
beta_g <- param_df%>%
  filter(notation == 'beta')%>%
  arrange(G_compartment)

beta_g <- beta_g$Reference_expected_value
#beta<-1

###### phi_h - Relative transmissibility of TB in HIV pops#########
phi_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'phi',
           HIV_compartment == h)
  phi_h[h] <- temp$Reference_expected_value
}

#### varepsilon_g - Fraction of new TB infections that are MDR-TB ####
varepsilon_g <- param_df%>%
  filter(notation == 'varepsilon')%>%
  arrange(G_compartment)

varepsilon_g <- varepsilon_g$Reference_expected_value

##### iota_r - Indicator for whether infection with given TB strain can occur while on IPT by DR compartment#####
iota_r <- param_df%>%
  filter(notation == 'iota')
iota_r <- iota_r$Reference_expected_value

##### zeta - Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection ######
zeta <- param_df%>%
  filter(notation == 'zeta')
zeta <-zeta$Reference_expected_value

#########Parameters that Describe TB progression ######

#### kappa_t_h_g_p - Rate of IPT initiation, per year ####
#currently averaged over gender, and only testing for policy 1
kappa_t_h_g <- array(data = 0, c(length(TB_SET), length(HIV_SET), length(G_SET)))

for (t in TB_SET){
  for (h in HIV_SET){
    for (g in G_SET){
      temp <- param_df%>%
        filter(P_compartment == 1,
               notation == 'kappa',
               TB_compartment == t,
               HIV_compartment == h,
               G_compartment == g)
      
      if (nrow(temp) == 1){
        kappa_t_h_g[t,h,g] <- temp$Reference_expected_value
      }
    }
  }
}

####### varpi_g - IPT adherence #######
#currently averaged over gender, and only testing for policy 1
varpi_g <- param_df%>%
  filter(P_compartment == 1, notation == 'varpi')%>%
  arrange(G_compartment)

varpi_g <- varpi_g$Reference_expected_value

##### omega - Rate of moving off of IPT, per year ####
omega <-param_df%>%filter(notation == 'omega')
omega <- omega$Reference_expected_value


######### pi_i_t - Base rates of TB progression  #####
pi_i_t <- array(data = 0, c(length(TB_SET), length(TB_SET)))

for (t_from in TB_SET){
  for (t_to in TB_SET){
    num_temp = (t_from*10) + t_to
    
    temp <- param_df%>%
      filter(notation == 'pi',
             TB_compartment == num_temp)
    
    if (nrow(temp) == 1){
      pi_i_t[t_from,t_to] <- temp$Reference_expected_value
    }
  }
}

#test impact of pi_i_t
pi_i_t[8,6] <- .02

#########theta_h - relative risk of TB progression###########
theta_h <-array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'theta',
           HIV_compartment == h)
  
  theta_h[h]<- temp$Reference_expected_value
}

###########gamma_r -indicator if DR compartment can move onto after IPT####
gamma_r <- c(1,0)

#######Parameters that describe HIV progression########
#####eta_i_h_g rate HIV transitions ######
#(currently averaged over gender)#

HIV_transitions_param_func<-function(yr){
  eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                              length(HIV_SET), 
                              length(G_SET)))
  for (g in G_SET){
    #get HIV progression CD4 counts
    temp <- param_df%>%
      filter(notation == 'eta',
             HIV_compartment == 23,
             P_compartment == 1,
             G_compartment == g)
    
    eta_i_h_g[2,3,g]<-temp$Reference_expected_value
    
    gender_name <-if_else(g == 1, 'Males', 'Females')
    
    temp2<-hiv_transition_df%>%
      filter(gender == gender_name,
             year == yr)
    
    eta_i_h_g[1,2,g]<-temp2$hiv_incidence
    eta_i_h_g[2,4,g]<-temp2$eta_24
    eta_i_h_g[3,4,g]<-temp2$eta_34
  }
  return(eta_i_h_g)
}
      
#########parameters for death and aging rates ###########

######mu_t_h_g - mortality rates ########

mort_param_func <-function(yr){
  mu_t_h_g <- array(0, dim = c(length(TB_SET), length(HIV_SET), length(G_SET)))
  for (t in TB_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- mort_df%>%
          filter(year == yr,
                 TB_compartment == t,
                 HIV_compartment == h,
                 G_compartment == g)
        
        mu_t_h_g[t,h,g] <- temp$total_mort
      }
    }
  }
  return(mu_t_h_g)
}

########alpha_in_t_h_g - Proportion of population that enters each compartment#####
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), length(DR_SET), length(HIV_SET), length(G_SET)))

for (t in TB_SET){
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- param_df%>%
          filter(notation == 'alpha^in',
                 TB_compartment == t, 
                 HIV_compartment == h, 
                 DR_compartment == r, 
                 G_compartment == g)
        
        if (nrow(temp) == 1){
          alpha_in_t_r_h_g[t,r,h,g] <- temp$Reference_expected_value
        }
      }
    }
  }
}

####### alpha_out - Rate of exit from the population ######
alpha_out <- param_df%>%filter(notation == 'alpha^out')
alpha_out <- alpha_out$Reference_expected_value


#############Pre-processing parameter equations, for ease of use in ode solver#######

####total_out_t_r_h - total amount leaving from compartment######

total_out_param_func<-function(mu_t_h_g){
  total_out_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET))
  count_temp <- 1
  for (t in TB_SET){
    for (r in DR_SET){
      for (h in HIV_SET){
        for (g in G_SET){
          total_out_t_r_h_g[count_temp] <- (mu_t_h_g[t,h,g]*(1-alpha_out))+
            ((1-mu_t_h_g[t,h,g])*alpha_out)+
            (mu_t_h_g[t,h,g]*alpha_out)
          count_temp <- count_temp + 1
        }
      }
    }
  }
  return(total_out_t_r_h_g)
}

#####N_init - total initial population in each compartment flattened in 1D array for ODE solver #####
#arrange pop init df, TB-->DR-->HIV-->G
pop_init_df<-pop_init_df%>%
  arrange(G_compartment)%>%
  arrange(HIV_compartment)%>%
  arrange(DR_compartment)%>%
  arrange(TB_compartment)


N_init <- pop_init_df$initialized_population_in_compartment
names(N_init) <- c(pop_init_df$compartment_id)

####N_t_r_h_g - matrix that identifies the location of compartment in 1D array#####
N_t_r_h_g_ref <- array(0, dim = c(length(TB_SET), length(DR_SET), length(HIV_SET), length(G_SET)))
count_temp <- 1

for (t in TB_SET){
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        N_t_r_h_g_ref[t,r,h,g] <- count_temp
        count_temp <- count_temp + 1
      }
    }
  }
}

#where the equations are stored
open_seir_model <- function(time, N_t_r_h_g, parms){
  
  dN_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET))
  names(dN_t_r_h_g) <- pop_init_df$dcompartment_id
  
  #calculate time varying parameters
  
  #Force of Infection Calculations
  FOI_1_g <- array(0, dim = length(G_SET))
  FOI_2_g <- array(0, dim = length(G_SET))
  
  for (g in G_SET){
    FOI_1_g[g]<-(beta_g[g]*(sum((phi_h)*N_t_r_h_g[N_t_r_h_g_ref[6, 1, HIV_SET,g]])/sum(N_t_r_h_g)))
  }
  
  for (g in G_SET){
    FOI_2_g[g] <-(varepsilon_g[g]*FOI_1_g[g])/(1-varepsilon_g[g])
  }

  FOI_r <- c(sum(FOI_1_g), sum(FOI_2_g))
  FOI <- sum(FOI_r)
  
  #write year parameter for parameters that change over time
  current_yr <-as.integer(start_yr+time)
  
  print(current_yr)
  
  #HIV transitions
  eta_i_h_g<-HIV_transitions_param_func(current_yr)
  
  #deaths
  mu_t_h_g<-mort_param_func(current_yr)
  total_out_t_r_h_g<-total_out_param_func(mu_t_h_g)
  
  #entries and exits from the population
  B <- sum(total_out_t_r_h_g*N_t_r_h_g)
  
  #######TB compartment 1 Equations #########
  #Set DR compartment to 1, since not applicable to drug resistant compartments#
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]<-(sum(alpha_in_t_r_h_g[1,1,h,g]*B) + #entries from births
                                 (omega*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) - #entries from off IPT 
                                 (total_out_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) - #exists from aging out and death
                                 (FOI*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])- #exists from TB infection 
                                 (kappa_t_h_g[1,h,g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) + #exists from on to IPT 
        (sum(eta_i_h_g[HIV_SET, h, g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,HIV_SET,g]])) - #entries into HIV compartment
        (sum(eta_i_h_g[h,HIV_SET, g])*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) #exit from HIV compartment
      )
    }
  }
  
  #############TB compartment 2 Equations########
  #Set DR compartment to 1, since not applicable to drug resistant compartments#
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]<-((kappa_t_h_g[1,h,g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])- #entries from on to IPT 
                                         (total_out_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])- #exists from aging out and death
                                         (omega*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])- #exits from off IPT 
                                         ((sum(iota_r*FOI_r))*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) +#exits from infection (diminished for IPT)
                                         (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[2,1,HIV_SET,g]])) - #entries into HIV compartment
                                         (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) #exit from HIV compartment
      )
    }
  }
  
  #############TB compartment 3 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]<-((alpha_in_t_r_h_g[3,r,h,g]*B) + #entries from births
                                         (FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])+ #infection from compartment 1 
                                         (iota_r[r]*FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])+ #infections from compartment 2 
                                         (zeta*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,DR_SET,h,g]]))+ #re-infection from compartment 4
                                         (zeta*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[7,DR_SET,h,g]]))+ #re-infection from compartment 7
                                         (zeta*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,DR_SET,h,g]]))- #re-infection from compartment 8
                                         (total_out_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #exists from aging out and death
                                         (pi_i_t[3,4]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #from recent to remote infection
                                         (kappa_t_h_g[3,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #from recent to on IPT 
                                         ((1/varpi_g[g])*theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from recent TB infection to active
                                         (sum(eta_i_h_g[HIV_SET,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,HIV_SET,g]])) - #entries into HIV compartment
                                         (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 4 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]<-((alpha_in_t_r_h_g[4,r,h,g]*B) + #entries from births
                                           (pi_i_t[3,4]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #from recent to remote infection
                                           (total_out_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) - #exists from aging out and death
                                           (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]])- #re-infection
                                           (kappa_t_h_g[4,h,g]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) - #onto IPT
                                           ((1/varpi_g[g])*theta_h[h]*pi_i_t[4,6]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #from remote TB infection to active 
                                           (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[4,r,HIV_SET,g]])) - #entries into HIV compartment
                                           (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 5 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]] <-((kappa_t_h_g[3,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from recent to on IPT  
                                          (kappa_t_h_g[4,h,g]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) - #onto IPT
                                          (total_out_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) - #exists from aging out and death
                                          ((1/varpi_g[g])*theta_h[h]*pi_i_t[5,6]*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) - #from TB infection on IPT to active
                                          (gamma_r[r]*omega*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) + #off IPT to after IPT
                                          (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[5,r,HIV_SET,g]])) - #entries into HIV compartment
                                          (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 6 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]] <- ((alpha_in_t_r_h_g[6,r,h,g]*B) + #entries from births
                                           ((1/varpi_g[g])*theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from recent TB infection to active 
                                           ((1/varpi_g[g])*theta_h[h]*pi_i_t[4,6]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #from remote TB infection to active 
                                           ((1/varpi_g[g])*theta_h[h]*pi_i_t[5,6]*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) + #from TB infection on IPT to active
                                           ((1/varpi_g[g])*theta_h[h]*pi_i_t[8,6]*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) - #from TB after on IPT to active
                                           (total_out_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) - #total out
                                           (pi_i_t[6,7]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) +#from active to recovered
                                           (sum(eta_i_h_g[HIV_SET, h, g]*N_t_r_h_g[N_t_r_h_g_ref[6,r,HIV_SET,g]])) - #entries into HIV compartment
                                           (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 7 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]] <- ((pi_i_t[6,7]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) - #from active to recovered
                                           (total_out_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) - #total out
                                           (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) + #re-infection from compartment 7
                                           (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[7,r,HIV_SET,g]])) - #entries into HIV compartment
                                           (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  #############TB compartment 8 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]] <- ((gamma_r[r]*omega*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) - #off IPT to after IPT
                                           ((1/varpi_g[g])*theta_h[h]*pi_i_t[8,6]*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) - #from TB after on IPT to active
                                           (total_out_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) - #total out
                                           (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]])+ 
                                           (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[8,r,HIV_SET,g]])) - #entries into HIV compartment
                                           (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) #exit from HIV compartment
        )
      }
    }
  }
  
  list(dN_t_r_h_g)
}

#Time Horizon and Evaluation intervals (1 month)
start_yr = 1990
end_yr = 2017
TT<-end_yr-start_yr
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval)

#feed in to solve, evaluation time intervals (TT_SET), 
#initial population (N_init)
#differential equations and dynamic parameter (open_seir_model)
#specify differential solver (lsoda)
# read params from global (NULL)

out_df<-as.data.frame(ode(times = TT_SET, y = N_init, 
                       func = open_seir_model, method = 'lsoda',
                       parms = NULL))

#out_melt<-melt(data = out, 
#              id.vars = c("time"))

out_df <-out_df%>%
  mutate(total_pop = rowSums(.[2:ncol(out_df)]))

out_df<- cbind(year = as.integer(start_yr+out_df$time), 
            beta_1 = rep(beta_g[1], times = nrow(out_df)), beta_2 = rep(beta_g[2], times = nrow(out_df)), 
            out_df)

out_df_melt <-melt(out_df,
                   id.vars = c("time", "year", "beta_1", "beta_2", "total_pop"))

out_df_melt <- cbind(out_df_melt, data.frame(do.call('rbind', strsplit(as.character(out_df_melt$variable),'_',fixed=TRUE))))

names(out_df_melt)[names(out_df_melt) == "X2"] <- "TB_compartment"
names(out_df_melt)[names(out_df_melt) == "X3"] <- "DR_compartment"
names(out_df_melt)[names(out_df_melt) == "X4"] <- "HIV_compartment"
names(out_df_melt)[names(out_df_melt) == "X5"] <- "G_compartment"
out_df_melt<-out_df_melt%>%select(-c('X1'))
out_df_melt$month <- as.integer(((out_df_melt$time%%1)*(12)+1))


#####Calibration Calculations######

mort_est <- array(data = 0, dim = nrow(out_df_melt))

#so I do not need to call on mort_param_func too many times
out_df_melt<-out_df_melt%>%
  arrange(year)

counter <-1
for (yr in start_yr:end_yr){
  temp<-out_df_melt%>%
    filter(year == yr)
  mu_t_h_g<-mort_param_func(yr)
  print(yr)
  for (row in 1:nrow(temp)){
    tb <- as.integer(out_df_melt[row, 'TB_compartment'])
    hiv <- as.integer(out_df_melt[row, 'HIV_compartment'])
    gender <-as.integer(out_df_melt[row, 'G_compartment'])
    pop<-as.double(out_df_melt[row, 'value'])
    mort_rate<-(mu_t_h_g[tb,hiv,gender]*(1/12))
    
    mort_est[counter]<-(mort_rate*pop)
    counter <- counter + 1
  }
}

out_df_melt$mort_est <- mort_est

TB_mort_summary_df<-out_df_melt%>%
  filter(TB_compartment ==6,
         HIV_compartment==1)%>%
  group_by(year)%>%
  summarise(total_mort=sum(mort_est))


#######graphs######

out_melt = melt(out, id.vars = c('time', 'total_pop'))
temp <- strsplit(as.character(out_melt$variable), "_")

TB_compartment_temp <- rep(0, times = length(temp))
DR_compartment_temp <- rep(0, times = length(temp))
HIV_compartment_temp <- rep(0, times = length(temp))
G_compartment_temp <- rep(0, times = length(temp))

for (n in 1:length(temp)){
  compartment_split <- temp[[n]]
  TB_compartment_temp[n] <-compartment_split[2]
  DR_compartment_temp[n] <- compartment_split[3]
  HIV_compartment_temp[n] <- compartment_split[4]
  G_compartment_temp[n] <- compartment_split[5]
}

out_melt$TB_compartment <- TB_compartment_temp
out_melt$DR_compartment <- DR_compartment_temp
out_melt$HIV_compartment <- HIV_compartment_temp
out_melt$G_compartment <- G_compartment_temp

#TB compartment grouped
out_melt_grouped_by_TB_compartment <- out_melt%>%
  group_by(TB_compartment, time)%>%
  summarise(total_pop = sum(value))

ggplot(out_melt_grouped_by_TB_compartment, aes(x = time, y = total_pop, group = TB_compartment, color = TB_compartment))+
  geom_line()+
  ylab('total in compartment')

ggplot(out_melt_grouped_by_TB_compartment%>%filter(TB_compartment=='6'), aes(x = time, y = total_pop))+
  geom_line()+
  ylab('total with Active TB')+
  ylim(0, 2500)

#TB HIV compartment grouped
out_melt_grouped_by_TB_HIV_compartment <- out_melt%>%
  group_by(TB_compartment, HIV_compartment, time)%>%
  summarise(total_pop = sum(value))

#par(mfrow=c(2,2))
HIV_TB_graphs <- ggplot(out_melt_grouped_by_TB_HIV_compartment%>%filter(TB_compartment == 6), 
                        aes(x = time, y = total_pop, 
                            group = TB_compartment, 
                            color = TB_compartment))+
  geom_line()+
  ylab('total in compartment')+
  facet_wrap(~HIV_compartment)+
  ggtitle('TB compartment populations by HIV compartment')
dev.off()
