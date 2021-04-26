#model March April 1
#TB/HIV/DR/G + mort changing over time
#+ HIV incidence and initiation changing overtime
#calibration calcs to HIV/TB deaths and TB only deaths
#changing birthrates overtime
#pop init corresponds to betas 

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops'), require, character.only=T)

#mort calibration test
#(1 = low, 2 = medium, 3 = high)
#TB_only, TB_HIV_CD4More, TB_HIV_CD4Less, TB_HIV_ART
mort_calib_TB_only_test <- 1:3
mort_calib_TB_HIV_CD4More_test <- 1:3
mort_calib_TB_HIV_CD4Less_test <- 1:3
mort_calib_TB_HIV_ART_test <- 1:3
  
#betas to test
beta_1_test<-c(14)
beta_2_test<-c(14)

#set in directory and out directory
#Make sure you have the epi_model_HIV_TB.Rproj open, otherwise 
#you will need to change the working directory manually.
indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/model_outputs/')

#read in data
setwd(indir)
param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
pop_init_df <- read.csv('pop_init_df.csv')
mort_df <- read.csv('mort_df.csv')
hiv_transition_df<-read.csv('hiv_transmission_df.csv')
birth_rate_df<-read.csv('birth_rate_df.csv')

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

####increased risk of re-infection for recovered populations
xi<-param_df%>%
  filter(notation == 'xi')
xi <-xi$Reference_expected_value

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

#########theta_h - relative risk of TB progression###########
theta_h <-array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'theta',
           HIV_compartment == h)
  
  theta_h[h]<- temp$Reference_expected_value
}

### upsilon_h - increased time infectious due to delayed treatment###
upsilon_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'upsilon',
           HIV_compartment == h)
  upsilon_h[h] <- temp$Reference_expected_value
}

###########gamma_r -indicator if DR compartment can move onto after IPT####
gamma_r <- c(1,0)

#######Parameters that describe HIV progression########
#####eta_i_h_g rate HIV transitions ######

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
    
    if(yr == 2018){
      yr = yr -1
    }
    
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
  
  if(yr == 2018){
    yr = yr -1
  }
  
  for (t in TB_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- mort_df%>%
          filter(year == yr,
                 TB_compartment == t,
                 HIV_compartment == h,
                 G_compartment == g)
        
        mu_t_h_g[t,h,g] <- temp$mort_rate
      }
    }
  }
  return(mu_t_h_g)
}

########alpha_in_t_h_g - Proportion of population that enters each compartment#####
birth_rate_param_func<-function(yr){
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                      length(DR_SET), 
                                      length(HIV_SET), 
                                      length(G_SET)))
if(yr == 2018){
  yr = yr -1
}

for (t in TB_SET){
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- (unique(birth_rate_df$prop_of_pop[(birth_rate_df$year == yr) &
                                                         (birth_rate_df$TB_compartment == t)&
                                                         (birth_rate_df$DR_compartment == r)&
                                                         (birth_rate_df$HIV_compartment == h)&
                                                         (birth_rate_df$G_compartment == g)]))
                                              
        
        alpha_in_t_r_h_g[t,r,h,g] <- temp
        
      }
    }
  }
}
return(alpha_in_t_r_h_g)
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
  
  #initiation delta array, +4 for tracking mort
  dN_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET)+4)
  
  #calculate time varying parameters
  
  #Force of Infection Calculations
  FOI_1_g <- array(0, dim = length(G_SET))
  FOI_2_g <- array(0, dim = length(G_SET))
  
  for (g in G_SET){
    FOI_1_g[g]<-(beta_g[g]*(sum((phi_h)*N_t_r_h_g[N_t_r_h_g_ref[6, 1, HIV_SET,g]])/
                              sum(N_t_r_h_g[1:128])))
  }
  
  for (g in G_SET){
    FOI_2_g[g] <-(varepsilon_g[g]*FOI_1_g[g])/(1-varepsilon_g[g])
  }

  FOI_r <- c(sum(FOI_1_g), sum(FOI_2_g))
  FOI <- sum(FOI_r)
  
  #write year parameter for parameters that change over time
  current_yr <-as.integer(start_yr+time)
  
  print(time)
  
  #HIV transitions
  eta_i_h_g<-HIV_transitions_param_func(current_yr)
  
  #alpha in proportions
  alpha_in_t_r_h_g<-birth_rate_param_func(current_yr)
  
  #deaths
  mu_t_h_g<-mort_param_func(current_yr)
  total_out_t_r_h_g<-total_out_param_func(mu_t_h_g)
  
  #track mort
  dN_t_r_h_g[129]<-mu_t_h_g[6,1,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1,1,1]]) #tb active, hiv -, males
  dN_t_r_h_g[130]<-mu_t_h_g[6,1,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1,1,2]]) #tb active, hiv -, females
  dN_t_r_h_g[131]<-sum(mu_t_h_g[6,2:3,1]*N_t_r_h_g[N_t_r_h_g_ref[6,1,2:3,1]])+ #ds #tb active, hiv +, females
    sum(mu_t_h_g[6,2:3,1]*N_t_r_h_g[N_t_r_h_g_ref[6,2,2:3,1]]) #mdr
  dN_t_r_h_g[132]<-sum(mu_t_h_g[6,2:3,2]*N_t_r_h_g[N_t_r_h_g_ref[6,1,2:3,2]])+ #ds #tb active, hiv +, females
    sum(mu_t_h_g[6,2:3,2]*N_t_r_h_g[N_t_r_h_g_ref[6,2,2:3,2]]) #mdr
  
  #entries and exits from the population
  B <- sum(total_out_t_r_h_g*N_t_r_h_g[1:128])
  
  #######TB compartment 1 Equations #########
  #Set DR compartment to 1, since not applicable to drug resistant compartments#
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]<-(
        (omega*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) - #entries from off IPT
          (FOI*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])- #exists from TB infection 
          (kappa_t_h_g[1,h,g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) + #exists from on to IPT
          sum(alpha_in_t_r_h_g[1,1,h,g]*B) - #entries from births
          (total_out_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) + #exists from aging out and death
          (sum(eta_i_h_g[HIV_SET, h, g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,HIV_SET,g]])) - #entries into HIV compartment
          (sum(eta_i_h_g[h,HIV_SET, g])*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) #exit from HIV compartment
      )
    }
  }
  
  #############TB compartment 2 Equations########
  #Set DR compartment to 1, since not applicable to drug resistant compartments#
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]<-(
        (kappa_t_h_g[1,h,g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])- #entries from on to IPT (and adherence)
          (omega*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])- #exits from off IPT
          ((sum(iota_r*FOI_r))*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) - #exits from infection (diminished for IPT)
          (total_out_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])+ #exists from aging out and death
          (sum(eta_i_h_g[HIV_SET, h,g]*N_t_r_h_g[N_t_r_h_g_ref[2,1,HIV_SET,g]])) - #entries into HIV compartment
          (sum(eta_i_h_g[h,HIV_SET,g])*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) #exit from HIV compartment
      )
    }
  }
  
  #############TB compartment 3 Equations########
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        dN_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]<-(
          (FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])+ #infection from compartment 1
            (iota_r[r]*FOI_r[r]*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])+ #infections from compartment 2
            (zeta*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,DR_SET,h,g]]))+ #re-infection from compartment 4
            (xi*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[7,DR_SET,h,g]]))+ #re-infection from compartment 7
            (zeta*FOI_r[r]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,DR_SET,h,g]]))- #re-infection from compartment 8
            (pi_i_t[3,4]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #from recent to remote infection
            (gamma_r[r]*kappa_t_h_g[3,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) - #from recent to on IPT (and adherence)
            (theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from recent TB infection to active
            (alpha_in_t_r_h_g[3,r,h,g]*B) - #entries from births
            (total_out_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #exists from aging out and death
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
        dN_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]<-(
          (pi_i_t[3,4]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]])- #from recent to remote infection
            (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]])- #re-infection
            (gamma_r[r]*kappa_t_h_g[4,h,g]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) - #onto IPT (and adherence)
            (theta_h[h]*pi_i_t[4,6]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #from remote TB infection to active
            (alpha_in_t_r_h_g[4,r,h,g]*B) - #entries from births
            (total_out_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #exists from aging out and death
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
        dN_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]] <-(
          (gamma_r[r]*kappa_t_h_g[3,h,g]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from recent to on IPT  
            (gamma_r[r]*kappa_t_h_g[4,h,g]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) - # remote to on IPT (and adherence)
            (theta_h[h]*pi_i_t[5,6]*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) - #from TB infection on IPT to active
            (omega*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) - #off of IPT
            (total_out_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) + #exists from aging out and death
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
        dN_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]] <- (
          (theta_h[h]*pi_i_t[3,6]*N_t_r_h_g[N_t_r_h_g_ref[3,r,h,g]]) + #from recent TB infection to active 
            (theta_h[h]*pi_i_t[4,6]*N_t_r_h_g[N_t_r_h_g_ref[4,r,h,g]]) + #from remote TB infection to active 
            (theta_h[h]*pi_i_t[5,6]*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) + #from TB infection on IPT to active
            (theta_h[h]*pi_i_t[8,6]*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) + #from TB after on IPT to active
            (pi_i_t[7,6]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) - #TB relapse rate
            (upsilon_h[h]*pi_i_t[6,7]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) + #from active to recovered
            (alpha_in_t_r_h_g[6,r,h,g]*B) - #entries from births
            (total_out_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) + #total out
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
        dN_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]] <- (
          (upsilon_h[h]*pi_i_t[6,7]*N_t_r_h_g[N_t_r_h_g_ref[6,r,h,g]]) - #from active to recovered
            (pi_i_t[7,6]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) - #relapse rate
            (xi*FOI*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) - #re-infection from compartment 7
            (total_out_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[7,r,h,g]]) + #total out
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
        dN_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]] <- (
          (omega*N_t_r_h_g[N_t_r_h_g_ref[5,r,h,g]]) - #off IPT to after IPT
            (theta_h[h]*pi_i_t[8,6]*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) - #from TB after on IPT to active
            (zeta*FOI*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]])- #reinfection
            (total_out_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[8,r,h,g]]) + #total out
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

out_df_all<-data.frame()
sim_id <-1

for (b1 in beta_1_test){
  for (b2 in beta_2_test){
    
    #to track how long it takes to solve
    start<-Sys.time()
    
    beta_g[1]<-b1
    beta_g[2]<-b2
    
    #pop_init with betas
    pop_init_df_temp<-pop_init_df%>%
      filter(beta_1 == b1,
             beta_2 == b2)
    
    #####N_init - total initial population in each compartment flattened in 1D array for ODE solver #####
    #arrange pop init df, TB-->DR-->HIV-->G
    pop_init_df_temp<-pop_init_df_temp%>%
      arrange(G_compartment)%>%
      arrange(HIV_compartment)%>%
      arrange(DR_compartment)%>%
      arrange(TB_compartment)
    
    N_init <- pop_init_df_temp$value
    
    #add in mort calc placeholders
    calibration_mort_states<-c('HIV_neg_male', 'HIV_neg_female', 'HIV_pos_male', 'HIV_pos_female')
    N_init <-c(N_init, 0,0,0,0)
    names(N_init) <- c(pop_init_df_temp$compartment_id, calibration_mort_states)

    out_df<-as.data.frame(ode(times = TT_SET, y = N_init,
                              func = open_seir_model, method = 'lsoda',
                              parms = NULL))

    out_df<- cbind(year = as.integer(start_yr+out_df$time), sim_id = rep(sim_id, times = nrow(out_df)),
                   beta_1 = rep(b1, times = nrow(out_df)), beta_2 = rep(b2, times = nrow(out_df)),
                   out_df)

    out_df_all<-rbind(out_df_all, out_df)

    sim_id <- sim_id + 1

    print(Sys.time()-start)

  }
}




#####Graphs that describe model states overtime####
state_prog_df<-out_df%>%
  select(-c(calibration_mort_states))

#melt out all df for easy manipulation
state_prog_df <-melt(state_prog_df, id.vars = c("time", "year",
                                               'sim_id', "beta_1", "beta_2"))
state_prog_df <- cbind(state_prog_df, 
                     data.frame(do.call('rbind', 
                                        strsplit(as.character(state_prog_df$variable),
                                                 '_',fixed=TRUE))))
names(state_prog_df)[names(state_prog_df) == "X2"] <- "TB_compartment"
names(state_prog_df)[names(state_prog_df) == "X3"] <- "DR_compartment"
names(state_prog_df)[names(state_prog_df) == "X4"] <- "HIV_compartment"
names(state_prog_df)[names(state_prog_df) == "X5"] <- "G_compartment"
state_prog_df<-state_prog_df%>%select(-c('X1'))
state_prog_df$month <- round(((state_prog_df$time%%1)*(12)+1),0)

outdir_state_prog <- paste0(here(),'/model_outputs/pre_calibration_tests/Apr23/beta14_lowmort')
setwd(outdir_state_prog)

######TB states########
TB_all_overtime<-state_prog_df%>%
  group_by(TB_compartment, time)%>%
  summarise(total_in_compartment = sum(value))%>%
  filter(TB_compartment != 'pop')%>%
  mutate(time = time+start_yr)

tb_all_graph<-ggplot(data = TB_all_overtime%>%filter(TB_compartment == 7), 
                     mapping = aes(x = time, 
                                   y = total_in_compartment, 
                                   color = TB_compartment))+
  geom_line()+
  labs(title = 'All TB compartments overtime')

png("tb_all_overtime_after_warmup.png")
print(tb_all_graph)
dev.off()

TB_recovered_overtime<-state_prog_df%>%
  filter(TB_compartment)%>%
  group_by(HIV_compartment, time)%>%
  summarise(total_in_compartment = sum(value))%>%
  mutate(time = time+start_yr)

ggplot(data = TB_recovered_overtime, 
       mapping = aes(x = time, 
                     y = total_in_compartment, 
                     color = HIV_compartment))+
  geom_line()+
  labs(title = 'TB compartments overtime')

tb_active_graph_summarised<-ggplot()+
  +
  labs(title = 'Active TB overtime')

png("tb_active_overtime_after_warmup.png")
print(tb_active_graph_summarised)
dev.off()

TB_prev_overtime_gender_HIV<-state_prog_df%>%
  group_by(time, G_compartment)%>%
  mutate(total_g_compartment = sum(value))%>%
  ungroup()%>%
  filter(TB_compartment == 6)%>%
  group_by(time, HIV_compartment, G_compartment)%>%
  mutate(total_pop_in_compartment = sum(value),
         time = time+1990,
         percent_pop_in_compartment = 
           total_pop_in_compartment/total_g_compartment,
         rate = percent_pop_in_compartment*100000)%>%
  ungroup()%>%
  group_by(year, HIV_compartment, G_compartment)%>%
  mutate()

for (g in G_SET){
  
  file_name <- paste0('ActiveTB_g_compartment_', g, '.png')
  
  gender_name = if_else(g == 1, 'Males', 'Females')
  
  graph_temp <-ggplot() +
    geom_line(TB_prev_overtime_gender_HIV%>%filter(G_compartment == g),
              mapping = aes(x = time, y = rate, 
                            color = HIV_compartment), size = 1)+
    labs(title = paste0('TB Prevalence, rate per 100K,\nby HIV compartment for ', gender_name))+
    scale_y_continuous(name=paste0("Rate per 100k, ", gender_name), limits=c(0, 600))+
    scale_x_continuous(name = 'Time', limits = c(1990, 2017), breaks=(seq(1990, 2017, 5)))+
    scale_color_manual(values=c("red", "#56B4E9", "purple", "green4"))+
    #theme(legend.position='none', text = element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  png(file_name)
  print(graph_temp)
  dev.off()
}


#####Calibration Calculations######
total_pop_in_gender_df<-state_prog_df%>%
  group_by(time, year, G_compartment)%>%
  summarise(total_gender_pop = sum(value))%>%
  ungroup()%>%
  group_by(year, G_compartment)%>%
  summarise(total_gender_pop = median(total_gender_pop))


mort_calibration_df<-out_df%>%
  select(c('year', 'time', 'HIV_neg_male', 'HIV_neg_female', 'HIV_pos_male', 'HIV_pos_female'))%>%
  group_by(year)%>%
  summarise(cum_mort_hiv_neg_male = max(HIV_neg_male),
         cum_mort_hiv_neg_female= max(HIV_neg_female),
         cum_mort_hiv_pos_male= max(HIV_pos_male),
         cum_mort_hiv_pos_female= max(HIV_pos_female))

mort_calibration_df<-melt(mort_calibration_df, id = c('year'))
mort_calibration_df<-mort_calibration_df%>%
  group_by(variable)%>%
  mutate(mort_in_year = value - lag(value))%>%
  ungroup()%>%
  mutate(mort_in_year = if_else(year == 1990, value, mort_in_year))%>%
  filter(year < 2018)%>%
  mutate(G_compartment = if_else(grepl('female', variable), '2', '1'))%>%
  left_join(total_pop_in_gender_df, by = c('G_compartment', 'year'))%>%
  mutate(mort_rate_per_100K = ((mort_in_year*100000)/total_gender_pop))

for (g in G_SET){
  
  file_name <- paste0('TB_moratlity_g_compartment_', g, '.png')
  
  gender = if_else(g == 1, 'Males', 'Females')
  
  graph_temp <- ggplot(data = mort_calibration_df%>%filter(G_compartment == g), 
                       mapping = aes(x = year, y = mort_rate_per_100K))+
    geom_line(aes(colour = variable), size = 1)+
    labs(title = paste0('Deaths, rate per 100K, ', gender))+
    #scale_y_continuous(name="rate", breaks=seq(from = 0, to = 300, by = 20))+
    scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
    scale_color_manual(values=c("green", 'purple'))+
    ylim(0, 600)
  
  png(file_name)
  print(graph_temp)
  dev.off()
}

#####overlap model outputs with calibration#####
setwd(indir)
GBD_calibration_df<-read.csv('calibration_rates_df.csv')
mort_calibration_df_all<-mort_calibration_df%>%
  select(c('year', 'variable', 'mort_rate_per_100K'))

calibration_group_name<-c()

for (n in 1:nrow(mort_calibration_df_all)){
  if(mort_calibration_df_all$variable[n] == 
     'cum_mort_hiv_neg_male'){
    calibration_group_name<-c(calibration_group_name, 'TB_only_Male')
  } else if (mort_calibration_df_all$variable[n] == 
             'cum_mort_hiv_neg_female'){
    calibration_group_name<-c(calibration_group_name, 'TB_only_Female')
  } else if (mort_calibration_df_all$variable[n] == 
             'cum_mort_hiv_pos_female'){
    calibration_group_name<-c(calibration_group_name, 'HIV/TB_coinfection_Female')
  } else {
    calibration_group_name<-c(calibration_group_name, 'HIV/TB_coinfection_Male')
  }
}

mort_calibration_df_all$calibration_group<-calibration_group_name
mort_calibration_df_all<-mort_calibration_df_all%>%
  rename(model_output = mort_rate_per_100K)%>%
  select(c('year', 'calibration_group', 'model_output'))%>%
  left_join(GBD_calibration_df, by = c('year','calibration_group'))%>%
  filter(year < 2017)

setwd(outdir_state_prog)
for (g in unique(mort_calibration_df_all$sex_name)){
  
  file_name <- paste0('TB_moratality_calibration', g, '.png')
  
  df_temp <- mort_calibration_df_all%>%filter(sex_name == g)
  df_temp_model<-df_temp
  df_temp_model$calibration_group<-paste0('Model--', 
                                          df_temp_model$calibration_group)
  df_temp_GBD<-df_temp
  df_temp_GBD$calibration_group<-paste0('GBD projections (expected)--', 
                                        df_temp_GBD$calibration_group)
  df1<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[1])
  df2<-df_temp%>%filter(calibration_group == unique(df_temp$calibration_group)[2])
  
  graph_temp <- ggplot()+
    geom_ribbon(data=df1,aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE,fill = "plum2")+
    geom_ribbon(data=df2,aes(x = year, ymin = min_rate, ymax = max_rate), 
                inherit.aes = FALSE,fill = "darkseagreen1")+
    geom_line(data = df_temp_model, aes(x = year, y = model_output, 
                  color = calibration_group),
                  linetype="dashed", size = 1.2)+
    geom_line(data = df_temp_GBD, aes(x = year, y = expected_rate, 
                  color = calibration_group))+
    scale_color_manual(values=c('green', 'purple', 'darkgreen', 'darkorchid4'))+
    labs(title = paste0('Deaths, rate per 100K, ', g))+
    scale_x_continuous(name = 'time', breaks=seq(from = 1990, to = 2017, by = 5))+
    scale_y_continuous(name = 'mortality rate', limits = c(0, 600), breaks=(seq(0, 600, 100)))+
    theme(legend.position='none', text = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  png(file_name, width = 480, height = 480)
  print(graph_temp)
  dev.off()
}







  
  
  #data.frame(diff(as.matrix(mort_df$cum_mort_hiv_pos_female)))
# 
# HIV_overtime<-out_df_melt%>%
#   mutate(HIV_temp = if_else(HIV_compartment > 1, value, 0))%>%
#   group_by(time, G_compartment)%>%
#   mutate(total_in_g_compartment = sum(value))%>%
#   group_by(time, G_compartment, total_in_g_compartment)%>%
#   summarise(total_hiv_pos = sum(HIV_temp))%>%
#   mutate(hiv_prev = total_hiv_pos/total_in_g_compartment)
# 
# #png("HIV_prev_overtime_reactivation.png")
# #hivplot1 <- ggplot(data = HIV_overtime, 
#                    mapping = aes(x = time, y = hiv_prev, color = G_compartment))+
#   geom_line()+
#   lims(y = c(0,1), x = c(0,27))+
#   labs(title = 'REACTIVATION (based on remote rates) - HIV Prev Over Time')
# print(hivplot1)
# dev.off()
# 
# 
# ART_coverage_overtime<-out_df_melt%>%
#   group_by(time, HIV_compartment, G_compartment)%>%
#   summarise(total = sum(value))%>%
#   filter(HIV_compartment > 1)%>%
#   group_by(time, G_compartment)%>%
#   mutate(total_PLHIV = sum(total))%>%
#   filter(HIV_compartment == 4)%>%
#   mutate(ART_coverage = total/total_PLHIV)
# 
# png("ART_coverage_overtime_reactivation.png")
# artcov_plot <- ggplot(data = ART_coverage_overtime, 
#                       mapping = aes(x = time, y = ART_coverage, color = G_compartment))+
#   geom_line()
# print(artcov_plot)
# dev.off()
# 
# 
# 
# 
# # 
# #   
# # 
# # 
# # 
# # 
# # #hiv compartments overtime test to compare with caras coverages
# # out_df_by_HIV<-out_df_melt%>%
# #   group_by(HIV_compartment, G_compartment, time)%>%
# #   summarise(total_hiv_pop = sum(value))
# # 
# # #testing hiv prevalence
# # ggplot(data = out_df_by_HIV, 
# #        mapping = aes(x = time, y = total_pop, fill = HIV_compartment))+
# #   geom_area()
# # 
# # 
# # #filter only active TB compartments, since that is what we are calibrating to
# # out_df_TB_active<-out_df_melt%>%
# #   filter(TB_compartment == 6)%>%
# #   mutate(calibration_group_name = if_else((HIV_compartment == 1 & G_compartment == 1),
# #                                      'TB_only_Male',
# #                                      if_else((HIV_compartment != 1 & G_compartment == 1),
# #                                              'HIV/TB_coinfection_Male',
# #                                              if_else((HIV_compartment == 1 & G_compartment == 2),
# #                                              'TB_only_Female',
# #                                              'HIV/TB_coinfection_Female'))),
# #          calibration_group_id = if_else((HIV_compartment == 1 & G_compartment == 1),
# #                                      1,
# #                                      if_else((HIV_compartment != 1 & G_compartment == 1),
# #                                              2,
# #                                              if_else((HIV_compartment == 1 & G_compartment == 2),
# #                                                      3,
# #                                                      4))))
# # 
# # #so I do not need to call on mort_param_func too many times
# # out_df_TB_active<-out_df_TB_active%>%
# #   arrange(year)
# # 
# # #testing TB prev overtime
# # test<-out_df_TB_active%>%
# #   group_by(G_compartment, HIV_compartment, time)%>%
# #   summarise(value = sum(value))%>%
# #   mutate(ID = paste0('HIV_', HIV_compartment, '_G_', G_compartment))
# # 
# # ggplot(data = test, mapping = aes(x = time, y = value, fill = ID))+
# #   geom_area()
# # 
# # 
# # TB_grouping_test_df<-out_df_melt%>%
# #   group_by(TB_compartment, time)%>%
# #   summarise(total_pop = sum(value))
# # 
# # 
# # ggplot(data = TB_grouping_test_df%>%filter(TB_compartment == 6), mapping = aes(x = time, y = total_pop,
# #                                                  fill = TB_compartment))+
# #   geom_area()
# # 
# # counter <-1
# # mort_est <- rep(0, times = nrow(out_df_TB_active))
# # 
# # for (yr in start_yr:end_yr){
# #   temp<-out_df_TB_active%>%
# #     filter(year == yr)
# #   mu_t_h_g<-mort_param_func(yr)
# #   for (row in 1:nrow(temp)){
# #     hiv <- as.integer(temp[row, 'HIV_compartment'])
# #     gender <-as.integer(temp[row, 'G_compartment'])
# #     pop<-as.double(temp[row, 'value'])
# #     mort_rate<-(mu_t_h_g[6,hiv,gender]*(1/12))
# #     
# #     mort_est[counter]<-(mort_rate*pop)
# #     counter <- counter + 1
# #   }
# # }
# # 
# # out_df_TB_active$mort_est <- mort_est
# # 
# # #group and combine data for calibration
# # setwd(indir)
# # calibration_rates_df<-read.csv('calibration_rates_df.csv')
# # 
# # calibration_df<-out_df_TB_active%>%
# #   group_by(calibration_group_id, calibration_group_name, 
# #            year)%>%
# #   summarise(model_rate = sum(mort_est))%>%
# #   left_join(calibration_rates_df, by = c('calibration_group_name', 'year'))%>%
# #   filter(year < 2017)
# #   #mutate(time = time+1990) #so that time graphs go from 1990 - 2017
# # 
# # calibration_df$within<-if_else((calibration_df$model_rate<=calibration_df$max_rate)&(calibration_df$model_rate>=calibration_df$min_rate),
# #                                1, 0)
# # calibration_df$diff <-calibration_df$model_rate-calibration_df$expected_rate
# # calibration_df$mse <-(calibration_df$diff)^2
# # 
# # 
# # 
