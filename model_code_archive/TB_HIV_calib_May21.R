#model May 21
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

#############Set in directory and out directory###########
###########Make sure epi_model_HIV_TB.Rproj is open, otherwise will need to change wd manually########
current_date <- gsub('-', '_', Sys.Date())
warmup_date<-'2021_05_24'
indir_params <- paste0(here(),'/param_files')
indir_calib_specific <- paste0(here(),'/model_outputs/calibration_data_sets/', warmup_date)
outdir<-paste0(here(),'/model_outputs/calibration/', current_date)

#will give warning if wd already exists
#ignore warning!
dir.create(file.path(outdir))

##########Read in data##########

#these data frames change depending on combin of calib parameters
setwd(indir_calib_specific)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')
pop_init_df <- read.csv('pop_init_df.csv')
mort_df <- read.csv('mort_df.csv')
mort_df<-mort_df%>%
  select(c(year, TB_compartment, HIV_compartment, G_compartment, mort_rate, level))
mort_df$mort_level_id<-paste0(mort_df$HIV_compartment, '-', mort_df$level)

#these data frames change depending on combin of calib parameters
setwd(indir_params)
param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
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

################ Define sets ###############################

#8 TB set description (TB)#
#1:Uninfected, not on IPT;
#2:Uninfected, on IPT; 
#3:LTBI, infected recently (within the past two-years)
#4: LTBI, infected remotely (more than two-years ago)
#5: LTBI, on IPT
#6: Active
#7: Recovered/Treated
#8: LTBI, after IPT

TB_SET<-1:8

#4 HIV compartments description (HIV)#
#1 : HIV Negative
#2 : HIV Positive CD4 > 200 - No ART
#3 : HIV Positive CD4 =<: 200 - No Art
#4 : HIV Positive - ART 

HIV_SET<-1:4

#2 Drug Resistance compartments description (DR)#
#1 : Drug Susceptible
#2 : Multi Drug Resistant

DR_SET<-1:2

#2 Gender compartments description (G)#
#1: Male
#2: Female

G_SET<-1:2

#######Parameter extraction########

######## Parameters that impact Force of Infection #######

#beta_g - Number of effective contacts for TB transmission per infectious year#
beta_g <- param_df%>%
  filter(notation == 'beta')%>%
  arrange(G_compartment)

beta_g <- beta_g$Reference_expected_value

#phi_h - Relative transmissibility of TB in HIV pops#
phi_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'phi',
           HIV_compartment == h)
  phi_h[h] <- temp$Reference_expected_value
}

#varepsilon_g - Fraction of new TB infections that are MDR-TB #
varepsilon_g <- param_df%>%
  filter(notation == 'varepsilon')%>%
  arrange(G_compartment)

varepsilon_g <- varepsilon_g$Reference_expected_value

#iota_r - Indicator for whether infection with given TB strain can occur while on IPT by DR compartment#
iota_r <- param_df%>%
  filter(notation == 'iota')
iota_r <- iota_r$Reference_expected_value

#zeta - Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection#
zeta <- param_df%>%
  filter(notation == 'zeta')
zeta <-zeta$Reference_expected_value

#########Parameters that Describe TB progression ######

#kappa_t_h_g_p - Rate of IPT initiation, per year#
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

#omega - Rate of moving off of IPT, per year #
omega <-param_df%>%filter(notation == 'omega')
omega <- omega$Reference_expected_value

#pi_i_t - Base rates of TB progression#
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

#theta_h - relative risk of TB progression#
theta_h <-array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'theta',
           HIV_compartment == h)
  
  theta_h[h]<- temp$Reference_expected_value
}

#upsilon_h - increased time infectious due to delayed treatment#
upsilon_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'upsilon',
           HIV_compartment == h)
  upsilon_h[h] <- temp$Reference_expected_value
}

#gamma_r -indicator if DR compartment can move onto after IPT#
gamma_r <- c(1,0)

#increased risk of progression after re-infection for recovered populations#
xi<-param_df%>%
  filter(notation == 'xi')
xi <-xi$Reference_expected_value

#######Parameters that describe HIV progression########
#eta_i_h_g rate HIV transitions#
eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                            length(HIV_SET), 
                            length(G_SET)))

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
      
#########Parameters for death and aging rates ###########
#mu_t_h_g - mortality rates#

mort_df_calib_temp<-mort_df #to store current calib mort parameters
mu_t_h_g <- array(0, dim = c(length(TB_SET), length(HIV_SET), length(G_SET)))

mort_param_func <-function(yr){
  mu_t_h_g <- array(0, dim = c(length(TB_SET), length(HIV_SET), length(G_SET)))
  
  for (t in TB_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- mort_df_calib_temp%>%
          filter(year == yr,
                 TB_compartment == t,
                 HIV_compartment == h,
                 G_compartment == g)
        
        mu_t_h_g[t,h,g] <- unique(temp$mort_rate)
      }
    }
  }
  return(mu_t_h_g)
}

#alpha_in_t_h_g - Proportion of population that enters each compartment#
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                      length(DR_SET), 
                                      length(HIV_SET), 
                                      length(G_SET)))

birth_rate_param_func<-function(yr){
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                      length(DR_SET), 
                                      length(HIV_SET), 
                                      length(G_SET)))

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

#alpha_out - Rate of exit from the population#
alpha_out <- param_df%>%filter(notation == 'alpha^out')
alpha_out <- alpha_out$Reference_expected_value


#############Pre-processing parameter equations, for ease of use in ode solver#######

#total_out_t_r_h - total amount leaving from compartment#
total_out_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET))

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

#######EQUATIONS THAT DESCRIBE TB AND HIV PROG IN DESOLVE######
tb_hiv_prog_calibration_model <- function(time, N_t_r_h_g, parms){
  
  #print(time)
  
  #initiation delta array, +4 for tracking mort
  dN_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET)+4)
  
  #time varying parameters
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
  
  if (current_yr > last_year){
    #HIV transitions
    eta_i_h_g<<-HIV_transitions_param_func(current_yr)
    
    #alpha in proportions
    alpha_in_t_r_h_g<<-birth_rate_param_func(current_yr)
    
    #deaths
    mu_t_h_g<<-mort_param_func(current_yr)
    total_out_t_r_h_g<<-total_out_param_func(mu_t_h_g)
    
    last_year<<-current_yr
  }
  
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


########Feed paramters into desolve########
#feed in to solve, evaluation time intervals (TT_SET), 
#initial population (N_init)
#differential equations and dynamic parameter (open_seir_model)
#specify differential solver (lsoda)
# read params from global (NULL)
mort_calib_array<-c()
setwd(outdir)
last_year<-1989

for(n in 1:nrow(sim_calibration_ref_df)){
  last_year<-1989
  print(n)
  #to track how long it takes to solve each run
  start<-Sys.time()
  sim_attributes<-sim_calibration_ref_df[n, ]
  sim_id_temp2<-as.integer(sim_attributes[1])
  beta_g<<-as.double(sim_attributes[2:3])
  mort_rate_levels_test<-sim_attributes[4:7]
  
  mort_df_calib_temp<<-mort_df%>%
    filter(mort_level_id %in% mort_rate_levels_test)
  
  #pop init (by betas)
  pop_init_df_temp<-pop_init_df%>%
    filter(sim_id == sim_id_temp2)
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
                            func = tb_hiv_prog_calibration_model, method = 'lsoda',
                            parms = NULL))
  
  out_df<- cbind(year = as.integer(start_yr+out_df$time), sim_id = rep(sim_id_temp2, times = nrow(out_df)),
                 out_df)
  
  write.csv(out_df, file = paste0('out_df_sim_id_', sim_id_temp2, '.csv'))
  
  #to track the current simulation
  print(Sys.time()-start)
}

