#Last updated Feb 25 2022
#warmup 1940-1990
#calibration period 1990-2017

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr',
         'bitops'), require, character.only=T)

####HYAK SPECIFIC####
start_sim<-1
end_sim<-2

#######Define working directories#######
indir<-paste0(here(), '/param_files')
outdir<-paste0(here(), '/calibration_runs/test')

##########Read in param files##########
setwd(indir)

#only keep sims in run
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')%>%
  filter(sim_id >= start_sim,
         sim_id <= end_sim)

#updated for each sim
sim_calibration_ref_df_current<-data.frame()

pop_init_df <- read.csv('pop_init_df.csv')
hiv_incidence_df<-read.csv('hiv_inc_df.csv')
art_coverage_df<-read.csv('art_coverage_df.csv')
art_prop_eligible_df<-read.csv('ART_prop_eligible_df.csv')
ipt_initiation_df<-read.csv('ipt_initiation_df.csv')
birth_perc_df<-read.csv('birth_perc_df.csv')
base_mort_df<-read.csv('base_mort_df.csv')

param_df <- read_excel('KZN_SA_model_parameters.xlsx', 
                       sheet = 'model_matched_parameters')

#clean dataframe column names for consistency
names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))
#make sure all compartments are integer type for proper indexing#
param_df$TB_compartment<-as.integer(param_df$TB_compartment)
param_df$DR_compartment<-as.integer(param_df$DR_compartment)
param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
param_df$G_compartment<-as.integer(param_df$G_compartment)

#make sure all values available are numeric
param_df$value<-as.numeric(param_df$value) #NAs okay

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

########INITIALIZE POPULATION AND ID LOCATION OF CURRENT POPULATION AMOUNTS######
pop_init_df<-pop_init_df%>%
  arrange(G_compartment)%>%
  arrange(HIV_compartment)%>%
  arrange(DR_compartment)%>%
  arrange(TB_compartment)

#initialize population values
N_init <- pop_init_df$total_pop
names(N_init)<-pop_init_df$compartment_id

#N_t_r_h_g - matrix that identifies the location of compartment in 1D array#
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

######## Parameters That Impact Force of Infection #######

#beta_g - Number of effective contacts for TB transmission per infectious year#
#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
beta_g<-array(0, dim = length(G_SET)) 

#pull range for effective contacts
beta_g_temp <- param_df%>%
  filter(notation == 'beta')%>%
  arrange(G_compartment)

beta_g_max <- beta_g_temp$max
beta_g_min<-beta_g_temp$min

rm(beta_g_temp)

#phi_h - Relative transmissibility of TB in HIV pops#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
phi_h <- array(0, dim = length(HIV_SET))

#except for baseline (not calibrated)
phi_h[1]<-1

#pull range for relative transmissibility
phi_h_max <- array(0, dim = length(HIV_SET))
phi_h_min <- array(0, dim = length(HIV_SET))

for (h in 2:4){
  phi_h_temp <- param_df%>%
    filter(notation == 'phi',
           HIV_compartment == h)
  phi_h_max[h] <- phi_h_temp$max
  phi_h_min[h] <- phi_h_temp$min
}

rm(phi_h_temp)

#varepsilon- Fraction of new TB infections that are MDR-TB #

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
varepsilon<-0

#MDR infection only calibrated
varepsilon_temp <- param_df%>%
  filter(notation == 'varepsilon')

#pull range for relative transmissibility
varepsilon_max <- varepsilon_temp$max
varepsilon_min <- varepsilon_temp$min

rm(varepsilon_temp)

#iota_r - Indicator for whether infection with given TB strain can occur while on IPT by DR compartment#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
iota_r<-array(0, dim = length(DR_SET))

#only calibrate for MDR strains (0 for DR strains)
iota_r_temp <- param_df%>%
  filter(notation == 'iota')%>%
  filter(DR_compartment == 2)

iota_r_max<-array(0, dim = length(DR_SET))
iota_r_max[2] <- iota_r_temp$max

iota_r_min<-array(0, dim = length(DR_SET))
iota_r_min[2] <- iota_r_temp$min

rm(iota_r_temp)

#increased risk of progression after re-infection for recovered populations#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
xi <- 0

xi_temp<-param_df%>%
  filter(notation == 'xi')

xi_min <-xi_temp$min
xi_max <-xi_temp$max

rm(xi_temp)

#zeta - Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
zeta<-0

zeta_temp <- param_df%>%
  filter(notation == 'zeta')

zeta_max <-zeta_temp$max
zeta_min <-zeta_temp$min

rm(zeta_temp)

#########Parameters that Describe TB progression ######

#omega - Rate of moving off of IPT, per year #

#not calibrated
#constant over all time intervals
omega_temp <-param_df%>%
  filter(notation == 'omega')
omega <- omega_temp$value
rm(omega_temp)

#pi_i_t - Base rates of TB progression#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
pi_i_t <- array(data = 0, c(length(TB_SET), length(TB_SET)))

pi_temp <- param_df%>%
  filter(notation == 'pi')

#EXCEPT base recovery rates (based off of duration of treatment)
#not calibrated
#constant over all time intervals
pi_temp67<-pi_temp%>%
  filter(TB_compartment == 67)

pi_i_t[6,7] <-pi_temp67$value

rm(pi_temp67)

#calibrated
#constant over all time intervals

#pi_3_4 recently infection period
#pi_3_6 base rates recent to active TB
#pi_4_6 base rates remote to active TB
#pi_5_6 base rates LTBI, on IPT to active TB
#pi_7_6 base rates relapse
#pi_8_6 base rates after LTBI, on IPT to active TB

pi_i_t_max <-array(data = 0, c(length(TB_SET), length(TB_SET)))
pi_i_t_min <-array(data = 0, c(length(TB_SET), length(TB_SET)))

for (t_from in TB_SET){
  for (t_to in TB_SET){
    num_temp = (t_from*10) + t_to
    
    temp <- param_df%>%
      filter(TB_compartment == num_temp)
    
    if (nrow(temp) == 1){
      pi_i_t_max[t_from,t_to] <- temp$max
      pi_i_t_min[t_from,t_to] <- temp$min
    }
  }
}

#theta_h - relative risk of TB progression#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
theta_h <-array(0, dim = length(HIV_SET))

theta_h[1]<-1 #reference value not calibrated

theta_h_min<-array(0, dim = length(HIV_SET))
theta_h_max<-array(0, dim = length(HIV_SET))

for (h in 2:4){
  theta_h_temp <- param_df%>%
    filter(notation == 'theta',
           HIV_compartment == h)
  
  theta_h_min[h]<- theta_h_temp$min
  theta_h_max[h]<- theta_h_temp$max
}

#upsilon_h - increased time infectious due to delayed treatment#

#calibrated and constant over all time intervals (updated in sim_calibration_update_fun)
#placeholder
upsilon_h <- array(0, dim = length(HIV_SET))

upsilon_h_min <- array(0, dim = length(HIV_SET))

upsilon_h_max <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  upsilon_h_temp <- param_df%>%
    filter(notation == 'upsilon',
           HIV_compartment == h)
  upsilon_h_min[h] <- upsilon_h_temp$min
  upsilon_h_max[h] <- upsilon_h_temp$max
}

#gamma_r -indicator if DR compartment can move onto after IPT#
#not calibrated
#not time varying
gamma_r <- c(1,0)

#kappa_t_h_g_p - Rate of IPT initiation, per year#

#policy parameter not calibrated#
#time varying
kappa_t_h_g <- array(data = 0, c(length(TB_SET), length(HIV_SET), length(G_SET)))
start_ipt_yr <- 2005

IPT_initiation_function<-function(yr){
  
  for (g in G_SET){
    
    IPT_init_temp<-ipt_initiation_df%>%
      filter(gender_id == g,
             year == yr,
             POLICY_ID == 1)
    
    IPT_init_temp1<-IPT_init_temp%>%filter(on_art == 'no')
    
    kappa_t_h_g[c(1,3,4), c(2,3), g]<-IPT_init_temp1$ipt_init_perc
    
    IPT_init_temp2<-IPT_init_temp%>%filter(on_art == 'yes')
    
    kappa_t_h_g[c(1,3,4), 4, g]<-IPT_init_temp2$ipt_init_perc
  }
  return(kappa_t_h_g)
}

#######Parameters that describe HIV progression########
#eta_i_h_g rate HIV transitions placeholder# 
eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                            length(HIV_SET), 
                            length(G_SET)))

#calibrated and changes overtime according to HIV_transitions_param_func
HIV_transitions_param_func<-function(yr, N_t_r_h_g){
  eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                              length(HIV_SET), 
                              length(G_SET)))
  
  #hiv incidence
  #calibrated
  #time varying
  #no hiv infections until 1980
  if (yr >= 1980){
    for (g in G_SET){
      hiv_incidence_temp<-hiv_incidence_df%>%
        filter(sex == if_else(g == 1, 'Male', 'Female'),
               year == yr)
      if (g == 1){
        eta_i_h_g[1,2,g]<-
          qunif(sim_calibration_ref_df_current$eta_12_1, 
                hiv_incidence_temp$min, 
                hiv_incidence_temp$max)
      } else {
        eta_i_h_g[1,2,g]<-
          qunif(sim_calibration_ref_df_current$eta_12_2, 
                hiv_incidence_temp$min, 
                hiv_incidence_temp$max)
      }
    } 
  } 
  
  #get HIV progression CD4 counts 
  #calibrated
  #constant
    hiv_prog_temp<- param_df%>%
      filter(notation == 'eta',
             HIV_compartment == 23,
             G_compartment == g)
    if(g == 1){
      eta_i_h_g[2,3,g]<-
        qunif(sim_calibration_ref_df_current$eta_23_1, 
              hiv_prog_temp$min, 
              hiv_prog_temp$max)
    } else {
      eta_i_h_g[2,3,g]<-
        qunif(sim_calibration_ref_df_current$eta_23_2, 
              hiv_prog_temp$min, 
              hiv_prog_temp$max)
    }
  
    #ART initiation rates
      
    #get ART initiation rates
    #before 2004 no ART is available
    if(yr >= 2004){
      
      #get current model states (V_h_g in the appendix)
      #time varying based on model states
      hiv_prev_current<-sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4), g]])
      n2_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 2, g]])/hiv_prev_current
      n3_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 3, g]])/hiv_prev_current
      n4_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 4, g]])/hiv_prev_current
      
      #get percent eligible (calibrated and time varying)
      prop_eligible<-0
      prop_eligible_temp<-art_prop_eligible_df%>%
        filter(gender_id == g,
               year == yr)
      
      #after 2016 all are eligible
      if (yr < 2016){
        if(g == 1){
          prop_eligible<-
            qunif(sim_calibration_ref_df_current$prop_eligible_males, 
                  prop_eligible_temp$min, 
                  prop_eligible_temp$max)
        } else {
          prop_eligible<-
            qunif(sim_calibration_ref_df_current$prop_eligible_females, 
                  prop_eligible_temp$min, 
                  prop_eligible_temp$max)
        }
      } else {
        prop_eligible<-1
      }
      
      #ART coverage (not calibrated and time varying)
      art_coverage_temp<-art_coverage_df%>%
        filter(sex == if_else(g == 1, 'Male', 'Female'),
               year == yr,
               policy_id == 1)
      
      #calculate general art initiation rate (eta_g^all in the appendix)
      art_coverage_temp<-art_coverage_temp$art_coverage
      
      #calculate overall art initation rate (all CD4 counts)
      art_initiation_rate_all_numerator<-art_coverage_temp-n4_prop_current_yr
      
      art_initiation_rate_all_denominator<-(prop_eligible*n2_prop_current_yr)+(n3_prop_current_yr)
      
      art_initiation_rate_all<-art_initiation_rate_all_numerator/
        art_initiation_rate_all_denominator
      
      eta_i_h_g[2,4,g]<-art_initiation_rate_all*prop_eligible
      eta_i_h_g[3,4,g]<-art_initiation_rate_all
    }
  return(eta_i_h_g)
}

  
#########Parameters for death and aging rates ###########
#mu_t_h_g - mortality rates#

#time varying and calibrated
#placeholder
mu_t_h_g <- array(0, dim = c(length(TB_SET), length(HIV_SET), length(G_SET)))

mu_hiv_increase_temp<-param_df%>%
  filter(notation == 'mu_hiv_increase')%>%
  arrange(HIV_compartment)

mu_tb_hiv_increase_temp<-param_df%>%
  filter(notation == 'mu_tb_hiv_increase')%>%
  arrange(HIV_compartment)

mu_hiv_increase_min<-mu_hiv_increase_temp$min
mu_hiv_increase_max<-mu_hiv_increase_temp$max
mu_tb_hiv_increase_max<-mu_tb_hiv_increase_temp$max
mu_tb_hiv_increase_min<-mu_tb_hiv_increase_temp$min

mort_param_func <-function(yr){
  
  mu_t_h_g <- array(0, dim = c(length(TB_SET), length(HIV_SET), length(G_SET)))
  
  #pull base mortality rates under current calibration and year
  for (g in G_SET){
    if (yr < 1990){
      base_mort_temp <- base_mort_df%>%
        filter(year == 1990,
               sex == if_else(g == 1, 'Male', 'Female'))
    } else {
      base_mort_temp <- base_mort_df%>%
        filter(year == yr,
               sex == if_else(g == 1, 'Male', 'Female'))
    }
    
    increased_hiv_only_mort <-sim_calibration_ref_df_current[31:33]
    increased_tb_hiv_mort <-sim_calibration_ref_df_current[34:37]
    
    if (g == 1){
      base_mort_calib_temp<-
        qunif(sim_calibration_ref_df_current$mort_base_1, 
              base_mort_temp$min, 
              base_mort_temp$max)
    } else {
      base_mort_calib_temp<-
        qunif(sim_calibration_ref_df_current$mort_base_2, 
              base_mort_temp$min, 
              base_mort_temp$max)
    }
    
    for (h in HIV_SET){
      #first update hiv only values
      for (t in c(1:5, 7:8)){
        if (h > 1){
          #pull increased mort
          increased_mort_temp<-qunif(p = as.numeric(increased_hiv_only_mort[h-1]), 
                                     min = as.numeric(mu_hiv_increase_min[h]), 
                                     max = as.numeric(mu_hiv_increase_max[h]))
          
          mu_t_h_g[t,h,g]<-increased_mort_temp*base_mort_calib_temp
        } else {
          #no hiv no tb
          mu_t_h_g[t,1,g]<-base_mort_calib_temp
        }
      }
      #tb and/or hiv values
      #pull increased mort
      increased_mort_temp<-
        qunif(p = as.numeric(increased_tb_hiv_mort[h]), 
              min = as.numeric(mu_tb_hiv_increase_min[h]), 
              max = as.numeric(mu_tb_hiv_increase_max[h]))
      mu_t_h_g[6,h,g]<-increased_mort_temp*base_mort_calib_temp
    }
  }
  return(mu_t_h_g)
}

#alpha_in_t_h_g - Proportion of population that enters each compartment#
#not calibrated
#not time varying
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                      length(DR_SET), 
                                      length(HIV_SET), 
                                      length(G_SET)))

for (t in TB_SET){
  for (r in DR_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- (unique(birth_perc_df$prop_of_pop[(birth_perc_df$TB_compartment == t)&
                                                    (birth_perc_df$DR_compartment == r)&
                                                    (birth_perc_df$HIV_compartment == h)&
                                                    (birth_perc_df$G_compartment == g)]))
        
        alpha_in_t_r_h_g[t,r,h,g] <- temp
      }
    }
  }
}

#alpha_out - Rate of exit from the population#
alpha_out <- param_df%>%filter(notation == 'alpha^out')
alpha_out <- alpha_out$value

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


#######EQUATIONS THAT DESCRIBE TB AND HIV PROG IN DESOLVE######
tb_hiv_prog_calibration_model <- function(time, N_t_r_h_g, parms){
  
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
    FOI_2_g[g] <-(varepsilon*FOI_1_g[g])/(1-varepsilon)
  }
  
  FOI_r <- c(sum(FOI_1_g), sum(FOI_2_g))
  FOI <- sum(FOI_r)
  
  #write year parameter for parameters that change over time
  current_yr <-as.integer(start_yr+time)
  print(current_yr)
  
  if(current_yr > last_year){
    if (current_yr > 1980){
      #HIV transitions start in 1980
        eta_i_h_g<<-HIV_transitions_param_func(current_yr, N_t_r_h_g)
    }
    
    #deaths
    mu_t_h_g<<-mort_param_func(current_yr)
    total_out_t_r_h_g<<-total_out_param_func(mu_t_h_g)
    
    #only start IPT initiations in 2005
    if(current_yr > 2004){
      kappa_t_h_g<<-IPT_initiation_function(current_yr)
    }
    last_year<<-current_yr
  }
  
  #track mort
  dN_t_r_h_g[129]<-mu_t_h_g[6,1,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,1,1]]) #tb active, hiv -, males
  
  dN_t_r_h_g[130]<-mu_t_h_g[6,1,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,1,2]]) #tb active, hiv -, females
  
  dN_t_r_h_g[131]<-mu_t_h_g[6,2,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,2,1]])+ #tb active, hiv+, CD4 > 200, males
    mu_t_h_g[6,3,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,3,1]])+ #tb active, hiv+, CD4 <= 200, males
    mu_t_h_g[6,4,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,4,1]]) #tb active, hiv+, on ART, males
  
  dN_t_r_h_g[132]<-mu_t_h_g[6,2,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,2,2]])+ #tb active, hiv+, CD4 > 200, females
    mu_t_h_g[6,3,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,3,2]])+ #tb active, hiv+, CD4 <= 200, females
    mu_t_h_g[6,4,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,4,2]]) #tb active, hiv+, on ART, females
  
  #track incidence
  #HIV-, males
  dN_t_r_h_g[133]<- (theta_h[1]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,1,1]])) + #from recent TB infection to active 
    (theta_h[1]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,1,1]])) + #from remote TB infection to active 
    (theta_h[1]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,1,1]])) + #from TB infection on IPT to active
    (theta_h[1]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,1,1]])) #from TB after on IPT to active
  
  #HIV- females
  dN_t_r_h_g[134]<- (theta_h[1]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,1,2]])) + #from recent TB infection to active 
    (theta_h[1]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,1,2]])) + #from remote TB infection to active 
    (theta_h[1]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,1,2]])) + #from TB infection on IPT to active
    (theta_h[1]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,1,2]])) #from TB after on IPT to active
  
  #HIV+, males
  #CD4 more
  dN_t_r_h_g[135]<- (theta_h[2]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,2,1]])) + #from recent TB infection to active 
    (theta_h[2]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,2,1]])) + #from remote TB infection to active 
    (theta_h[2]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,2,1]])) + #from TB infection on IPT to active
    (theta_h[2]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,2,1]])) + #from TB after on IPT to active
    #CD4 less
    (theta_h[3]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,3,1]])) + #from recent TB infection to active 
    (theta_h[3]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,3,1]])) + #from remote TB infection to active 
    (theta_h[3]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,3,1]])) + #from TB infection on IPT to active
    (theta_h[3]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,3,1]])) + #from TB after on IPT to active
    #on ART
    (theta_h[4]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,4,1]])) + #from recent TB infection to active 
    (theta_h[4]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,4,1]])) + #from remote TB infection to active 
    (theta_h[4]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,4,1]])) + #from TB infection on IPT to active
    (theta_h[4]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,4,1]])) #from TB after on IPT to active
  
  #HIV+, females
  #CD4 more
  dN_t_r_h_g[136]<-(theta_h[2]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,2,2]])) + #from recent TB infection to active 
    (theta_h[2]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,2,2]])) + #from remote TB infection to active 
    (theta_h[2]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,2,2]])) + #from TB infection on IPT to active
    (theta_h[2]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,2,2]])) + #from TB after on IPT to active
    #CD4 less
    (theta_h[3]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,3,2]])) + #from recent TB infection to active 
    (theta_h[3]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,3,2]])) + #from remote TB infection to active 
    (theta_h[3]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,3,2]])) + #from TB infection on IPT to active
    (theta_h[3]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,3,2]])) + #from TB after on IPT to active
    #on ART
    (theta_h[4]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,4,2]])) + #from recent TB infection to active 
    (theta_h[4]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,4,2]])) + #from remote TB infection to active 
    (theta_h[4]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,4,2]])) + #from TB infection on IPT to active
    (theta_h[4]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,4,2]])) #from TB after on IPT to active
  
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
start_yr = 1940
end_yr = 1941
TT<-end_yr-start_yr
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval)


########Feed paramters into desolve########
setwd(outdir)
last_year<-start_yr-1

for(n in sim_calibration_ref_df$sim_id){
  last_year<-start_yr-1
  print(n)
  #if want to track how long it takes to solve each run
  #start<-Sys.time()
  
  #set calibrated not time varying parameters
  sim_calibration_ref_df_current<<-sim_calibration_ref_df%>%filter(sim_id == n)
  
  beta_g[1]<-qunif(sim_calibration_ref_df_current$beta_1, 
                    min = beta_g_min[1], max = beta_g_max[1])
  beta_g[2]<-qunif(sim_calibration_ref_df_current$beta_2, 
                    min = beta_g_min[2], max = beta_g_max[2])
  phi_h[2]<-qunif(sim_calibration_ref_df_current$phi_2, 
                   min = phi_h_min[2], max = phi_h_max[2])
  phi_h[3]<-qunif(sim_calibration_ref_df_current$phi_3, 
                   min = phi_h_min[3], max = phi_h_max[3])
  phi_h[4]<-qunif(sim_calibration_ref_df_current$phi_4, 
                   min = phi_h_min[4], max = phi_h_max[4])
  varepsilon<-qunif(sim_calibration_ref_df_current$varepsilon,
                     min = varepsilon_min, max = varepsilon_max)
  iota_r[2]<-qunif(sim_calibration_ref_df_current$iota_2,
                    min = iota_r_min[2], max = iota_r_max[2])
  xi<-qunif(sim_calibration_ref_df_current$xi,
             min = xi_min, max = xi_max)
  zeta<-qunif(p = as.numeric(sim_calibration_ref_df_current$zeta),
               min = as.numeric(zeta_min), max = as.numeric(zeta_max))
  pi_i_t[3,4]<-qunif(sim_calibration_ref_df_current$pi_34,
                      min = pi_i_t_min[3,4], max = pi_i_t_max[3,4])
  pi_i_t[3,6]<-qunif(sim_calibration_ref_df_current$pi_36,
                      min = pi_i_t_min[3,6], max = pi_i_t_max[3,6])
  pi_i_t[4,6]<-qunif(sim_calibration_ref_df_current$pi_46,
                      min = pi_i_t_min[4,6], max = pi_i_t_max[4,6])
  pi_i_t[5,6]<-qunif(sim_calibration_ref_df_current$pi_56,
                      min = pi_i_t_min[5,6], max = pi_i_t_max[5,6])
  pi_i_t[7,6]<-qunif(sim_calibration_ref_df_current$pi_76,
                      min = pi_i_t_min[7,6], max = pi_i_t_max[7,6])
  pi_i_t[8,6]<-qunif(sim_calibration_ref_df_current$pi_86,
                      min = pi_i_t_min[8,6], max = pi_i_t_max[8,6])
  theta_h[2]<-qunif(sim_calibration_ref_df_current$theta_2,
                     min = theta_h_min[2], max = theta_h_max[2])
  theta_h[3]<-qunif(sim_calibration_ref_df_current$theta_3,
                     min = theta_h_min[3], max = theta_h_max[3])
  theta_h[4]<-qunif(sim_calibration_ref_df_current$theta_4,
                     min = theta_h_min[4], max = theta_h_max[4])
  upsilon_h[1]<-qunif(sim_calibration_ref_df_current$upsilon_1,
                       min = upsilon_h_min[1], max = upsilon_h_max[1])
  upsilon_h[2]<-qunif(sim_calibration_ref_df_current$upsilon_2,
                       min = upsilon_h_min[2], max = upsilon_h_max[2])
  upsilon_h[3]<-qunif(sim_calibration_ref_df_current$upsilon_3,
                       min = upsilon_h_min[3], max = upsilon_h_max[3])
  upsilon_h[4]<-qunif(sim_calibration_ref_df_current$upsilon_4,
                       min = upsilon_h_min[4], max = upsilon_h_max[4])
  
  #reset hiv transitions
  eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                              length(HIV_SET), 
                              length(G_SET)))
  
  #####N_init - total initial population in each compartment flattened in 1D array for ODE solver #####
  #arrange pop init df, TB-->DR-->HIV-->G
  pop_init_df_temp<-pop_init_df%>%
    arrange(G_compartment)%>%
    arrange(HIV_compartment)%>%
    arrange(DR_compartment)%>%
    arrange(TB_compartment)
  
  N_init <- pop_init_df_temp$total_pop
  
  #add in mort calc placeholders
  calibration_mort_states<-c('TB_mort_HIV_neg_male', 'TB_mort_HIV_neg_female', 
                             'TB_mort_HIV_pos_male', 'TB_mort_HIV_pos_female')
  calibration_incidence<-c('TB_incidence_HIV_neg_male', 'TB_incidence_HIV_neg_female', 
                           'TB_incidence_HIV_pos_male', 'TB_incidence_HIV_pos_female')
  N_init <-c(N_init, 0,0,0,0,0,0,0,0)
  names(N_init) <- c(pop_init_df_temp$compartment_id, 
                     calibration_mort_states,
                     calibration_incidence)
  
  out_df<-as.data.frame(ode(times = TT_SET, y = N_init,
                            func = tb_hiv_prog_calibration_model, method = 'lsoda',
                            parms = NULL))
  
  out_df<- cbind(year = as.integer(start_yr+out_df$time), sim_id = rep(n, times = nrow(out_df)),
                 out_df)
  
  write.csv(out_df, file = paste0('out_df_sim_id_', n, '.csv'), row.names = FALSE)
  
  #if want to track time in the current simulation
  #print(Sys.time()-start)
  #time_track<-c(Sys.time()-start, time_track)
}

