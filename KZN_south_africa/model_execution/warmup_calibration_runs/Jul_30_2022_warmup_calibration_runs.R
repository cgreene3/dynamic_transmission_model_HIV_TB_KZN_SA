#last updated july 30 2022
#1940-2017

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here'), require, character.only=T)

sim_id_current_eval<-1

#Time Horizon and Evaluation intervals (1 month)
start_yr <- 1940
current_yr <- start_yr #tracks for dynamic params
last_year <- current_yr - 1
end_yr <- 2018
TT<-end_yr-start_yr
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval) #tau

#location where input parameters are
indir<-paste0(here(), '/param_files/input_parameters')

#location where want outputs
outdir<-paste0(here(), '/calib_model_test')

#these data frames change depending on combin of calib parameters
setwd(indir)
sim_calibration_ref_df<-read.csv('sim_calibration_ref_df.csv')%>%
  filter(sim_id == sim_id_current_eval)
pop_init_df <- read.csv('pop_init_df_1940.csv')
mort_rate_df <- read.csv('base_mort_df.csv')
birth_perc_df<-read.csv('birth_perc_df_overtime.csv')
param_df <- read_excel("KZN_SA_model_parameters.xlsx", sheet = 'model_matched_parameters')
ipt_initiation_df <-read.csv('ipt_initiation_df.csv')
art_coverage_df <-read.csv('art_coverage_df.csv')
art_prop_eligible_df<-read.csv('art_prop_eligible_df.csv')
hiv_incidence_df<-read.csv('hiv_inc_df.csv')

#clean dataframe column names for consistency
names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))
names(pop_init_df)<-str_replace_all(names(pop_init_df), c(" " = "_" , "-" = "_" ))

#make sure all compartments are integer type for proper indexing#
param_df$TB_compartment<-as.integer(param_df$TB_compartment)
param_df$DR_compartment<-as.integer(param_df$DR_compartment)
param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
param_df$G_compartment<-as.integer(param_df$G_compartment)

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

#######Parameter Extraction########

######## Parameters that impact force of infection #######

#beta_g - Number of effective contacts for TB transmission per infectious year#
beta_g <- array(0, dim = length(G_SET))

for (g in G_SET){
  temp_param_df <- param_df%>%
    filter(notation == 'beta',
           G_compartment == g)
  
  temp_sim_df<-sim_calibration_ref_df%>%
    select(paste0("beta_",g))
  
  beta_g[g] <- as.numeric((temp_sim_df[1]*
                  (temp_param_df$max[1]-temp_param_df$min[1]))+
    temp_param_df$min[1])
}


#phi_h - Relative transmissibility of TB in HIV pops#
phi_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp_param_df <- param_df%>%
    filter(notation == 'phi',
           HIV_compartment == h)
  
  if(temp_param_df$calibrated == "yes"){
    temp_sim_df<-sim_calibration_ref_df%>%
      select(paste0("phi_",h,"."))
    phi_h[h] <- as.numeric((temp_sim_df[1]*
                              (temp_param_df$max[1]-temp_param_df$min[1]))+
                             temp_param_df$min[1])
  } else {
    phi_h[h] <- as.numeric(temp_param_df$value)
  }
  
}

#varepsilon_g - Fraction of new TB infections that are MDR-TB #
varepsilon_g <- array(0, dim = length(G_SET))

for (g in G_SET){
  temp_param_df <- param_df%>%
    filter(notation == 'varepsilon',
           G_compartment == g)
  
  temp_sim_df<-sim_calibration_ref_df%>%
    select(paste0("varepsilon_",g))
  
  varepsilon_g[g] <- as.numeric((temp_sim_df[1]*
                             (temp_param_df$max[1]-temp_param_df$min[1]))+
                            temp_param_df$min[1])
}

#iota_r - Indicator for whether infection with given TB strain can occur while on IPT by DR compartment#
iota_r<- array(0, dim = length(DR_SET))

for (r in DR_SET){
  temp_param_df <- param_df%>%
    filter(notation == 'iota',
           DR_compartment == r)
  
  if(temp_param_df$calibrated == "yes"){
    temp_sim_df<-sim_calibration_ref_df%>%
      select(paste0("iota_",r,"."))
    iota_r[r] <- as.numeric((temp_sim_df[1]*
                              (temp_param_df$max[1]-temp_param_df$min[1]))+
                             temp_param_df$min[1])
  } else {
    iota_r[r] <- as.numeric(temp_param_df$value)
  }
}

#zeta - Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection#
temp_param_df <- param_df%>%
  filter(notation == 'zeta')
temp_sim_df<-sim_calibration_ref_df%>%
  select("zeta_")

zeta <-as.numeric((temp_sim_df[1]*
                     (temp_param_df$max[1]-temp_param_df$min[1]))+
                    temp_param_df$min[1])

#increased risk of re-infection for recovered populations#
temp_param_df <- param_df%>%
  filter(notation == 'xi')
temp_sim_df<-sim_calibration_ref_df%>%
  select("xi_")

xi <-as.numeric((temp_sim_df[1]*
                   (temp_param_df$max[1]-temp_param_df$min[1]))+
                  temp_param_df$min[1])

#########Parameters that describe TB progression ######

#omega - Rate of moving off of IPT, per year #
omega <-param_df%>%filter(notation == 'omega')
omega <- as.numeric(omega$value)

#pi_i_t - Base rates of TB progression#
pi_i_t <- array(data = 0, c(length(TB_SET), length(TB_SET)))

for (t_from in TB_SET){
  for (t_to in TB_SET){
    num_temp = (t_from*10) + t_to
    
    temp_param_df <- param_df%>%
      filter(notation == 'pi',
             TB_compartment == num_temp)
    
    if (nrow(temp_param_df) == 1){
      if(temp_param_df$calibrated == "yes"){
        temp_sim_df<-sim_calibration_ref_df%>%
          select(paste0("pi_",num_temp,"."))
        
        pi_i_t[t_from,t_to] <- as.numeric((temp_sim_df[1]*
                                             (temp_param_df$max[1]-temp_param_df$min[1]))+
                                            temp_param_df$min[1])
        
      } else{
        pi_i_t[t_from,t_to] <- as.numeric(temp_param_df$value)
      }
    }
  }
}

#theta_h - relative risk of TB progression#
theta_h <-array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp_param_df <- param_df%>%
    filter(notation == 'theta',
           HIV_compartment == h)
  
  if(temp_param_df$calibrated == "yes"){
    temp_sim_df<-sim_calibration_ref_df%>%
      select(paste0("theta_",h,"."))
    theta_h[h] <- as.numeric((temp_sim_df[1]*
                              (temp_param_df$max[1]-temp_param_df$min[1]))+
                             temp_param_df$min[1])
  } else {
    theta_h[h] <- as.numeric(temp_param_df$value)
  }
  
}

#upsilon_h - increased time infectious due to delayed treatment#
upsilon_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp_param_df <- param_df%>%
    filter(notation == 'upsilon',
           HIV_compartment == h)
  
  if(temp_param_df$calibrated == "yes"){
    temp_sim_df<-sim_calibration_ref_df%>%
      select(paste0("upsilon_",h,"."))
    upsilon_h[h] <- as.numeric((temp_sim_df[1]*
                                (temp_param_df$max[1]-temp_param_df$min[1]))+
                               temp_param_df$min[1])
  } else {
    upsilon_h[h] <- as.numeric(temp_param_df$value)
  }
  
}

#gamma_r -indicator if DR compartment can move onto after IPT#
gamma_r <- c(1,0)


#kappa_t_h_g_p - Rate of IPT initiation, per year#

#policy parameter not calibrated#
#time varying
kappa_t_h_g <- array(data = 0, c(length(TB_SET), length(HIV_SET), length(G_SET)))

ipt_initiation_function<-function(yr){
  
  kappa_t_h_g <- array(data = 0, c(length(TB_SET), length(HIV_SET), length(G_SET)))
  
  for (g in G_SET){
    
    IPT_init_temp<-ipt_initiation_df%>%
      filter(gender_id == g,
             year == yr,
             POLICY_ID == 1)
    
    kappa_t_h_g[c(1,3,4), 4, g]<-IPT_init_temp$ipt_init_perc
  }
  return(kappa_t_h_g)
}

#########Parameters that describe HIV progression ######
eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                            length(HIV_SET), 
                            length(G_SET)))

#hiv incidence factor for calibration
hiv_incidence_factor_g<-array(0, dim = length(G_SET))
for (g in G_SET){
  
  temp_sim_df<-sim_calibration_ref_df%>%
    select(paste0("eta.baselinefactor_",g))
  
  hiv_incidence_factor_g[g] <- temp_sim_df[[1]]
}

#changes 
HIV_transitions_param_func<-function(yr, N_t_r_h_g){
  
  eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                              length(HIV_SET), 
                              length(G_SET)))
  
  #hiv incidence
  #calibrated
  #time varying
  for (g in G_SET){
    hiv_incidence_temp<-hiv_incidence_df%>%
      filter(sex == if_else(g == 1, 'Male', 'Female'),
             year == yr)
    
    eta_i_h_g[1,2,g]<-as.numeric((hiv_incidence_factor_g[g]*
                                    (hiv_incidence_temp$max[1]-hiv_incidence_temp$min[1]))+
                                   hiv_incidence_temp$min[1])
  }
  
  #CD4 decline
  for (g in G_SET){
    temp_param_df <- param_df%>%
      filter(notation == 'eta',
             HIV_compartment == 23,
             G_compartment == g)
    
    temp_sim_df<-sim_calibration_ref_df%>%
      select(paste0("eta_23.",g))
    
    eta_i_h_g[2,3,g] <- as.numeric((temp_sim_df[1]*
                                      (temp_param_df$max[1]-temp_param_df$min[1]))+
                                     temp_param_df$min[1])
  }
  
  if(yr >= 2004){
    
    #get current model states (V_h_g in the appendix)
    #time varying based on model states
    hiv_prev_current<-sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4), g]])
    n2_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 2, g]])/hiv_prev_current
    n3_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 3, g]])/hiv_prev_current
    n4_prop_current_yr <- sum(N_t_r_h_g[N_t_r_h_g_ref[TB_SET, DR_SET, 4, g]])/hiv_prev_current
    
    #get percent eligible (time varying)
    prop_eligible_temp<-art_prop_eligible_df%>%
      filter(gender_id == g,
             year == yr)
    prop_eligible<-prop_eligible_temp$percent_hiv_compartment_2_eligible
    
    
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
#alpha_out - Rate of exit from the population#
alpha_out <- param_df%>%filter(notation == 'alpha.out')
alpha_out <- alpha_out$value

##placeholder for mu_t_h_g
mu_t_h_g<-array(0, dim=c(length(TB_SET), 
                               length(HIV_SET), 
                               length(G_SET)))

#non-disease mort
mu_baseline_factor_g<-array(0, dim = length(G_SET))
for (g in G_SET){

  temp_sim_df<-sim_calibration_ref_df%>%
    select(paste0("mu.baselinefactor_",g))
  
  mu_baseline_factor_g[g] <- temp_sim_df[[1]]
}

#HIV only increased mort
mu_risk_HIV_other_h<-array(0, dim = length(HIV_SET))
for (h in HIV_SET){
  
  #hiv neg is baseline mort
  if(h > 1){
    
    temp_sim_df<-sim_calibration_ref_df%>%
      select(paste0("risk.other_",h,"."))
    
    temp_param_df <- param_df%>%
      filter(notation == 'risk.other',
             HIV_compartment == h)
    
    mu_risk_HIV_other_h[h] <- as.numeric((temp_sim_df[1]*
                                            (temp_param_df$max[1]-temp_param_df$min[1]))+
                                           temp_param_df$min[1])
  } else {
    mu_risk_HIV_other_h[h] <-1
  }
    
}

#TB/HIV increased mort
mu_risk_HIV_TB_h<-array(0, dim = length(HIV_SET))
for (h in HIV_SET){
  
  temp_sim_df<-sim_calibration_ref_df%>%
    select(paste0("risk.TB_",h,"."))
  
  temp_param_df <- param_df%>%
    filter(notation == 'risk.TB',
           HIV_compartment == h)
  
  mu_risk_HIV_TB_h[h] <- as.numeric((temp_sim_df[1]*
                                          (temp_param_df$max[1]-temp_param_df$min[1]))+
                                         temp_param_df$min[1])
}


#mu_t_h_g calcs- mortality rates#
mort_param_func <-function(yr){
  
  mu_t_h_g <- array(0, dim = c(length(TB_SET),
                               length(HIV_SET), 
                               length(G_SET)))
  
  for (g in G_SET){
    
    if(yr < 1990){
      mort_base_temp <- mort_rate_df%>%
        filter(year == 1990,
               sex == if_else(g == 1, 'Male', 'Female'))
      
      mort_base_temp <- as.numeric((mu_baseline_factor_g[g]*
                                      (mort_base_temp$max[1]-mort_base_temp$min[1]))+
                                     mort_base_temp$min[1])
      
    } else if (yr < 2018) {
      mort_base_temp <- mort_rate_df%>%
        filter(year == yr,
               sex == if_else(g == 1, 'Male', 'Female'))
      
      mort_base_temp <- as.numeric((mu_baseline_factor_g[g]*
                                      (mort_base_temp$max[1]-mort_base_temp$min[1]))+
                                     mort_base_temp$min[1])
    } else {
      mort_base_temp <- mort_rate_df%>%
        filter(year == 2018,
               sex == if_else(g == 1, 'Male', 'Female'))
      
      mort_base_temp <- as.numeric((mu_baseline_factor_g[g]*
                                      (mort_base_temp$max[1]-mort_base_temp$min[1]))+
                                     mort_base_temp$min[1])
    }
    
    for (h in HIV_SET){
      #first update hiv only values
      for (t in c(1:5, 7:8)){
        if (h > 1){
          mu_t_h_g[t,h,g]<-mu_risk_HIV_other_h[h]*mort_base_temp
        } else {
          #no hiv no tb
          mu_t_h_g[t,1,g]<-mort_base_temp
        }
      }
      #tb and/or hiv values
      mu_t_h_g[6,h,g]<-mu_risk_HIV_TB_h[h]*mort_base_temp
    }
  }
  return(mu_t_h_g)
}

#alpha_in_t_h_g - Proportion of population that enters each compartment each year
alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                      length(DR_SET), 
                                      length(HIV_SET), 
                                      length(G_SET)))

aging_in_param_func <-function(yr){
  alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                        length(DR_SET), 
                                        length(HIV_SET), 
                                        length(G_SET)))
  if(yr < 1990){
    birth_perc_df_temp<-birth_perc_df%>%
    filter(year == 1989)
  } else if (yr > 2018){
    birth_perc_df_temp<-birth_perc_df%>%
      filter(year == 2017)
  } else {
    birth_perc_df_temp<-birth_perc_df%>%
      filter(year == yr)
    }
  
  for (t in TB_SET){
    for (r in DR_SET){
      for (h in HIV_SET){
        for (g in G_SET){
          birth_perc_df_temp2<-birth_perc_df_temp%>%
            filter(TB_compartment == t,
                   DR_compartment == r,
                   HIV_compartment == h,
                   G_compartment == g)
          
          alpha_in_t_r_h_g[t,r,h,g] <- birth_perc_df_temp2$prop_of_pop
        }
      }
    }
  }
  return(alpha_in_t_r_h_g)
}

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
  
  #initiation delta array, extra for tracking model outputs
  dN_t_r_h_g <- array(0, dim = 139) 
  
  #####time varying parameters#####
  #Force of Infection Calculations#
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
  print(current_yr)
  
  #HIV transitions start in 1980
  if (current_yr > 1980){
    eta_i_h_g<<-HIV_transitions_param_func(current_yr, N_t_r_h_g)
  }
  
  #parameters that change yearly
  if(current_yr > last_year){
    
    #deaths
    mu_t_h_g<<-mort_param_func(current_yr)
    total_out_t_r_h_g<<-total_out_param_func(mu_t_h_g)
    
    #entries prop
    alpha_in_t_r_h_g<<-aging_in_param_func(current_yr)
    
    #only start IPT initiations in 2005
    if(current_yr > 2004){
      kappa_t_h_g<<-ipt_initiation_function(current_yr)
    }
    last_year<<-current_yr
  }
  
  #####Model outputs#####

  ##health metric calculations
  #hiv prev calculated after dsolve complete
  
  Tb_inc_neg_male_loc<-129
  Tb_inc_neg_female_loc<-130
  Tb_inc_pos_male_loc<-131
  Tb_inc_pos_female_loc<-132
  
  Tb_mort_neg_male_loc<-133
  Tb_mort_neg_female_loc<-134
  Tb_mort_pos_male_loc<-135
  Tb_mort_pos_female_loc<-136
  
  O_mort_male_loc<-137
  O_mort_female_loc<-138
  
  #treatment rates
  #ART calculated after desolve complete
  IPT_loc<-139
  
  #TB incidence metrics
  #HIV-, males (TBinc_1^HIV+)
  dN_t_r_h_g[Tb_inc_neg_male_loc]<-dN_t_r_h_g[129]<- (theta_h[1]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,1,1]])) + #from recent TB infection to active 
    (theta_h[1]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,1,1]])) + #from remote TB infection to active 
    (theta_h[1]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,1,1]])) + #from TB infection on IPT to active
    (theta_h[1]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,1,1]])) #from TB after on IPT to active
  
  #HIV- females
  dN_t_r_h_g[Tb_inc_neg_female_loc]<- (theta_h[1]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,1,2]])) + #from recent TB infection to active 
    (theta_h[1]*pi_i_t[4,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[4,1:2,1,2]])) + #from remote TB infection to active 
    (theta_h[1]*pi_i_t[5,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[5,1:2,1,2]])) + #from TB infection on IPT to active
    (theta_h[1]*pi_i_t[8,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[8,1:2,1,2]])) #from TB after on IPT to active
  
  #HIV+, males
  #CD4 more
  dN_t_r_h_g[Tb_inc_pos_male_loc]<- (theta_h[2]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,2,1]])) + #from recent TB infection to active 
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
  dN_t_r_h_g[Tb_inc_pos_female_loc]<-(theta_h[2]*pi_i_t[3,6]*sum(N_t_r_h_g[N_t_r_h_g_ref[3,1:2,2,2]])) + #from recent TB infection to active 
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
  
  #track TB mort
  #tb active, hiv -, males
  dN_t_r_h_g[Tb_mort_neg_male_loc]<-mu_t_h_g[6,1,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,1,1]]) 
  #tb active, hiv -, females
  dN_t_r_h_g[Tb_mort_neg_female_loc]<-mu_t_h_g[6,1,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,1,2]]) 
  #tb active, hiv +, males
  dN_t_r_h_g[Tb_mort_pos_male_loc]<-mu_t_h_g[6,2,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,2,1]])+ #tb active, hiv+, CD4 > 200, males
    mu_t_h_g[6,3,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,3,1]])+ #tb active, hiv+, CD4 <= 200, males
    mu_t_h_g[6,4,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,4,1]]) #tb active, hiv+, on ART, males
  #tb active, hiv +, females
  dN_t_r_h_g[Tb_mort_pos_female_loc]<-mu_t_h_g[6,2,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,2,2]])+ #tb active, hiv+, CD4 > 200, females
    mu_t_h_g[6,3,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,3,2]])+ #tb active, hiv+, CD4 <= 200, females
    mu_t_h_g[6,4,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[6,1:2,4,2]]) #tb active, hiv+, on ART, females
  
  #track HIV other mort
  #base mortality is same for TB compartments no TB
  #so look up from TB compartment 1 by HIV status and gender since by TB (inactive) are all the same
  
  dN_t_r_h_g[O_mort_male_loc]<-mu_t_h_g[1,2,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1,2,3,4,5,7,8),1:2,2,1]])+ #NO tb active, hiv+, CD4 > 200, males
    mu_t_h_g[1,3,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1,2,3,4,5,7,8),1:2,3,1]])+ #NO tb, hiv+, CD4 <= 200, males
    mu_t_h_g[1,4,1]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1,2,3,4,5,7,8),1:2,4,1]]) #NO tb, hiv+, on ART, males
  
  dN_t_r_h_g[O_mort_female_loc]<-mu_t_h_g[1,2,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1,2,3,4,5,7,8),1:2,2,2]])+ #no tb, hiv+, CD4 > 200, females
    mu_t_h_g[1,3,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1,2,3,4,5,7,8),1:2,3,2]])+ #no tb, hiv+, CD4 <= 200, females
    mu_t_h_g[1,4,2]*sum(N_t_r_h_g[N_t_r_h_g_ref[c(1,2,3,4,5,7,8),1:2,4,2]]) #no tb, hiv+, on ART, females
  
  
  #Treatment Rates
  dN_t_r_h_g[IPT_loc]<-kappa_t_h_g[1, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[1, DR_SET, 4, 1]])+ #males from uninfected
    kappa_t_h_g[1, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[1, DR_SET, 4, 2]])+ #females from uninfected
    kappa_t_h_g[3, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[3, DR_SET, 4, 1]])+ #males from recent
    kappa_t_h_g[3, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[3, DR_SET, 4, 2]])+#females from recent
    kappa_t_h_g[4, 4, 1]*sum(N_t_r_h_g[N_t_r_h_g_ref[4, DR_SET, 4, 1]])+ #males from remote
    kappa_t_h_g[4, 4, 2]*sum(N_t_r_h_g[N_t_r_h_g_ref[4, DR_SET, 4, 2]]) #males from remote
  
  #entries and exits from the population
  B <- sum(total_out_t_r_h_g*N_t_r_h_g[1:128])
  
  #######TB compartment 1 Equations #########
  #Set DR compartment to 1, since not applicable to drug resistant compartments#
  for (h in HIV_SET){
    for (g in G_SET){
      dN_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]<-(
        sum(alpha_in_t_r_h_g[1,1,h,g]*B) - #entries from births
          (total_out_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) + #exists from aging out and death
          (omega*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) - #entries from off IPT
          (kappa_t_h_g[1,h,g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]]) - #exists from on to IPT
          (FOI*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])+ #exists from TB infection 
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
        -(total_out_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])+ #exists from aging out and death
          (kappa_t_h_g[1,h,g]*N_t_r_h_g[N_t_r_h_g_ref[1,1,h,g]])- #entries from on to IPT
          (omega*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]])- #exits from off IPT
          ((sum(iota_r*FOI_r))*N_t_r_h_g[N_t_r_h_g_ref[2,1,h,g]]) + #exits from infection (diminished for IPT)
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

########Feed paramters into desolve########
setwd(outdir)

#####N_init - total initial population in each compartment flattened in 1D array for ODE solver #####
#arrange pop init df, TB-->DR-->HIV-->G
pop_init_df<-pop_init_df%>%
  arrange(G_compartment)%>%
  arrange(HIV_compartment)%>%
  arrange(DR_compartment)%>%
  arrange(TB_compartment)

N_init <- pop_init_df$total_pop

#add in mort calc placeholders
model_output_names<-c('Tb_inc_neg_male',
                      'Tb_inc_neg_female',
                      'Tb_inc_pos_male',
                      'Tb_inc_pos_female',
                      'Tb_mort_neg_male',
                      'Tb_mort_neg_female',
                      'Tb_mort_pos_male',
                      'Tb_mort_pos_female',
                      'O_mort_male',
                      'O_mort_female',
                      'total_IPT_init')


N_init <-c(N_init, rep(0, times = length(model_output_names)))

names(N_init) <- c(pop_init_df$compartment_id, 
                   model_output_names)

out_df<-as.data.frame(ode(times = TT_SET, y = N_init,
                          func = tb_hiv_prog_calibration_model, method = 'lsoda',
                          parms = NULL))

out_df$hiv_prev_male<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),1])])
out_df$hiv_prev_female<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, c(2,3,4),2])])
out_df$ART_coverage<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, 4,G_SET])])

out_df$total_male_pop<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, HIV_SET,1])])
out_df$total_female_pop<-rowSums(out_df[,1+c(N_t_r_h_g_ref[TB_SET, DR_SET, HIV_SET,1])])

out_df<- cbind(year = as.integer(start_yr+out_df$time), sim_id = rep(sim_id_current_eval,
                                                                     times = nrow(out_df)),
               out_df)

write.csv(out_df, file = paste0('out_df_sim_id_', sim_id_current_eval, '.csv'), row.names = FALSE)
