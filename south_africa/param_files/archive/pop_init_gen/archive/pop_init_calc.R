#pop init warmup period 1900-1990
#1900-1980 (Gender, TB, and DR only)
#1980 introduce HIV

#mort params set to 1990 estimated values
#birth rate params set to 1990 values 
#estimated pop init (before warm up period) set at 1990 values
#birth rate params and estimated pop init are...
#summed over HIV compartments and assigned to HIV-, before 1980

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'readr'), require, character.only=T)

#Make sure you have the epi_model_HIV_TB.Rproj open (top right)
current_date <- gsub('-', '_', Sys.Date())
outdir <- paste0(here(),'/model_outputs/calibration/', current_date, '/ref_data')
indir <- paste0(here(),'/param_files')
dir.create(file.path(outdir))

#######Create calibration ref df#########
#mort calibration test
#(1 = low, 2 = high)
#TB_only, TB_HIV_CD4More, TB_HIV_CD4Less, TB_HIV_ART

#Test low/High
# mort_calib_TB_only_test <- c('L', 'H')
# mort_calib_TB_HIV_CD4More_test <- c('L', 'H')
# mort_calib_TB_HIV_CD4Less_test <- c('L', 'H')
# mort_calib_TB_HIV_ART_test <- c('L', 'H')

#Test Low/Medium/High
mort_calib_TB_only_test <- c('L', 'M', 'H')
mort_calib_TB_HIV_CD4More_test <- c('L', 'M', 'H')
mort_calib_TB_HIV_CD4Less_test <- c('L', 'M', 'H')
mort_calib_TB_HIV_ART_test <- c('L', 'M', 'H')

#Test High
# mort_calib_TB_only_test <- c('H')
# mort_calib_TB_HIV_CD4More_test <- c('H')
# mort_calib_TB_HIV_CD4Less_test <- c('H')
# mort_calib_TB_HIV_ART_test <- c('H')

#betas to test
beta_1_test<-seq(from = 6, to = 10, by = 0.5)
beta_2_test<-seq(from = 6, to = 10, by = 0.5)

#create a sim id for each combination
sim_calib_id<-c()
b1_calib_id<-c()
b2_calib_id<-c()
m1_calib_id<-c()
m2_calib_id<-c()
m3_calib_id<-c()
m4_calib_id<-c()

sim_id_temp<-1

for (b1 in beta_1_test){
  for (b2 in beta_2_test){
    for (m1 in mort_calib_TB_only_test){
      for (m2 in mort_calib_TB_HIV_CD4More_test){
        for (m3 in mort_calib_TB_HIV_CD4Less_test){
          for (m4 in mort_calib_TB_HIV_ART_test){
            sim_calib_id<-c(sim_calib_id, sim_id_temp)
            b1_calib_id<-c(b1_calib_id, b1)
            b2_calib_id<-c(b2_calib_id, b2)
            m1_calib_id<-c(m1_calib_id, paste0('1-', m1)) #to id proper
            m2_calib_id<-c(m2_calib_id, paste0('2-', m2)) #calibration params
            m3_calib_id<-c(m3_calib_id, paste0('3-', m3))
            m4_calib_id<-c(m4_calib_id, paste0('4-', m4))
            sim_id_temp<-sim_id_temp+1
          }
        }
      }
    }
  }
}

sim_calibration_ref_df<-data.frame(sim_calib_id, b1_calib_id, b2_calib_id,
                                   m1_calib_id, m2_calib_id, m3_calib_id, m4_calib_id)

setwd(outdir)
write.csv(sim_calibration_ref_df, 'sim_calibration_ref_df.csv', row.names = FALSE)

#remove everything but df
rm(list=setdiff(ls(), c("sim_calibration_ref_df",
                "indir",
                "outdir")))

#read in data
setwd(indir)
param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
mort_df <- read.csv('mort_df.csv')
#hiv_transition_df is different from calib code
hiv_transition_df<-read_excel("hiv_param_gen/hiv_input_gen_data.xlsx") 
birth_rate_df<-read.csv('birth_rate_df.csv')

#read GBD pop_estimates
setwd(paste0(here(), '/param_files/GBD_pop_estimates'))
pop_df<-read.csv('pop_df.csv')

#clean dataframe column names for consistency
names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))

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

#pop init 1990
#establish compartment and dcompartment ids
#TB factors 
prop_recent<-0.049  #proprotion of LTBI that is recent
prop_remote<-0.951 #proportion remote
prop_IPT<-0#.01
prop_LTBI<-.49
LTBI_recent<-(prop_LTBI*prop_recent)*(1-prop_IPT)
LTBI_remote<-(prop_LTBI*prop_remote)*(1-prop_IPT)
LTBI_on_IPT<-(prop_LTBI*prop_IPT)
TB_active_male<-.01
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
g_adj<-c()

#calculate gender init adjustments from GBD data
pop_df<-pop_df%>%
  group_by(year)%>%
  mutate(total_pop = sum(expected_total_pop))%>%
  ungroup()%>%
  mutate(percent_gender_pop = expected_total_pop/total_pop)%>%
  group_by(sex_id)%>%
  summarise(gender_percent = mean(percent_gender_pop))%>%
  arrange(sex_id)

gender_adjustments <- pop_df$gender_percent

for (tb in TB_SET){
  for (dr in DR_SET){
    for (g in G_SET){
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
      
      #gender adjustments
      g_adj<-c(g_adj, gender_adjustments[g])
    }
  }
}


pop_init_df<-data.frame(TB_compartment, DR_compartment, G_compartment, 
                        tb_adj, dr_adj, g_adj)

rm(tb_adj, dr_adj, g_adj, TB_compartment, DR_compartment, G_compartment)

pop_init_df<-pop_init_df%>%
  mutate(total_adj = round(tb_adj*dr_adj*g_adj,4),
         prop_of_pop = round(total_adj/sum(total_adj),5))%>%
  mutate(prop_of_pop = prop_of_pop*(1/sum(prop_of_pop)))%>% #fix rounding errors make sure add to 1
  mutate(total_pop = prop_of_pop*100000)%>%
  select(c('TB_compartment', 'DR_compartment', 'G_compartment', 'total_pop'))

#set the other HIV compartments to 0
pop_init_df_tempHIV<-rbind(pop_init_df, pop_init_df, pop_init_df)
pop_init_df_tempHIV$HIV_compartment <- rep(2:4, each = nrow(pop_init_df))
pop_init_df_tempHIV$total_pop<-rep(0, times = nrow(pop_init_df_tempHIV)) 

#combine dfs
pop_init_df$HIV_compartment<-rep(1, times = nrow(pop_init_df))
pop_init_df<-rbind(pop_init_df, pop_init_df_tempHIV)

pop_init_df<- pop_init_df%>%
  mutate(compartment_id = paste0("N_", TB_compartment, 
                                 "_", DR_compartment ,
                                 "_", HIV_compartment, 
                                 "_", G_compartment),
         dcompartment_id = paste0("dN_", TB_compartment,
                                  "_", DR_compartment ,"_",
                                  HIV_compartment, "_",
                                  G_compartment))

pop_init_df<-pop_init_df%>%
  arrange(G_compartment)%>%
  arrange(HIV_compartment)%>%
  arrange(DR_compartment)%>%
  arrange(TB_compartment)

N_init <- pop_init_df$total_pop
names(N_init)<-pop_init_df$compartment_id

######## PARAMETERS THAT IMPACT FORCE OF INFECTION #######

######### beta_g - Number of effective contacts for TB transmission per infectious year######

#calibrated parameter in sim_calibration_ref_df

###### phi_h - Relative transmissibility of TB in HIV pops#########
phi_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'phi',
           HIV_compartment == h)
  phi_h[h] <- temp$Reference_expected_value
}

### upsilon_h - increased time infectious due to delayed treatment###
upsilon_h <- array(0, dim = length(HIV_SET))

for (h in HIV_SET){
  temp <- param_df%>%
    filter(notation == 'upsilon',
           HIV_compartment == h)
  upsilon_h[h] <- temp$Reference_expected_value
}

####increased risk of re-infection for recovered populations
xi<-param_df%>%
  filter(notation == 'xi')
xi <-xi$Reference_expected_value

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

###########gamma_r -indicator if DR compartment can move onto after IPT####
gamma_r <- c(1,0)

#######Parameters that describe HIV progression########
#####eta_i_h_g rate HIV transitions ######
eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                            length(HIV_SET), 
                            length(G_SET)))

HIV_transitions_param_func<-function(yr){
  eta_i_h_g <- array(0, dim=c(length(HIV_SET), 
                              length(HIV_SET), 
                              length(G_SET)))
  
  if (yr < 1980){
    return(eta_i_h_g)
  } else {
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
      eta_i_h_g[1,2,g]<-temp2$hiv_incidence #NO ART in warmup period
      }
    return(eta_i_h_g)
    }
  
}

#########parameters for death and aging rates ###########

######mu_t_h_g - mortality rates ########
#set at 1990 mortality rates#
#testing two levels
#N_MORT_LEVELS = 2
#MORT_LEVELS<-c('L', 'H')
N_MORT_LEVELS = 3
MORT_LEVELS<-c('L', 'M', 'H')

mu_t_h_g_LEVEL <- array(0, dim = c(length(TB_SET), 
                                   length(HIV_SET), 
                                   length(G_SET),
                                   N_MORT_LEVELS))

mu_t_h_g <- array(0, dim = c(length(TB_SET), 
                                   length(HIV_SET), 
                                   length(G_SET)))

for (l in MORT_LEVELS){
  for (t in TB_SET){
    for (h in HIV_SET){
      for (g in G_SET){
        temp <- mort_df%>%
          filter(year == 1990,
                 TB_compartment == t,
                 HIV_compartment == h,
                 G_compartment == g,
                 level == l)
        
        #level_temp<-if_else(l == 'L', 1, 2)
        level_temp<-1
        
        mu_t_h_g_LEVEL[t,h,g,level_temp] <- unique(temp$mort_rate)
      }
    }
  }
}

########alpha_in_t_h_g - Proportion of population that enters each compartment#####
birth_rate_param_func<-function(yr){
  alpha_in_t_r_h_g <- array(data = 0, c(length(TB_SET), 
                                        length(DR_SET), 
                                        length(HIV_SET), 
                                        length(G_SET)))
    
    for (t in TB_SET){
      for (r in DR_SET){
        for (h in HIV_SET){
          for (g in G_SET){
            if (yr > 1980){
            temp <- (unique(birth_rate_df$prop_of_pop[(birth_rate_df$year == 1990) &
                                                        (birth_rate_df$TB_compartment == t)&
                                                        (birth_rate_df$DR_compartment == r)&
                                                        (birth_rate_df$HIV_compartment == h)&
                                                        (birth_rate_df$G_compartment == g)]))
            
            
            alpha_in_t_r_h_g[t,r,h,g] <- temp
            
            } else {
            temp <- birth_rate_df%>%
              filter(year == 1990,
                     TB_compartment == t,
                     DR_compartment == r,
                     G_compartment == g)
            alpha_in_t_r_h_g[t,r,1,g]<-sum(temp$prop_of_pop) #sum over HIV compartments
          }
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



#add in mort calc placeholders
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
pop_init_HIV_TB_model <- function(time, N_t_r_h_g, parms){
  
  #initiation delta array, +4 for tracking mort
  dN_t_r_h_g <- array(0, dim = length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET))
  
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
  
  #HIV transitions
  if (current_yr > last_year){
    eta_i_h_g<<-HIV_transitions_param_func(current_yr)
    
    #alpha in proportions
    alpha_in_t_r_h_g<<-birth_rate_param_func(current_yr)
    
    #deaths/total out (set at 1990 rates)
    #mu_t_h_g<-mort_param_func(current_yr)
    #total_out_t_r_h_g<-total_out_param_func(mu_t_h_g)
    
    last_year<<-current_yr
  }
  
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
start_yr = 1900
end_yr = 1990
TT<-end_yr-start_yr
time_interval <- 1/12
TT_SET <- seq(from = 0, to = TT, by = time_interval)

#feed in to solve, evaluation time intervals (TT_SET), 
#initial population (N_init)
#differential equations and dynamic parameter (open_seir_model)
#specify differential solver (lsoda)
# read params from global (NULL)

pop_init_df_out_all<-data.frame()
sim_id <-1
last_year<-start_yr-1

for(n in 1:nrow(sim_calibration_ref_df)){
    
    last_year<-start_yr-1
    sim_id <-n
    print(sim_id)
    #to track how long it takes to solve
    start<-Sys.time()
    sim_attributes<-sim_calibration_ref_df[n, ]
    sim_id<-as.integer(sim_attributes[1])
    beta_g<<-as.double(sim_attributes[2:3])
    mort_rate_levels_test<-sim_attributes[4:7]
    #mort_rate_levels_num<-lapply(mort_rate_levels_test, function(x) if_else(grepl('L', x), 1, 2))
    mort_rate_levels_num<-lapply(mort_rate_levels_test, function(x) if_else(grepl('L', x), 1, 1))
    
    mu_t_h_g[6,1,G_SET]<-mu_t_h_g_LEVEL[6, 1, G_SET, unlist(mort_rate_levels_num[1])]
    mu_t_h_g[6,2,G_SET]<-mu_t_h_g_LEVEL[6, 2, G_SET, unlist(mort_rate_levels_num[2])]
    mu_t_h_g[6,3,G_SET]<-mu_t_h_g_LEVEL[6, 3, G_SET, unlist(mort_rate_levels_num[3])]
    mu_t_h_g[6,4,G_SET]<-mu_t_h_g_LEVEL[6, 4, G_SET, unlist(mort_rate_levels_num[4])]
    
    out_df<-as.data.frame(ode(times = TT_SET, y = N_init, 
                              func = pop_init_HIV_TB_model, method = 'lsoda',
                              parms = NULL))
    
    out_df<- cbind(sim_id = rep(sim_id, times = nrow(out_df)), 
                   out_df)
    
    pop_init_df<-out_df%>%
      filter(time == 90)%>%
      select(-c('time'))
    
    pop_init_df_out <-reshape2::melt(pop_init_df, id.vars = c('sim_id'))
    pop_init_df_out <- cbind(pop_init_df_out,
                             data.frame(do.call('rbind',
                                                strsplit(as.character(pop_init_df_out$variable),
                                                         '_',fixed=TRUE))))
    names(pop_init_df_out)[names(pop_init_df_out) == "X2"] <- "TB_compartment"
    names(pop_init_df_out)[names(pop_init_df_out) == "X3"] <- "DR_compartment"
    names(pop_init_df_out)[names(pop_init_df_out) == "X4"] <- "HIV_compartment"
    names(pop_init_df_out)[names(pop_init_df_out) == "X5"] <- "G_compartment"
    pop_init_df_out<-pop_init_df_out%>%select(-c('X1'))
    
    pop_init_df_out_all<-rbind(pop_init_df_out_all, pop_init_df_out)
    
    print(Sys.time()-start)
    
}

setwd(outdir)
write.csv(pop_init_df_out_all, 'pop_init_df.csv', row.names = FALSE)

