##Ref HIV_TB_Model_Appendix_Jun3.pdf##
#Last Modified: Jun 17

#clean workspace
rm(list = ls())
gc()

#Import Libraries#
library(dplyr)
library(deSolve)
library(readxl)
library(stringr)

#Define input directory
indir<-("~/Google Drive/HIV:TB/param_file")

#Define output directory - Set to JR files
outdir<-("~/Google Drive/HIV:TB/test_outputs/Jun17")

set.seed(1)

#read parameter file
setwd(indir)
param_df <- read_excel("Epi parameters_June_2020+fakedata.xlsx")

#remove all spaces from column names
names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))

#####make sure all compartments are integer type for proper indexing#####
param_df$TB_compartment<-as.integer(param_df$TB_compartment)
param_df$DR_compartment<-as.integer(param_df$DR_compartment)
param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
param_df$G_compartment<-as.integer(param_df$G_compartment)
param_df$P_compartment<-as.integer(param_df$P_compartment)

#Set of Scenarios (S)
#S_SET<-1:10

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

#TB subsets - groups of compartments
TB_SUBSET_UNINFECTED<-1:2
TB_SUBSET_IPT<-c(2,5)
TB_SUBSET_LTBI<-3:5
TB_SUBSET_IPT_INIT_FROM<-c(1,3,4)

#TB subsets - individual compartments
TB_SUBSET_UNINFECTED_NOIPT<-1
TB_SUBSET_UNINFECTED_IPT<-2
TB_SUBSET_INFECTED_RECENT<-3
TB_SUBSET_INFECTED_REMOTE<-4
TB_SUBSET_INFECTED_IPT<-5
TB_SUBSET_ACTIVE<-6
TB_SUBSET_RECOVERED<-7
TB_SUBSET_AFTERIPT<-8

#DR - Set of TB Drug Resistance States (DR)
#1: Drug-susceptible (DS)
#2: Multidrug resistant (MDR)

DR_SET<-1:2

#TB subsets - individual compartments
DR_SUBSET_DS<-1
DR_SUBSET_MDR<-2

#4 HIV compartments (HIV)#
#1 : HIV Negative
#2 : HIV Positive CD4 > 200 - No ART
#3 : HIV Positive CD4 =<: 200 - No Art
#4 : HIV Positive - ART 

HIV_SET<-1:4

#HIV subset groups
HIV_SUBSET_POS_NOART<-2:3

#HIV subsets - individual compartments
HIV_SUBSET_NEG <- 1
HIV_SUBSET_POS_CD4LESS <-2
HIV_SUBSET_POS_CD4MORE <-3
HIV_SUBSET_POS_ART <- 4

#Genders (G): Male : 1, Female : 2
G_SET<-1:2

#Gender subsets - individual compartments
G_SUBSET_M<-1
G_SUBSET_F<-2

#Set of policies (P): 
#1 : Standard of care only, no ART or IPT delivery (base)
#2: Standard of care and community-based ART delivery only, no IPT delivery
#3: Standard of care and community-based ART and IPT delivery
P_SET<-1:3

P_SUBSET_BASE <-1
P_SUBSET_NO_IPT_DELIVERY <-2
P_SUBSET_IPT_DELIVERY <- 3

#Time Horizon (5 years)
TT<-5
time_interval <- 1/12
TT_SET <- c(1:(12*TT))*time_interval

######create compartment ids and reference matrix to go between multidementional arrays and 1D array####
n_compartments <- length(TB_SET)*length(DR_SET)*length(HIV_SET)*length(G_SET)
comparment_id<-rep(NA, n_compartments)
n = 0

lapply(TB_SET, function(t){
  lapply(DR_SET, function(r){
    lapply(HIV_SET, function(h){
      lapply(G_SET, function(g){
        n <<- n+1
        comparment_id[n] <<- paste0("N_", t, "_", r, "_", h, "_", g)
      })
    })
  })
})

N_t_r_h_g_ref <- seq(1: n_compartments)
N_t_r_h_g_ref <-array(data=N_t_r_h_g_ref, dim = c(length(TB_SET),
                                                  length(DR_SET),
                                                  length(HIV_SET),
                                                  length(G_SET)
))

#create delta compartment ids
d_comparment_id<-paste0('d', comparment_id)

##########POPULATION##########

#set initial populations
N_t_r_h_g_init <-array(data = 0, dim = c(length(TB_SET),
                                          length(DR_SET),
                                          length(HIV_SET),
                                          length(G_SET)))

#set init populations here at TT = 1
N_t_r_h_g_init[TB_SET, DR_SET, HIV_SET, G_SET] = 10000 # will change to pop init when determined

#make into 1D array for ODE
compartments_init <- c(N_t_r_h_g_init) #under all policies the initial pop is the same so can just set init population using any aribritrary policy
names(compartments_init) <- comparment_id

#create an matrix for recording FOI over time
lambda_r_g_tau_p <-array(data=0, dim = c(length(DR_SET),
                                         length(G_SET),
                                         length(TT_SET),
                                         length(P_SET)
))

#############MODEL PARAMETERS#############

#########Parameters that impact force of infection#######

#Number of effective contact for TB transmission per infectious year for gender g
beta_g <- array(data = 0, dim = length(G_SET))

#filter dataframe for beta params
beta_params <- param_df%>%
  filter(notation == 'beta')

lapply(G_SET, function (g){
  temp <- beta_params%>%filter(G_compartment == g)
  beta_g[g] <<- temp$Reference_expected_value
})

#rm dataframe for beta params
rm(beta_params)

#Relative transmissibility of TB in populations living in HIV compartment h
phi_h <-array(data = 0, dim = length(HIV_SET))

#filter dataframe for phi params
phi_params <- param_df%>%
  filter(notation == 'phi')

lapply(HIV_SET, function(h){
  temp <- phi_params%>%filter(HIV_compartment == h)
  phi_h[h] <<- temp$Reference_expected_value
})

rm(phi_params)

#Fraction of new TB infections that are MDR-TB
varepsilon_g <-array(data = 0, dim = length(G_SET))

#filter dataframe for varepsilon params
varepsilon_params <- param_df%>%
  filter(notation == 'varepsilon')

lapply(G_SET, function(g){
  temp <- varepsilon_params%>%filter(G_compartment == g)
  varepsilon_g[g] <<- temp$Reference_expected_value
})

rm(varepsilon_params)

#Indicator for whether infection with given TB strain can occur while on IPT for populations in DR compartment r
iota_r <-array(data = 0, dim = length(DR_SET))

#filter dataframe for iota params
iota_params <- param_df%>%
  filter(notation == 'iota')

lapply(DR_SET, function(r){
  temp <- iota_params%>%filter(DR_compartment == r)
  iota_r[r] <<- temp$Reference_expected_value
})

rm(iota_params)

#Indicator the diminished force of infection due to partially-protective effects of IPT after moving off of IPT for populations with LTBI
#filter dataframe for upsilon params
upsilon_params <- param_df%>%
  filter(notation == 'upsilon')

upsilon <- upsilon_params$Reference_expected_value

rm(upsilon_params)

#Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection
#filter dataframe for zeta params
zeta_params <- param_df%>%
  filter(notation == 'zeta')

zeta <- zeta_params$Reference_expected_value

rm(zeta_params)

###########Parameters that describe TB progression##############
#Rate of IPT initiation from TB compartment t and HIV compartment h for gender g under policy p, per year
#0 where not applicable

kappa_t_h_g_p <- array(data = 0, dim = c(length(TB_SET),
                                         length(HIV_SET),
                                         length(G_SET),
                                         length(P_SET)))

#filter dataframe for kappa params
kappa_params <- param_df%>%
  filter(notation == 'kappa')

lapply(TB_SUBSET_IPT_INIT_FROM, function(t){
  lapply(HIV_SET, function(h){
    lapply(G_SET, function(g){
      lapply(P_SET, function(p){
        temp <- kappa_params%>%
          filter(TB_compartment == t,
                 HIV_compartment == h,
                 G_compartment == g,
                 P_compartment == p)
        
        kappa_t_h_g_p[t,h,g,p] <<- temp$Reference_expected_value
      })
    })
  })
})

rm(kappa_params)

#Rate rate of moving off of IPT from TB compartment t under policy p, per year
#0 where not applicable
omega_t_p<- array(data = 0, dim = c(length(TB_SET),
                                    length(P_SET)
                                    ))

#filter dataframe for omega params
omega_params <- param_df%>%
  filter(notation == 'omega')

lapply(TB_SUBSET_IPT, function(t){
  lapply(P_SET, function(p){
    temp <- omega_params%>%
          filter(TB_compartment == t,
                 P_compartment == p)
    
    omega_t_p[t,p] <<- temp$Reference_expected_value
    
  })
})

rm(omega_params)

#Base rates of TB progression of infected populations from TB compartment i to TB compartment t, per year (i, t) ∈ T B (set to zero where not applicable and not included in equations, as shown in Figure 2 and Section 2 respectively)

#filter dataframe for pi params
pi_params <- param_df%>%
  filter(notation == 'pi')

pi_i_t <- array(data = 0, dim = c(length(TB_SET),
                            length(TB_SET)
                            ))

lapply(1:nrow(pi_params), function(x){
  i <- as.integer(pi_params$TB_compartment[x]/10)
  t <- pi_params$TB_compartment[x]%%10
  #pi_i_t[i,t] <<- pi_params$Reference_expected_value[x]
  pi_i_t[i,t] <<- .05
})

rm(pi_params)

#elative risk for TB progression from LTBI to Active for HIV compartment h
theta_h <- array(data = 0, dim = length(HIV_SET))

#filter dataframe for theta params
theta_params <- param_df%>%
  filter(notation == 'theta')

lapply(HIV_SET, function(h){
  temp<-theta_params%>%filter(HIV_compartment == h)
  theta_h[h] <<-temp$Reference_expected_value
})

rm(theta_params)

#1 if in drug-susceptible, DR compartment r ∈ 1 ⊂ DR, 0 if in MDR-TB, DR com- partment r ∈ 2 ⊂ DR to indicate that populations with MDR-TB cannot move into LTBI after IPT
gamma_r <- c(0,1)

##########Parameters that describe HIV progression########
#Rate of populations moving from HIV compartment i to HIV compartment h, per year, (i, h) ∈ HIV under policy p ∈ P #(set to zero where not applicable)

#filter dataframe for eta params
eta_params <- param_df%>%
  filter(notation == 'eta')

eta_i_h <-array(data = 0, dim = c(length(HIV_SET),
                                  length(HIV_SET)))

lapply(1:nrow(eta_params), function(x){
  i <- as.integer(eta_params$HIV_compartment[x]/10)
  h <- eta_params$HIV_compartment[x]%%10
  
  eta_i_h[i,h] <<- eta_params$Reference_expected_value[x]
  #eta_i_h[i,h] <<- .05
})

rm(eta_params)

######Parameters for death and aging rates#######

#Rate of entry into the population due to aging into TB compartment t, HIV com- partment h and gender compartment g, per year, ∀t ∈ {1,3,4,5} ⊂ TB,h ∈ {1,2} ⊂ HIV,g ∈ G
rho_t_h_g <-array(data = 0, dim = c(length(TB_SET),
                                    length(HIV_SET),
                                    length(G_SET)))

#filter dataframe for rho params
rho_params <- param_df%>%
  filter(notation == 'rho')

lapply(c(TB_SUBSET_UNINFECTED_NOIPT, TB_SUBSET_LTBI), function(t){
  lapply(c(HIV_SUBSET_NEG, HIV_SUBSET_POS_NOART), function(h){
    lapply(G_SET, function(g){
      temp <- rho_params%>%
        filter(TB_compartment == t,
               HIV_compartment == h,
               G_compartment == g)
      #rho_t_h_g[t,h,g] <<-temp$Reference_expected_value
      rho_t_h_g[t,h,g] <<- .05
    })
  })
})

rm(rho_params)

#Mortality rates from populations in TB compartment t and HIV compartment h and gender compartment g, per year, ∀t ∈ T B, ∀h ∈ HIV, g ∈ G
mu_t_h_g<-array(data = 0, dim = c(length(TB_SET),
                                  length(HIV_SET),
                                  length(G_SET)))

#filter dataframe for mu params
mu_params <- param_df%>%
  filter(notation == 'mu')

lapply(TB_SET, function(t){
  lapply(HIV_SET, function(h){
    lapply(G_SET, function(g){
      
      temp <- mu_params%>%
        filter(TB_compartment == t,
               HIV_compartment == h,
               G_compartment == g)
      
      mu_t_h_g[t,h,g] <<-temp$Reference_expected_value
    })
  })
})

rm(mu_params)

#Rate of exit from the population due to aging
alpha <- 1/(65-15)


########Force of Infection (FOI) CALCULATIONS#####
FOI_DS <- function(active_pop, total_pop){
  n = 0
  lambda_1_1_numerator<-0
  lambda_1_2_numerator<-0
  
  lapply(HIV_SET, function(h){
    lapply(G_SET, function(g){
      n = n+1
      if (g == 1){
        lambda_1_1_numerator <<- phi_h[h]*active_pop[n] + lambda_1_1_numerator
      }
      else{
        lambda_1_2_numerator <<- phi_h[h]*active_pop[n] + lambda_1_2_numerator
      }
    })
  })
  
  lambda_1_1 <- beta_g[1]*(lambda_1_1_numerator/total_pop)
  lambda_1_2 <- beta_g[2]*(lambda_1_2_numerator/total_pop)
  
  return(c(unname(lambda_1_1), unname(lambda_1_2)))
}
  
FOI_MDR <- function(lambda_1_g){
  lambda_2_g <- c(0,0)
  
  lapply(G_SET, function(g){
    lambda_2_g[g] <<- (varepsilon_g[g]*lambda_1_g[g])/(1-varepsilon_g[g])
  })
  
  return(lambda_2_g)
}


######TB/SEIR model equations#####

#set current policy getting tested 
#we will eventually just run over all policies so we can evaluate against eachother
#but for now we will do one at a time
policy_id <- 1 

tau_itr <- 0 #iterator for recording current time

seir <- function(time, compartment_pop, parameters) {
  with(as.list(c(compartment_pop, parameters)), {
    
    # force of infection 
    active_pop <- compartment_pop[(5*2*4*2+1):((5*2*4*2)+(4*2))] # total active TB compartment populations
    total_pop <- sum(compartment_pop) #sum total population
    
    lambda_1_g <- FOI_DS(active_pop, total_pop)
    lambda_2_g <- FOI_MDR(lambda_1_g)
    
    #record FOI for evaluation
    tau_itr <- tau_itr + 1 #id loc in TT set getting evaluated for proper recording
    lambda_r_g_tau_p[DR_SUBSET_DS, G_SET, tau_itr, policy_id] <<- lambda_1_g
    lambda_r_g_tau_p[DR_SUBSET_MDR, G_SET, tau_itr, policy_id] <<- lambda_2_g
    
    lambda_r <-c(sum(lambda_1_g), sum(lambda_2_g))
    lambda <- sum(lambda_r)
    
    delta_pop <- rep(0, n_compartments) #net gains/losses in pop in compartment
    names(delta_pop)<-d_comparment_id
    
    
    ### TB compartment 1##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[1, r, h, g]
          
          #since DS is not applicable for unifected compartments assign all DR == 1
          if (r == 1){
            
            #look up relative compartment pop locations in 1D array
            other_hiv_ids <- HIV_SET[-h]
            other_hiv_loc1 <- 
            other_hiv_loc2 <- N_t_r_h_g_ref[1, r, other_hiv_ids[2], g]
            other_hiv_loc3 <- N_t_r_h_g_ref[1, r, other_hiv_ids[3], g]
            
            tb_compartment_2_loc <- N_t_r_h_g_ref[2, r, h, g]
            
            #calculate net gains/losses in compartment
            delta_pop[current_compartment_loc] <<- rho_t_h_g[1,h,g]*total_pop + 
              eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
              eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
              eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
              omega_t_p[1,policy_id]*compartment_pop[tb_compartment_2_loc] -
              #losses
              ((alpha + 
                 mu_t_h_g[1,h,g] + 
                 lambda + 
                 kappa_t_h_g_p[1,h,g,policy_id] + 
                 sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
              
          }
          
          
        })
      })
    })
        
    ### TB compartment 2##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[2, r, h, g]
          
          #since DS is not applicable for unifected compartments
          if (r == 1){
            
            #look up relative compartment pop locations in 1D array 
            other_hiv_ids <- HIV_SET[-h]
            other_hiv_loc1 <- N_t_r_h_g_ref[2, r, other_hiv_ids[1], g]
            other_hiv_loc2 <- N_t_r_h_g_ref[2, r, other_hiv_ids[2], g]
            other_hiv_loc3 <- N_t_r_h_g_ref[2, r, other_hiv_ids[3], g]
            
            tb_compartment_1_loc <- N_t_r_h_g_ref[1, r, h, g]
            
            #calculate net gains/losses in compartment
            delta_pop[current_compartment_loc] <<- eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
              eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
              eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
              kappa_t_h_g_p[1, h, g, policy_id]*compartment_pop[tb_compartment_1_loc] -
              #losses
              ((alpha + 
                  mu_t_h_g[2,h,g] +
                  omega_t_p[2, policy_id] + 
                  iota_r[DR_SUBSET_DS]*lambda_r[DR_SUBSET_DS] +
                  iota_r[DR_SUBSET_MDR]*lambda_r[DR_SUBSET_MDR] +
                  sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
          }
          
          
        })
      })
    })

    
    ### TB compartment 3##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[3, r, h, g]
          
          #look up relative compartment pop locations in 1D array  
          other_hiv_ids <- HIV_SET[-h]
          other_hiv_loc1 <- N_t_r_h_g_ref[3, r, other_hiv_ids[1], g]
          other_hiv_loc2 <- N_t_r_h_g_ref[3, r, other_hiv_ids[2], g]
          other_hiv_loc3 <- N_t_r_h_g_ref[3, r, other_hiv_ids[3], g]
          
          tb_compartment_1_loc<-N_t_r_h_g_ref[1,1,h,g] #we ref DR compartment 1 since we assume unifected is in that compartment
          tb_compartment_2_loc<-N_t_r_h_g_ref[2,1,h,g]
          tb_compartment_4_loc<-N_t_r_h_g_ref[4,r,h,g]
          tb_compartment_7_loc<-N_t_r_h_g_ref[7,r,h,g]
          tb_compartment_8_loc<-N_t_r_h_g_ref[8,r,h,g]
          
          #calculate net gains/losses in compartment
          delta_pop[current_compartment_loc] <<- rho_t_h_g[3,h,g]*total_pop + 
            eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
            eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
            eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
            lambda_r[r]*compartment_pop[tb_compartment_1_loc] +
            iota_r[r]*lambda_r[r]*compartment_pop[tb_compartment_2_loc]+
            ((zeta*lambda_r[r])*sum(compartment_pop[c(tb_compartment_4_loc, tb_compartment_7_loc)]))+
            ((upsilon*lambda_r[r])*compartment_pop[tb_compartment_8_loc]) -
            #losses
            ((alpha + 
                mu_t_h_g[3,h,g] + 
                pi_i_t[3,4] + 
                kappa_t_h_g_p[3, h, g, policy_id] +
                theta_h[h]*pi_i_t[3,6] +
                sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
        })
      })
    })
    
    ### TB compartment 4##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[4, r, h, g]
          
          #look up relative compartment pop locations in 1D array 
          other_hiv_ids <- HIV_SET[-h]
          other_hiv_loc1 <- N_t_r_h_g_ref[4, r, other_hiv_ids[1], g]
          other_hiv_loc2 <- N_t_r_h_g_ref[4, r, other_hiv_ids[2], g]
          other_hiv_loc3 <- N_t_r_h_g_ref[4, r, other_hiv_ids[3], g]
          
          tb_compartment_3_pop_loc <- N_t_r_h_g_ref[3,r,h,g]
          
          #calculate net gains/losses in compartment
          delta_pop[current_compartment_loc] <<- rho_t_h_g[4,h,g]*total_pop +
            eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
            eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
            eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
            pi_i_t[3,4]*compartment_pop[tb_compartment_3_pop_loc] -
            #losses
            ((alpha + 
                mu_t_h_g[4,h,g] + 
                zeta*lambda_r[r] +
                kappa_t_h_g_p[4,h,g,policy_id]+
                theta_h[h]*pi_i_t[4,6]+
                sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
        })
      })
    })
    
    ### TB compartment 5##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[5, r, h, g]
          
          #look up relative compartment pop locations in 1D array 
          other_hiv_ids <- HIV_SET[-h]
          other_hiv_loc1 <- N_t_r_h_g_ref[5, r, other_hiv_ids[1], g]
          other_hiv_loc2 <- N_t_r_h_g_ref[5, r, other_hiv_ids[2], g]
          other_hiv_loc3 <- N_t_r_h_g_ref[5, r, other_hiv_ids[3], g]
          
          tb_compartment_3_pop_loc <- N_t_r_h_g_ref[3, r, h, g]
          tb_compartment_4_pop_loc <- N_t_r_h_g_ref[4, r, h, g]
          
          #calculate net gains/losses in compartment
          delta_pop[current_compartment_loc] <<- eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
            eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
            eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
            kappa_t_h_g_p[3,h,g,policy_id]*compartment_pop[tb_compartment_3_pop_loc] +
            kappa_t_h_g_p[4,h,g,policy_id]*compartment_pop[tb_compartment_4_pop_loc] -
            ((alpha + 
                mu_t_h_g[5,h,g] +
                theta_h[h]*pi_i_t[5,6] +
                gamma_r[r]*omega_t_p[5,policy_id] +
                sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
        })
      })
    })
    
    ### TB compartment 6##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[6, r, h, g]
          
          #look up relative compartment pop locations in 1D array 
          other_hiv_ids <- HIV_SET[-h]
          other_hiv_loc1 <- N_t_r_h_g_ref[6, r, other_hiv_ids[1], g]
          other_hiv_loc2 <- N_t_r_h_g_ref[6, r, other_hiv_ids[2], g]
          other_hiv_loc3 <- N_t_r_h_g_ref[6, r, other_hiv_ids[3], g]
          
          tb_compartment_3_pop_loc <- N_t_r_h_g_ref[3, r, h, g]
          tb_compartment_4_pop_loc <- N_t_r_h_g_ref[4, r, h, g]
          tb_compartment_5_pop_loc <- N_t_r_h_g_ref[3, r, h, g]
          
          #calculate net gains/losses in compartment
          delta_pop[current_compartment_loc] <<- rho_t_h_g[6,h,g]*total_pop + 
            eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
            eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
            eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
            theta_h[h]*pi_i_t[3,6]*compartment_pop[tb_compartment_3_pop_loc]+
            theta_h[h]*pi_i_t[4,6]*compartment_pop[tb_compartment_4_pop_loc]+
            theta_h[h]*pi_i_t[5,6]*compartment_pop[tb_compartment_5_pop_loc] -
            #losses
            ((alpha +
                mu_t_h_g[6,h,g] +
                pi_i_t[6,7] +
                sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
        })
      })
    })
    
    ### TB compartment 7##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[7, r, h, g]
          
          #look up relative compartment pop locations in 1D array 
          other_hiv_ids <- HIV_SET[-h]
          other_hiv_loc1 <- N_t_r_h_g_ref[7, r, other_hiv_ids[1], g]
          other_hiv_loc2 <- N_t_r_h_g_ref[7, r, other_hiv_ids[2], g]
          other_hiv_loc3 <- N_t_r_h_g_ref[7, r, other_hiv_ids[3], g]
          
          tb_compartment_6_pop_loc <- N_t_r_h_g_ref[6,r,h,g]
          
          #calculate net gains/losses in compartment
          delta_pop[current_compartment_loc] <<- eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
            eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
            eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
            pi_i_t[6,7]*compartment_pop[tb_compartment_6_pop_loc] -
            #losses
            ((alpha + 
                mu_t_h_g[7,h,g] + 
                zeta*lambda_r[r] +
                sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
        })
      })
    })
    
    ### TB compartment 8##
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #lookup current compartment loc in 1D array
          current_compartment_loc <- N_t_r_h_g_ref[8, r, h, g]
          
          #look up relative compartment pop locations in 1D array  
          other_hiv_ids <- HIV_SET[-h]
          other_hiv_loc1 <- N_t_r_h_g_ref[8, r, other_hiv_ids[1], g]
          other_hiv_loc2 <- N_t_r_h_g_ref[8, r, other_hiv_ids[2], g]
          other_hiv_loc3 <- N_t_r_h_g_ref[8, r, other_hiv_ids[3], g]
          
          tb_compartment_5_pop_loc <- N_t_r_h_g_ref[5,r,h,g]
          
          #calculate net gains/losses in compartment
          delta_pop[current_compartment_loc] <<- eta_i_h[other_hiv_ids[1], h]*compartment_pop[other_hiv_loc1] +
            eta_i_h[other_hiv_ids[2], h]*compartment_pop[other_hiv_loc2] + 
            eta_i_h[other_hiv_ids[3], h]*compartment_pop[other_hiv_loc3] +
            gamma_r[r]*omega_t_p[8, policy_id]*compartment_pop[tb_compartment_5_pop_loc] -
            ((alpha + 
                mu_t_h_g[8,h,g] + 
                upsilon*lambda_r[r] + 
                sum(eta_i_h[h, other_hiv_ids]))*compartment_pop[current_compartment_loc])
        })
      })
    })

    return(list(c(delta_pop)))
  })
}

out <- ode(y = compartments_init, 
           times = TT_SET, 
           func = seir, 
           parms = c())

out <- as.data.frame(out)
