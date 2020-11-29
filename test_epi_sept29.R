##Ref HIV_TB_Model_Appendix_Sept25.pdf##
#Last Modified: Sept 25

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 'readxl', 'stringr', 'reshape2', 'ggplot2', 'varhandle', 'here'), require, character.only=T)

#define input and output directories
indir<-'param_files/'
outdir<-'test_outputs/epi_model'

#set working directory
setwd(here(indir))

#create parameter list and for ode function to be updated bby parameter list gen and pop init gen functions
parameter_list <- list() #parameters passed to ode function
N_t_r_h_g_init <-c() #initial population states
compartment_names <-c()

#define set of compartments
TB_SET <-1:8
DR_SET <-1:2
HIV_SET <-1:4
G_SET <-1:2

#number of compartments
n_compartments= 8*2*4*2

#define set of policies
P_SET <-1:3

#create matrix to look up location of compartment in 1-D array
#8 TB compartments, 2 DR compartments, 4 HIV compartments, 2 gender compartments
N_t_r_h_g_ref <-array(data = 0, dim = c(8,2,4,2))

#generate time inputs into model
TT<-5
time_interval <- 1/12
TT_SET <- c(1:(12*TT))*time_interval

policy_id <-1

######PARAMETER GEN FUNCTION
parameter_list_gen <- function(){
  
  #read in data
  param_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'model_matched_parameters')
  
  #clean dataframe
  #remove all spaces from column names
  names(param_df)<-str_replace_all(names(param_df), c(" " = "_" , "-" = "_" ))
  
  #make sure all compartments are integer type for proper indexing#
  param_df$TB_compartment<-as.integer(param_df$TB_compartment)
  param_df$DR_compartment<-as.integer(param_df$DR_compartment)
  param_df$HIV_compartment<-as.integer(param_df$HIV_compartment)
  param_df$G_compartment<-as.integer(param_df$G_compartment)
  param_df$P_compartment<-as.integer(param_df$P_compartment)

  #############MODEL PARAMETERS#############
  
  #########Parameters that impact force of infection#######
  
  #Number of effective contact for TB transmission per infectious year for gender g
  beta_g <- array(data = 0, dim = length(G_SET))
  
  #filter dataframe for beta params
  beta_params <- param_df%>%
    filter(notation == 'beta')
  
  lapply(1:nrow(beta_params), function(x){
    value <- beta_params$Reference_expected_value[x]
    g <- beta_params$G_compartment[x]
    
    beta_g[g] <<- value
  })
  
  #add to parameter list
  parameter_list[['beta_g']]<<-beta_g
  
  #Relative transmissibility of TB in populations living in HIV compartment h
  phi_h <-array(data = 0, dim = length(HIV_SET))
  
  #filter dataframe for phi params
  phi_params <- param_df%>%
    filter(notation == 'phi')
  
  lapply(1:nrow(phi_params), function(x){
    value <- phi_params$Reference_expected_value[x]
    h <- phi_params$HIV_compartment[x]
    
    phi_h[h]<<-value
  })
  
  #add to parameter list
  parameter_list[['phi_h']]<<-phi_h
  
  #Fraction of new TB infections that are MDR-TB
  varepsilon_g <-array(data = 0, dim = length(G_SET))
  
  #filter dataframe for varepsilon params
  varepsilon_params <- param_df%>%
    filter(notation == 'varepsilon')
  
  lapply(1:nrow(varepsilon_params), function(x){
    value <- varepsilon_params$Reference_expected_value[x]
    g <- varepsilon_params$G_compartment[x]
    
    varepsilon_g[g]<<- value
  })
  
  #add to parameter list
  parameter_list[['varepsilon_g']]<<-varepsilon_g
  
  #Indicator for whether infection with given TB strain can occur while on IPT for populations in DR compartment r
  iota_r <-array(data = 0, dim = length(DR_SET))
  
  #filter dataframe for iota params
  iota_params <- param_df%>%
    filter(notation == 'iota')
  
  lapply(1:nrow(iota_params), function(x){
    value <- iota_params$Reference_expected_value[x]
    r <- iota_params$DR_compartment[x]
    
    iota_r[r] <<- value
  })
  
  #add to parameter list
  parameter_list[['iota_r']]<<-iota_r
  
  #Indicator the diminished force of infection due to partially-protective effects of IPT after moving off of IPT for populations with LTBI
  #filter dataframe for upsilon params
  upsilon_params <- param_df%>%
    filter(notation == 'upsilon')
  
  parameter_list[['upsilon']]<<-upsilon_params$Reference_expected_value
  
  #Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection
  #filter dataframe for zeta params
  zeta_params <- param_df%>%
    filter(notation == 'zeta')
  
  parameter_list[['zeta']]<<-zeta_params$Reference_expected_value
  
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
  
  lapply(1:nrow(kappa_params), function(x){
    value <- kappa_params$Reference_expected_value[x]
    t <- kappa_params$TB_compartment[x]
    h <- kappa_params$HIV_compartment[x]
    g <- kappa_params$G_compartment[x]
    p <- kappa_params$P_compartment[x]
    
    kappa_t_h_g_p[t,h,g,p] <<- value
  })
  
  #add to parameter list
  parameter_list[['kappa_t_h_g_p']]<<-kappa_t_h_g_p
  
  #Rate rate of moving off of IPT from TB compartment t under policy p, per year
  #0 where not applicable
  
  #filter dataframe for omega params
  omega_params <- param_df%>%
    filter(notation == 'omega')
  
  parameter_list[['omega']] <<- omega_params$Reference_expected_value
  
  #Base rates of TB progression of infected populations from TB compartment i to TB compartment t, per year (i, t) ∈ T B (set to zero where not applicable and not included in equations, as shown in Figure 2 and Section 2 respectively)
  
  pi_i_t <- array(data = 0, dim = c(length(TB_SET),
                            length(TB_SET)
                            ))
  
  #filter dataframe for pi params
  pi_params <- param_df%>%
    filter(notation == 'pi')
  
  lapply(1:nrow(pi_params), function(x){
    i <- as.integer(pi_params$TB_compartment[x]/10)
    t <- pi_params$TB_compartment[x]%%10
    value <<-pi_params$Reference_expected_value[x]
    pi_i_t[i,t]<<-value
  })
  
  #add to parameter list
  parameter_list[['pi_i_t']]<<-pi_i_t
  
  #relative risk for TB progression from LTBI to Active for HIV compartment h
  theta_h <- array(data = 0, dim = length(HIV_SET))
  
  #filter dataframe for theta params
  theta_params <- param_df%>%
    filter(notation == 'theta')
  
  lapply(1:nrow(theta_params), function(x){
    value <- theta_params$Reference_expected_value[x]
    h <- theta_params$HIV_compartment[x]

    theta_h[h] <<-value
  })
  
  parameter_list[['theta_h']]<<-theta_h
  
  #1 if in drug-susceptible, DR compartment r ∈ 1 ⊂ DR, 0 if in MDR-TB, DR com- partment r ∈ 2 ⊂ DR to indicate that populations with MDR-TB cannot move into LTBI after IPT
  
  parameter_list[['gamma_r']]<<- c(0,1)
  
  #IPT adherence by gender and policy
  varpi_g_p <- array(data = 0, dim = c(length(G_SET), length(P_SET)))
  
  #filter dataframe for varpi params
  varpi_params <- param_df%>%
    filter(notation == 'varpi')
  
  lapply(1:nrow(varpi_params), function(x){
    value <-varpi_params$Reference_expected_value[x]
    g <- varpi_params$G_compartment[x]
    p <- varpi_params$P_compartment[x]
    varpi_g_p[g,p]<<- value
    
  })
  
  parameter_list[['varpi_g_p']]<<-varpi_g_p
  
  ##########Parameters that describe HIV progression########
  #Rate of populations moving from HIV compartment i to HIV compartment h, per year, (i, h) ∈ HIV under policy p ∈ P #(set to zero where not applicable)
  
  eta_i_h_g_p <-array(data = 0, dim = c(length(HIV_SET),
                                  length(HIV_SET),
                                  length(G_SET),
                                  length(P_SET)))
  
  #filter dataframe for eta params
  eta_params <- param_df%>%
    filter(notation == 'eta')
  
  lapply(1:nrow(eta_params), function(x){
    i <- as.integer(eta_params$HIV_compartment[x]/10)
    h <- eta_params$HIV_compartment[x]%%10
    g <- eta_params$G_compartment[x]
    p <- eta_params$P_compartment[x]
    value <-eta_params$Reference_expected_value[x]
    
    eta_i_h_g_p [i,h,g,p] <<- value
  })
  
  parameter_list[['eta_i_h_g_p']]<<-eta_i_h_g_p
  
  ######Parameters for death and aging rates#######
  
  #Rate of entry into the population due to aging into TB compartment t, HIV com- partment h and gender compartment g, per year, ∀t ∈ {1,3,4,5} ⊂ TB,h ∈ {1,2} ⊂ HIV,g ∈ G
  alpha_in_t_r_h_g <-array(data = 0, dim = c(length(TB_SET),
                                           length(DR_SET),
                                           length(HIV_SET),
                                           length(G_SET)))
  
  #filter dataframe for rho params
  alpha_in_params <- param_df%>%
    filter(notation == 'alpha^in')
  
  lapply(1:nrow(alpha_in_params), function(x){
    t <- alpha_in_params$TB_compartment[x]
    r <- alpha_in_params$DR_compartment[x]
    h <- alpha_in_params$HIV_compartment[x]
    g <- alpha_in_params$G_compartment[x]
    value <- as.double(alpha_in_params$Reference_expected_value[x])
    
    alpha_in_t_r_h_g[t,r,h,g] <<- value
    
  })
  
  parameter_list[['alpha_in_t_r_h_g']]<<-alpha_in_t_r_h_g
  
  
  #Mortality rates from populations in TB compartment t and HIV compartment h and gender compartment g, per year, ∀t ∈ T B, ∀h ∈ HIV, g ∈ G
  mu_t_h_g<-array(data = 0, dim = c(length(TB_SET),
                                  length(HIV_SET),
                                  length(G_SET)))
  
  #filter dataframe for mu params
  mu_params <- param_df%>%
    filter(notation == 'mu')
  
  lapply(1:nrow(mu_params), function(x){
    
    t <- mu_params$TB_compartment[x]
    h <- mu_params$HIV_compartment[x]
    g <- mu_params$G_compartment[x]
    value <- mu_params$Reference_expected_value[x]
        
    mu_t_h_g[t,h,g] <<-value
    })
  
  parameter_list[['mu_t_h_g']]<<-mu_t_h_g
  
  #Rate of exit from the population due to aging
  alpha_out_params <- param_df%>%
    filter(Matched_Model_Parameter == 'alpha^out')
  
  parameter_list[['alpha_out']]<<-alpha_out_params$Reference_expected_value

}

parameter_list_gen()

pop_init_gen<-function(){
  
  #read in pop init data frame
  pop_init_df <- read_excel("Epi_model_parameters.xlsx", sheet = 'pop_init')
  
  #data clean pop init
  names(pop_init_df)<-str_replace_all(names(pop_init_df), c(" " = "_" , "-" = "_" ))
  pop_init_df$TB_compartment<-as.integer(pop_init_df$TB_compartment)
  pop_init_df$DR_compartment<-as.integer(pop_init_df$DR_compartment)
  pop_init_df$HIV_compartment<-as.integer(pop_init_df$HIV_compartment)
  pop_init_df$G_compartment<-as.integer(pop_init_df$G_compartment)
  
  compartment_names_temp <- c()
  N_t_r_h_g_init_temp <-c()
  n <-1
  
  lapply(1:nrow(pop_init_df), function(x){
    t <-pop_init_df$TB_compartment[x]
    r <-pop_init_df$DR_compartment[x]
    h <-pop_init_df$HIV_compartment[x]
    g <-pop_init_df$G_compartment[x]
    value <- pop_init_df$initialized_population_in_compartment[x]
    
    N_t_r_h_g_init_temp <<- c(N_t_r_h_g_init_temp, value)
    compartment_names_temp <<- c(compartment_names_temp, paste0("N_", t, "_", r, "_", h, "_", g))
    N_t_r_h_g_ref[t,r,h,g] <<- n
    n <<- n + 1
    
  })
  
  names(N_t_r_h_g_init_temp)<-compartment_names_temp
  N_t_r_h_g_init <<- N_t_r_h_g_init_temp
  compartment_names <<- compartment_names_temp
  
}

pop_init_gen()


#################2 Model Equaions############
########2.1 aging in calcs##########
aging_in_calcs<-function(compartment_pop){
  N_in_total_temp <- parameter_list$alpha_out*sum(compartment_pop)
  lapply(TB_SET, function(t){
    lapply(HIV_SET, function(h){
      lapply(G_SET, function(g){
        N_in_total_temp <<- N_in_total_temp + (parameter_list$mu_t_h_g[t,h,g]*sum(compartment_pop[N_t_r_h_g_ref[t,DR_SET,h,g]]))
    })
  })
})
  return(N_in_total_temp)
}


######2.2 Force of infection calculations#########
FOI_DS <- function(active_pop, total_pop){
  
  n = 1
  
  numerator_MALE <- rep(0, times = length(HIV_SET))
  numerator_FEMALE <- rep(0, times = length(HIV_SET))
  
  
  lapply(G_SET, function(g){
    lapply(HIV_SET, function(h){
      
      if (g == 1){
        numerator_MALE[h] <<- parameter_list$phi_h[h] * active_pop[n]
        
      } else {
        numerator_FEMALE[h] <<- parameter_list$phi_h[h]*active_pop[n]
      }
      
      n <<- n + 1
      
    })
  })
  
  lambda_1_1 <- parameter_list$beta_g[1]*(sum(numerator_MALE)/total_pop)
  lambda_1_2 <- parameter_list$beta_g[2]*(sum(numerator_FEMALE)/total_pop)
  
  return(c(lambda_1_1, lambda_1_2)) 
}
  
FOI_MDR <- function(lambda_1_g){
  lambda_2_g <- c(0,0)
  
  lapply(G_SET, function(g){
    lambda_2_g[g] <<- (parameter_list$varepsilon_g[g]*lambda_1_g[g])/(1-parameter_list$varepsilon_g[g])
  })
  
  return(lambda_2_g)
}

seir_model <- function(time, compartment_pop, parameters) {
  with(as.list(c(compartment_pop, parameters)), {
    
    #calculate time dependent parameters that are impacted by compartment pop
    
    #aging in
    N_in_total <- aging_in_calcs(compartment_pop)
    
    #FOI
    total_pop <- sum(compartment_pop)
    active_pop <- sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,HIV_SET,G_SET]])
    lambda_1_g <- FOI_DS(active_pop, total_pop) #equation 0a
    lambda_2_g <- FOI_MDR(lambda_1_g) #equation 0b
    lambda_r <-c(sum(lambda_1_g), sum(lambda_2_g)) #equation 0c
    lambda <- sum(lambda_r) #equation 0d
    
    #initiate array for recording changes in each compartment
    delta_pop <- rep(0, n_compartments) #net gains/losses in pop in compartment
    names(delta_pop)<-paste0('d', compartment_names)
    
    ################ TB compartment 1 (Uninfected no IPT) #########################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #since DS is not applicable for unifected compartments
          if (r == 1){
          
          #ID other HIV locations (i) for movements between HIV states
          other_hiv_ids <- HIV_SET[-h]
          
          #Equation 1 
          delta_pop[N_t_r_h_g_ref[1,r,h,g]] <<- (parameter_list$alpha_in_t_r_h_g[1,r,h,g]*N_in_total) +
            sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[1, r, other_hiv_ids, g]]) + 
            (parameter_list$omega*compartment_pop[N_t_r_h_g_ref[2, r, h, g]]) -
            ((parameter_list$alpha_out + 
                parameter_list$mu_t_h_g[1,h,g] +
                lambda + 
                parameter_list$kappa_t_h_g_p[1,h,g,policy_id] +
                sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
               compartment_pop[N_t_r_h_g_ref[1,r,h,g]])
          }
          else {
            delta_pop[N_t_r_h_g_ref[1,r,h,g]] <<-0
          }
        })
      })
    })
    
    ################## TB compartment 2 (uninfected no IPT)###################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          if (r == 1){
            
            #ID other HIV locations (i) for movements between HIV states
            other_hiv_ids <- HIV_SET[-h]
            
            delta_pop[N_t_r_h_g_ref[2,r,h,g]] <<- sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*compartment_pop[N_t_r_h_g_ref[2, r, other_hiv_ids, g]]) +
              (parameter_list$kappa_t_h_g_p[1,h,g, policy_id]*compartment_pop[N_t_r_h_g_ref[1,r,h,g]]) -
              ((parameter_list$alpha_out +
                 parameter_list$mu_t_h_g[2,h,g] +
                 sum(parameter_list$iota_r*lambda_r) +
                 sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
              compartment_pop[N_t_r_h_g_ref[2,r,h,g]])
              
          }
          else {
            delta_pop[N_t_r_h_g_ref[2,r,h,g]] <<-0
          }
        })
      })
    })
    
    ################# TB compartment 3 (infected recently) #######################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #ID other HIV locations (i) for movements between HIV states
          other_hiv_ids <- HIV_SET[-h]
          
          delta_pop[N_t_r_h_g_ref[3,r,h,g]] <<- parameter_list$alpha_in_t_r_h_g[3,r,h,g]*N_in_total +
            sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[3, r, other_hiv_ids, g]]) + 
            (lambda_r[r]*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET, h, g]])) +
            (parameter_list$zeta*lambda_r[r]*sum(compartment_pop[N_t_r_h_g_ref[c(4,7),DR_SET, h, g]])) +
            (parameter_list$upsilon*lambda_r[r]*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET, h, g]])) -
            ((parameter_list$alpha_out+
               parameter_list$mu_t_h_g[3,h,g]+
               parameter_list$pi_i_t[3,4] +
               parameter_list$kappa_t_h_g_p[3,h,g,policy_id] +
               (parameter_list$phi_h[h]*parameter_list$varpi_g_p[g,policy_id]*parameter_list$pi_i_t[3,6])+
               sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
               compartment_pop[N_t_r_h_g_ref[3,r,h,g]])
        })
      })
    })
    
    ####################### TB compartment 4 (infected remotely) #################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #ID other HIV locations (i) for movements between HIV states
          other_hiv_ids <- HIV_SET[-h]
          
          delta_pop[N_t_r_h_g_ref[4,r,h,g]] <<- parameter_list$alpha_in_t_r_h_g[4,r,h,g]*N_in_total +
            sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[4, r, other_hiv_ids, g]]) + 
            (parameter_list$pi_i_t[3,4]*compartment_pop[N_t_r_h_g_ref[3,r,h,g]]) -
            ((parameter_list$alpha_out)+
               parameter_list$mu_t_h_g[4,h,g]+
               parameter_list$zeta*lambda +
               parameter_list$kappa_t_h_g_p[4,h,g, policy_id]+
               (parameter_list$phi_h[h]*parameter_list$varpi_g_p[g,policy_id]*parameter_list$pi_i_t[4,6])+
               sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id])*
               compartment_pop[N_t_r_h_g_ref[4,r,h,g]])
        })
      })
    })
    
    ############### TB compartment 5 (LTBI, on IPT) ###################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #ID other HIV locations (i) for movements between HIV states
          other_hiv_ids <- HIV_SET[-h]
          
          delta_pop[N_t_r_h_g_ref[5,r,h,g]] <<- sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[5, r, other_hiv_ids, g]]) + 
            (parameter_list$kappa_t_h_g_p[3,h,g, policy_id]*compartment_pop[N_t_r_h_g_ref[3,r,h,g]]) +
            (parameter_list$kappa_t_h_g_p[4,h,g, policy_id]*compartment_pop[N_t_r_h_g_ref[4,r,h,g]]) -
            ((parameter_list$alpha_out)+
               (parameter_list$mu_t_h_g[5,h,g])+
               (parameter_list$phi_h[h]*parameter_list$varpi_g_p[g,policy_id]*parameter_list$pi_i_t[5,6])+
               (parameter_list$gamma_r[r]*parameter_list$omega)+
               (sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
               compartment_pop[N_t_r_h_g_ref[5,r,h,g]])
        })
      })
    })
    
    ############## TB compartment 6 (Active TB) #########################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #ID other HIV locations (i) for movements between HIV states
          other_hiv_ids <- HIV_SET[-h]
          
          delta_pop[N_t_r_h_g_ref[6,r,h,g]] <<- parameter_list$alpha_in_t_r_h_g[6,r,h,g]*N_in_total +
            sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[6, r, other_hiv_ids, g]]) + 
            (parameter_list$phi_h[h]*parameter_list$varpi_g_p[g,policy_id]*
               parameter_list$pi_i_t[3,6]*compartment_pop[N_t_r_h_g_ref[3,r,h,g]])+ 
            (parameter_list$phi_h[h]*parameter_list$varpi_g_p[g,policy_id]*
               parameter_list$pi_i_t[4,6]*compartment_pop[N_t_r_h_g_ref[4,r,h,g]])+
            (parameter_list$phi_h[h]*parameter_list$varpi_g_p[g,policy_id]*
               parameter_list$pi_i_t[5,6]*compartment_pop[N_t_r_h_g_ref[5,r,h,g]]) -
            (parameter_list$alpha_out +
               parameter_list$mu_t_h_g[6,h,g] +
               parameter_list$pi_i_t[6,7]+
               (sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
               compartment_pop[N_t_r_h_g_ref[6,r,h,g]])
          
        })
      })
    })
    
    ############### TB compartment 7####################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #look up relative compartment pop locations in 1D array 
          other_hiv_ids <- HIV_SET[-h]
          
          delta_pop[N_t_r_h_g_ref[7,r,h,g]] <<- sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[7, r, other_hiv_ids, g]]) +
            (parameter_list$pi_i_t[6,7]*compartment_pop[N_t_r_h_g_ref[6,r,h,g]]) -
            ((parameter_list$alpha_out) +
               (parameter_list$mu_t_h_g[7,h,g])+
               (parameter_list$zeta*lambda)+
               (sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
               compartment_pop[N_t_r_h_g_ref[7,r,h,g]])
          
        })
      })
    })
    
    ################## TB compartment 8#################
    lapply(DR_SET, function(r){
      lapply(HIV_SET, function(h){
        lapply(G_SET, function(g){
          
          #look up relative compartment pop locations in 1D array  
          other_hiv_ids <- HIV_SET[-h]
          
          delta_pop[N_t_r_h_g_ref[8,r,h,g]] <<- sum(parameter_list$eta_i_h_g_p[other_hiv_ids, h, g, policy_id]*
                  compartment_pop[N_t_r_h_g_ref[8, r, other_hiv_ids, g]]) +
            (parameter_list$gamma_r[r]*parameter_list$omega*compartment_pop[N_t_r_h_g_ref[5,r,h,g]]) -
            ((parameter_list$alpha_out) +
                (parameter_list$mu_t_h_g[8,h,g]) +
               (parameter_list$upsilon*lambda)+
               (sum(parameter_list$eta_i_h_g_p[h, other_hiv_ids, g, policy_id]))*
               compartment_pop[N_t_r_h_g_ref[8,r,h,g]])
          
        })
      })
    })
    
    return(list(c(delta_pop)))
  })
}


out_all_df <- data.frame(matrix(ncol = n_compartments+2, nrow = 0))
colnames(out_all_df) <- c('time', compartment_id, 'policy_id')

lapply(P_SET, function(p){
  #lapply(seq(.001, .009, .001), function(b1){
    #lapply(seq(.001, .009, .001), function(b2){
      
      #beta_g<<-c(b1, b2)
  
      #set policy id
      policy_id <<- p
      
      #reset tau itr
      tau_itr <<- 1
      
      out <- ode(y = N_t_r_h_g_init, 
               times = TT_SET, 
               func = seir_model, 
               parms = NULL,
               method = "lsodar")
      
      out <- as.data.frame(out)
      
      out$policy_id <- rep(policy_id, length(TT_SET))
      #out$beta_1<- rep(b1, length(TT_SET))
      #out$beta_2<- rep(b2, length(TT_SET))
      
      out_all_df <<- rbind(out_all_df, out)
  
    })
  #})
#})

#calculate N(t)/total population at each timestep
out_all_df <- out_all_df%>% mutate(total_pop = rowSums(.[1:129]))
out_all_df$policy_id <-as.factor(out_all_df$policy_id)

