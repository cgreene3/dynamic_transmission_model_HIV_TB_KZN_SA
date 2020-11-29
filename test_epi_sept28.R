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
parameter_list <- c() #parameters passed to ode function
N_t_r_h_g_init <-c() #initial population states
compartment_names <-c()
d_compartment_names<-paste0('d', compartment_names)

#TB SET
TB_SET <-1:8

#DR SET
DR_SET <-1:2

#HIV SET
HIV_SET <-1:4

#G SET
G_SET <-1:2

n_compartments= 8*2*4*2

#create matrix to look up location of compartment in 1-D array
#8 TB compartments, 2 DR compartments, 4 HIV compartments, 2 gender compartments
N_t_r_h_g_ref <-array(data = 0, dim = c(8,2,4,2))

#generate time inputs into model
TT<-5
time_interval <- 1/12
TT_SET <- c(1:(12*TT))*time_interval

######PARAMETER GEN FUNCTION
parameter_list_gen <- function(policy_id){
  
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
  beta_g <- c()
  beta_g_names <- c()
  
  #filter dataframe for beta params
  beta_params <- param_df%>%
    filter(notation == 'beta')
  
  lapply(1:nrow(beta_params), function(x){
    value <- beta_params$Reference_expected_value[x]
    g <- beta_params$G_compartment[x]
    
    beta_g <<- c(beta_g, value)
    beta_g_names <<- c(beta_g_names, paste0('beta_', g))
  })
  
  #make named list
  names(beta_g) <- beta_g_names
  
  #add to parameter list
  parameter_list<<-c(parameter_list, beta_g)
  
  #Relative transmissibility of TB in populations living in HIV compartment h
  phi_h <-c()
  phi_h_names<-c()
  
  #filter dataframe for phi params
  phi_params <- param_df%>%
    filter(notation == 'phi')
  
  lapply(1:nrow(phi_params), function(x){
    value <- phi_params$Reference_expected_value[x]
    h <- phi_params$HIV_compartment[x]
    
    phi_h <<- c(phi_h, value)
    phi_h_names <<- c(phi_h_names, paste0('phi_', h))
  })
  
  names(phi_h) <- phi_h_names
  
  #add to parameter list
  parameter_list<<-c(parameter_list, phi_h)
  
  #Fraction of new TB infections that are MDR-TB
  varepsilon_g <-c()
  varepsilon_g_names <-c()
  
  #filter dataframe for varepsilon params
  varepsilon_params <- param_df%>%
    filter(notation == 'varepsilon')
  
  lapply(1:nrow(varepsilon_params), function(x){
    value <- varepsilon_params$Reference_expected_value[x]
    g <- varepsilon_params$G_compartment[x]
    
    varepsilon_g <<- c(varepsilon_g, value)
    varepsilon_g_names <<- c(varepsilon_g_names, paste0('varepsilon_', g))
  })
  
  
  names(varepsilon_g) <- varepsilon_g_names
  
  #add to parameter list
  parameter_list<<-c(parameter_list, varepsilon_g)
  
  #Indicator for whether infection with given TB strain can occur while on IPT for populations in DR compartment r
  iota_r <-c()
  iota_r_names <-c()
  
  #filter dataframe for iota params
  iota_params <- param_df%>%
    filter(notation == 'iota')
  
  lapply(1:nrow(iota_params), function(x){
    value <- iota_params$Reference_expected_value[x]
    r <- iota_params$DR_compartment[x]
    
    iota_r <<- c(iota_r, value)
    iota_r_names <<- c(iota_r_names, paste0('iota_', r))
  })
  
  names(iota_r) <- iota_r_names
  
  #add to parameter list
  parameter_list<<-c(parameter_list, iota_r)
  
  #Indicator the diminished force of infection due to partially-protective effects of IPT after moving off of IPT for populations with LTBI
  #filter dataframe for upsilon params
  upsilon_params <- param_df%>%
    filter(notation == 'upsilon')
  
  upsilon <- upsilon_params$Reference_expected_value
  names(upsilon) <- 'upsilon'
  
  #add to parameter list
  parameter_list<<-c(parameter_list, upsilon)
  
  #Indicator that diminishes force of infection due to the partially-protective effect of LTBI infection and acquiring a new TB infection
  #filter dataframe for zeta params
  zeta_params <- param_df%>%
    filter(notation == 'zeta')
  
  zeta <- zeta_params$Reference_expected_value
  names(zeta) <- 'zeta'
  
  #add to parameter list
  parameter_list<<-c(parameter_list, zeta)
  
  ###########Parameters that describe TB progression##############
  #Rate of IPT initiation from TB compartment t and HIV compartment h for gender g under policy p, per year
  #0 where not applicable
  
  kappa_t_h_g_p <- c()
  
  kappa_t_h_g_p_names <- c()
  
  #filter dataframe for kappa params
  kappa_params <- param_df%>%
    filter(notation == 'kappa')%>%
    filter(P_compartment == policy_id)
  
  lapply(1:nrow(kappa_params), function(x){
    value <- kappa_params$Reference_expected_value[x]
    t <- kappa_params$TB_compartment[x]
    h <- kappa_params$HIV_compartment[x]
    g <- kappa_params$G_compartment[x]
    
    kappa_t_h_g_p <<- c(kappa_t_h_g_p, value)
    kappa_t_h_g_p_names <<- c(kappa_t_h_g_p_names, paste0('kappa_',t,"_",h,"_",g))
  })
  
  names(kappa_t_h_g_p)<-kappa_t_h_g_p_names
  
  #add to parameter list
  parameter_list<<-c(parameter_list, kappa_t_h_g_p)
  
  #Rate rate of moving off of IPT from TB compartment t under policy p, per year
  #0 where not applicable
  
  #filter dataframe for omega params
  omega_params <- param_df%>%
    filter(notation == 'omega')
  
  omega <- omega_params$Reference_expected_value
  names(omega) <-'omega'
  
  #add to parameter list
  parameter_list<<-c(parameter_list, omega)
  
  #Base rates of TB progression of infected populations from TB compartment i to TB compartment t, per year (i, t) ∈ T B (set to zero where not applicable and not included in equations, as shown in Figure 2 and Section 2 respectively)
  
  #filter dataframe for pi params
  pi_params <- param_df%>%
    filter(notation == 'pi')
  
  pi_i_t <- c()
  pi_i_t_names <- c()
  
  lapply(1:nrow(pi_params), function(x){
    i <- as.integer(pi_params$TB_compartment[x]/10)
    t <- pi_params$TB_compartment[x]%%10
    pi_i_t <<- c(pi_i_t, pi_params$Reference_expected_value[x])
    pi_i_t_names <<- c(pi_i_t_names, paste0('pi_', i, "_", t))
  })
  
  names(pi_i_t) <- pi_i_t_names
  
  #add to parameter list
  parameter_list<<-c(parameter_list, pi_i_t)
  
  #relative risk for TB progression from LTBI to Active for HIV compartment h
  theta_h <- c()
  theta_h_names <-c()
  
  #filter dataframe for theta params
  theta_params <- param_df%>%
    filter(notation == 'theta')
  
  lapply(1:nrow(theta_params), function(x){
    value <- theta_params$Reference_expected_value[x]
    h <- theta_params$HIV_compartment[x]

    theta_h <<-c(theta_h, value)
    theta_h_names<<-c(theta_h_names, paste0('theta_', h))
  })
  
  names(theta_h) <- theta_h_names
  
  parameter_list<<-c(parameter_list, theta_h)
  
  #1 if in drug-susceptible, DR compartment r ∈ 1 ⊂ DR, 0 if in MDR-TB, DR com- partment r ∈ 2 ⊂ DR to indicate that populations with MDR-TB cannot move into LTBI after IPT
  gamma_r <- c(0,1)
  names(gamma_r)<-c('gamma_1', 'gamma_2')
  
  parameter_list<<-c(parameter_list, gamma_r)
  
  #IPT adherence by gender and policy
  varpi_g_p <- c()
  varpi_g_p_names <-c()
  
  #filter dataframe for varpi params
  varpi_params <- param_df%>%
    filter(notation == 'varpi')%>%
    filter(P_compartment == policy_id)
  
  lapply(1:nrow(varpi_params), function(x){
    value <-varpi_params$Reference_expected_value[x]
    g <- varpi_params$G_compartment[x]
    
    varpi_g_p<<-c(varpi_g_p, value)
    varpi_g_p_names<<-c(varpi_g_p_names, paste0("varpi_", g))
    
  })
  
  names(varpi_g_p) <- varpi_g_p_names
  
  parameter_list<<-c(parameter_list, varpi_g_p)
  
  ##########Parameters that describe HIV progression########
  #Rate of populations moving from HIV compartment i to HIV compartment h, per year, (i, h) ∈ HIV under policy p ∈ P #(set to zero where not applicable)
  
  #filter dataframe for eta params
  eta_params <- param_df%>%
    filter(notation == 'eta')%>%
    filter(P_compartment == policy_id)
  
  eta_i_h_g_p <-c()
  eta_i_h_g_p_names <-c()
  
  lapply(1:nrow(eta_params), function(x){
    i <- as.integer(eta_params$HIV_compartment[x]/10)
    h <- eta_params$HIV_compartment[x]%%10
    g <- eta_params$G_compartment[x]
    
    eta_i_h_g_p <<- c(eta_i_h_g_p, eta_params$Reference_expected_value[x])
    eta_i_h_g_p_names <<- c(eta_i_h_g_p_names, paste0('eta_', i, "_", h, "_", g))
  })
  
  names(eta_i_h_g_p)<-eta_i_h_g_p_names
  
  parameter_list<<-c(parameter_list, eta_i_h_g_p)
  
  ######Parameters for death and aging rates#######
  
  #Rate of entry into the population due to aging into TB compartment t, HIV com- partment h and gender compartment g, per year, ∀t ∈ {1,3,4,5} ⊂ TB,h ∈ {1,2} ⊂ HIV,g ∈ G
  alpha_in_t_r_h_g <-c()
  alpha_in_t_r_h_g_names <-c()
  
  #filter dataframe for rho params
  alpha_in_params <- param_df%>%
    filter(notation == 'alpha^in')
  
  lapply(1:nrow(alpha_in_params), function(x){
    t <- alpha_in_params$TB_compartment[x]
    r <- alpha_in_params$DR_compartment[x]
    h <- alpha_in_params$HIV_compartment[x]
    g <- alpha_in_params$G_compartment[x]
    value <- as.double(alpha_in_params$Reference_expected_value[x])
    
    alpha_in_t_r_h_g <<- c(alpha_in_t_r_h_g, value)
    alpha_in_t_r_h_g_names <<- c(alpha_in_t_r_h_g_names, paste0("alpha_in_", t, "_", r, "_", h, "_", g))
    
  })
  
  names(alpha_in_t_r_h_g) <- alpha_in_t_r_h_g_names
  
  parameter_list<<-c(parameter_list, alpha_in_t_r_h_g)
  
  
  #Mortality rates from populations in TB compartment t and HIV compartment h and gender compartment g, per year, ∀t ∈ T B, ∀h ∈ HIV, g ∈ G
  mu_t_h_g<-c()
  mu_t_h_g_names<-c()
  
  #filter dataframe for mu params
  mu_params <- param_df%>%
    filter(notation == 'mu')
  
  lapply(1:nrow(mu_params), function(x){
    
    t <- mu_params$TB_compartment[x]
    h <- mu_params$HIV_compartment[x]
    g <- mu_params$G_compartment[x]
    value <- mu_params$Reference_expected_value[x]
        
    mu_t_h_g <<-c(mu_t_h_g, value)
    mu_t_h_g_names <<- c(mu_t_h_g_names, paste0("mu_", t, "_", h, "_", g))
    })
  
  names(mu_t_h_g) <- mu_t_h_g_names
  
  parameter_list<<-c(parameter_list, mu_t_h_g)
  
  #Rate of exit from the population due to aging
  alpha_out_params <- param_df%>%
    filter(Matched_Model_Parameter == 'alpha^out')
  
  alpha_out <- alpha_out_params$Reference_expected_value
  names(alpha_out) <- 'alpha_out'
  
  parameter_list<<-c(parameter_list, alpha_out)

}

parameter_list_gen(policy_id = 1)

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


#########SEIR model equations###########
seir_model <- function(time, compartment_pop, parameters) {
  with(as.list(c(compartment_pop, parameters)), {
  
  ##########calculate N_in (equation 0)##############
  N_in = alpha_out*(sum(compartment_pop))+
    #all tb compartments, hiv compartment 1 (hiv-), gender compartment 1 (male)
    mu_1_1_1*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,1,1]])+
    mu_2_1_1*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,1,1]])+
    mu_3_1_1*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,1,1]])+
    mu_4_1_1*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,1,1]])+
    mu_5_1_1*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,1,1]])+
    mu_6_1_1*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,1,1]])+
    mu_7_1_1*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,1,1]])+
    mu_8_1_1*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,1,1]])+
    #all tb compartments, hiv compartment 2 (hiv+, not on ART, CD4 >= 200), gender compartment 1 (male)
    mu_1_2_1*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,2,1]])+
    mu_2_2_1*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,2,1]])+
    mu_3_2_1*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,2,1]])+
    mu_4_2_1*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,2,1]])+
    mu_5_2_1*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,2,1]])+
    mu_6_2_1*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,2,1]])+
    mu_7_2_1*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,2,1]])+
    mu_8_2_1*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,2,1]])+
    #all tb compartments, hiv compartment 3 (hiv+, not on ART, CD4 < 200), gender compartment 1 (male)
    mu_1_3_1*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,3,1]])+
    mu_2_3_1*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,3,1]])+
    mu_3_3_1*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,3,1]])+
    mu_4_3_1*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,3,1]])+
    mu_5_3_1*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,3,1]])+
    mu_6_3_1*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,3,1]])+
    mu_7_3_1*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,3,1]])+
    mu_8_3_1*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,3,1]])+
    #all tb compartments, hiv compartment 4 (hiv+, on ART), gender compartment 1 (male)
    mu_1_4_1*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,4,1]])+
    mu_2_4_1*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,4,1]])+
    mu_3_4_1*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,4,1]])+
    mu_4_4_1*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,4,1]])+
    mu_5_4_1*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,4,1]])+
    mu_6_4_1*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,4,1]])+
    mu_7_4_1*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,4,1]])+
    mu_8_4_1*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,4,1]])+
    #all tb compartments, hiv compartment 1 (hiv-), gender compartment 2 (female)
    mu_1_1_2*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,1,2]])+
    mu_2_1_2*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,1,2]])+
    mu_3_1_2*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,1,2]])+
    mu_4_1_2*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,1,2]])+
    mu_5_1_2*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,1,2]])+
    mu_6_1_2*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,1,2]])+
    mu_7_1_2*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,1,2]])+
    mu_8_1_2*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,1,2]])+
    #all tb compartments, hiv compartment 2 (hiv+, not on ART, CD4 >= 200), gender compartment 2 (female)
    mu_1_2_2*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,2,1]])+
    mu_2_2_2*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,2,1]])+
    mu_3_2_2*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,2,1]])+
    mu_4_2_2*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,2,1]])+
    mu_5_2_2*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,2,1]])+
    mu_6_2_2*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,2,1]])+
    mu_7_2_2*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,2,1]])+
    mu_8_2_2*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,2,1]])+
    #all tb compartments, hiv compartment 3 (hiv+, not on ART, CD4 < 200), gender compartment 2 (female)
    mu_1_3_2*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,3,2]])+
    mu_2_3_2*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,3,2]])+
    mu_3_3_2*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,3,2]])+
    mu_4_3_2*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,3,2]])+
    mu_5_3_2*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,3,2]])+
    mu_6_3_2*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,3,2]])+
    mu_7_3_2*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,3,2]])+
    mu_8_3_2*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,3,2]])+
    #all tb compartments, hiv compartment 4 (hiv+, on ART), gender compartment 2 (female)
    mu_1_4_2*sum(compartment_pop[N_t_r_h_g_ref[1,DR_SET,4,2]])+
    mu_2_4_2*sum(compartment_pop[N_t_r_h_g_ref[2,DR_SET,4,2]])+
    mu_3_4_2*sum(compartment_pop[N_t_r_h_g_ref[3,DR_SET,4,2]])+
    mu_4_4_2*sum(compartment_pop[N_t_r_h_g_ref[4,DR_SET,4,2]])+
    mu_5_4_2*sum(compartment_pop[N_t_r_h_g_ref[5,DR_SET,4,2]])+
    mu_6_4_2*sum(compartment_pop[N_t_r_h_g_ref[6,DR_SET,4,2]])+
    mu_7_4_2*sum(compartment_pop[N_t_r_h_g_ref[7,DR_SET,4,2]])+
    mu_8_4_2*sum(compartment_pop[N_t_r_h_g_ref[8,DR_SET,4,2]])
  
  #############calculate FOI##################
  
  #Equation 0a- calculates the force of infection for DS populations for males and females
  #Equation 0a - Males
  lambda_1_1 <- beta_1*((parameter_list[3:6] * 
                           compartment_pop[N_t_r_h_g_ref[6,1,HIV_SET,1]])
                        /sum(compartment_pop)) #note parameter_list[3:6] is phi values
  
  #Equation 0a - Females
  lambda_1_2 <- beta_2*((parameter_list[3:6] * 
                           compartment_pop[N_t_r_h_g_ref[6,1,HIV_SET,2]])
                        /sum(compartment_pop)) #note parameter_list[3:6] is phi values
  
  #Equation 0b - calculates the force of infection for MDR-TB populations
  #Equation 0b - males
  lambda_2_1 <- (varepsilon_1*lambda_1_1)/(1-varepsilon_1)
  
  #equation 0b - females
  lambda_2_2 <- (varepsilon_2*lambda_1_2)/(1-varepsilon_2)
  
  #Equation (0c) calculates the force of infection by DR compartment.
  #Equation 0c
  lambda_1 <- sum(lambda_1_1+lambda_1_2)
  lambda_2 <- sum(lambda_2_1+lambda_2_2)
  
  #Equation (0d) calculates the total force of infection
  lambda <- sum(lambda_1, lambda_2)
  
  ##############TB COMPARTMENT EQUATIONS################
  
  #create named array for derivative functions
  delta_pop <- rep(0, n_compartments) #net gains/losses in pop in compartment
  names(delta_pop)<-d_compartment_names
  
  #######Uninfected, not on IPT- TB state 1###############
  
  #TB = 1, DR = 1, HIV = 1, G =1
  delta_pop[1] <- (alpha_in_1_1_1_1*N_in)+
    (omega*N_2_1_1_1)-
    ((alpha_out + mu_1_1_1 + lambda + kappa_1_1_1 + eta_1_2_1)*N_1_1_1_1)
  
  #TB = 1, DR = 2, HIV = 1, G = 1
  delta_pop[2] <- 0
    
  #TB = 1, DR = 1 , HIV = 2, G = 1
  delta_pop[3] <- (alpha_in_1_1_2_1*N_in)+
    (eta_1_2_1 * N_1_1_1_1)+
    (omega*N_2_1_2_1)-
    ((alpha_out + mu_1_2_1 + lambda + kappa_1_2_1 + eta_2_3_1 + eta_2_4_1)*N_1_1_2_1)
  
  delta_pop[4]<-0
  
  #TB = 1, DR = 1, HIV = 3, G = 1
  delta_pop[5] <- (alpha_in_1_1_3_1*N_in)+
    (eta_2_3_1 * N_1_1_2_1)+
    (omega*N_2_1_3_1)-
    ((alpha_out + mu_1_3_1 + lambda + kappa_1_3_1 +eta_3_4_1)*N_1_1_3_1)
  
  delta_pop[6]<-0
  
  #TB = 1, DR = 1, HIV = 4, G = 1
  delta_pop[7] <- (alpha_in_1_1_4_1*N_in)+
    (eta_2_4_1 * N_1_1_2_1)+
    (eta_3_4_1 * N_1_1_3_1)+
    (omega*N_2_1_4_1)-
    ((alpha_out + mu_1_4_1 + lambda + kappa_1_4_1)*N_1_1_4_1)
  
  delta_pop[8]<-0
  
  #TB = 1, DR = 1, HIV = 1, G =2
  delta_pop[9] <- (alpha_in_1_1_1_2*N_in)+
    (omega*N_2_1_1_2)-
    ((alpha_out + mu_1_1_2 + lambda + kappa_1_1_2 + eta_1_2_2)*N_1_1_1_1)
  
  #TB = 1, DR = 2, HIV = 1, G = 2
  delta_pop[10] <- 0
    
  #TB = 1, DR = 1 , HIV = 2, G = 2
  delta_pop[11] <- (alpha_in_1_1_2_2*N_in)+
    (eta_1_2_2 * N_1_1_1_2)+
    (omega*N_2_1_2_2)-
    ((alpha_out + mu_1_2_2 + lambda + kappa_1_2_2 + eta_2_3_2 + eta_2_4_2)*N_1_1_2_2)
  
  delta_pop[12]<-0
  
  #TB = 1, DR = 1, HIV = 3, G = 2
  delta_pop[13] <- (alpha_in_1_1_3_2*N_in)+
    (eta_2_3_2 * N_1_1_2_2)+
    (omega*N_2_1_3_2)-
    ((alpha_out + mu_1_3_2 + lambda + kappa_1_3_2 +eta_3_4_2)*N_1_1_3_2)
  
  delta_pop[14]<-0
  
  #TB = 1, DR = 1, HIV = 4, G = 2
  delta_pop[15] <- (alpha_in_1_1_4_2*N_in)+
    (eta_2_4_2 * N_1_1_2_2)+
    (eta_3_4_2 * N_1_1_3_2)+
    (omega*N_2_1_4_2)-
    ((alpha_out + mu_1_4_2 + lambda + kappa_1_4_2)*N_1_1_4_2)
  
  delta_pop[16]<-0
  
  #######Uninfected, on IPT - TB state 2###############
  
  #TB = 2, DR = 1, HIV = 1, G =1
  delta_pop[17] <- (kappa_1_1_1*N_1_1_1_1) -
    ((alpha_out + mu_2_1_1 + omega + (iota_1*lambda_1) + (iota_2*lambda_2)+ eta_1_2_1)*N_2_1_1_1)
  
  #TB = 2, DR = 2, HIV = 1, G = 1
  delta_pop[18] <- 0
    
  #TB = 2, DR = 1 , HIV = 2, G = 1
  delta_pop[19] <- (eta_1_2_1 * N_2_1_1_1)+
    (kappa_1_2_1*N_1_1_2_1) -
    ((alpha_out + mu_2_2_1 + omega + (iota_1*lambda_1) + (iota_2*lambda_2)+ eta_2_3_1 + eta_2_4_1)*N_2_1_2_1)
  
  delta_pop[20]<-0
  
  #TB = 2, DR = 1, HIV = 3, G = 1
  delta_pop[21] <- (eta_2_3_1 * N_2_1_2_1)+
    (kappa_1_3_1*N_1_1_3_1) -
    ((alpha_out + mu_2_3_1 + + omega + (iota_1*lambda_1) + (iota_2*lambda_2) + eta_3_4_1)*N_2_1_3_1)
  
  delta_pop[22]<-0
  
  #TB = 2, DR = 1, HIV = 4, G = 1
  delta_pop[23] <- (eta_2_4_1 * N_2_1_2_1)+
    (eta_3_4_1 * N_2_1_3_1)+
    (kappa_1_4_1*N_1_1_4_1) -
    ((alpha_out + mu_2_4_1 + + omega + (iota_1*lambda_1) + (iota_2*lambda_2))*N_2_1_4_1)
  
  delta_pop[24]<-0
  
  #TB = 2, DR = 1, HIV = 1, G =2
  delta_pop[25] <- (kappa_1_1_2*N_1_1_1_2) -
    ((alpha_out + mu_2_1_2 + omega + (iota_1*lambda_1) + (iota_2*lambda_2)+ eta_1_2_2)*N_2_1_1_2)
  
  #TB = 2, DR = 2, HIV = 1, G = 2
  delta_pop[26] <- 0
    
  #TB = 2, DR = 1 , HIV = 2, G = 2
  delta_pop[27] <- (eta_1_2_2 * N_2_1_1_2)+
    (kappa_1_2_2*N_1_1_2_2) -
    ((alpha_out + mu_2_2_2 + omega + (iota_1*lambda_1) + (iota_2*lambda_2)+ eta_2_3_2 + eta_2_4_2)*N_2_1_2_2)
  
  delta_pop[28]<-0
  
  #TB = 2, DR = 1, HIV = 3, G = 2
  delta_pop[29] <- (eta_2_3_2 * N_2_1_2_2)+
    (kappa_1_3_2*N_1_1_3_2) -
    ((alpha_out + mu_2_3_2 + omega + (iota_1*lambda_1) + (iota_2*lambda_2) + eta_3_4_2)*N_2_1_3_2)
  
  delta_pop[30]<-0
  
  #TB = 2, DR = 1, HIV = 4, G = 2
  delta_pop[31] <- (eta_2_4_2 * N_2_1_2_2)+
    (eta_3_4_2 * N_2_1_3_2)+
    (kappa_1_4_2*N_1_1_4_2) -
    ((alpha_out + mu_2_4_2 + omega + (iota_1*lambda_1) + (iota_2*lambda_2))*N_2_1_4_2)
  
  delta_pop[32]<-0
  
  #######LTBI, infected recently - TB state 3###############
  
  #TB = 3, DR = 1, HIV = 1, G = 1
  delta_pop[33]<-(alpha_in_3_1_1_1*N_in)+
    (lambda_1*N_1_1_1_1)+
    (iota_1*lambda_1*N_2_1_1_1)+
    ((zeta*lambda_1)*(N_4_1_1_1 + N_7_1_1_1))+
    (upsilon*lambda_1*N_8_1_1_1)-
    ((alpha_out + mu_3_1_1 + pi_3_4 + kappa_3_1_1 + (theta_1*varpi_g*pi_3_6) + eta_1_2_1)*N_3_1_1_1)
  
  #TB = 3, DR = 2, HIV = 1, G = 1
  delta_pop[33]<-(alpha_in_3_2_1_1*N_in)+
    (lambda_2*N_1_1_1_1)+
    (iota_1*lambda_2*N_2_1_1_1)+
    ((zeta*lambda_2)*(compartment_pop[N_t_r_h_g_ref[c(4,7),DR_SET,1,1]]))+
    (upsilon*lambda_1*N_8_1_1_1)-
    ((alpha_out + mu_3_1_1 + pi_3_4 + kappa_3_1_1 + (theta_1*varpi_g*pi_3_6) + eta_1_2_1)*N_3_2_1_1)
  
  })
}
