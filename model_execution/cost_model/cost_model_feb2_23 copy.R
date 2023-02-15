#calculates costs for each program

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('dplyr', 'deSolve', 
         'readxl', 'stringr', 
         'reshape2', 'ggplot2', 
         'varhandle', 'here', 'extrafont', 'RColorBrewer'), require, character.only=T)

library(here)
indir_outputs<-paste0(here(), '/results/economic_analysis/program_outputs_economic_analysis')
indir_input_cost_data<-paste0(here(), '/param_files/cost_parameters')
indir_disability_weight_data<-paste0(here(), '/param_files/disability_weights')
indir_delay_data<-paste0(here(), '/results/calibration_analysis')
outdir_results<-paste0(here(), '/results/ICER')

####Define Sets and Set ID location ####

#define sets
TB_SET<-1:8
DR_SET<-1:2
HIV_SET<-1:4
G_SET<-1:2

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

#######Read in data#########
#read in cost data
setwd(indir_input_cost_data)
input_costs_ref_df<-read.csv('input_costs_ref_df.csv')
cost_sensitivity_scenarios_df<-read.csv('sensitivity_scenarios.csv')

#read in disability weights
setwd(indir_disability_weight_data)
disability_weights_ref_df<-read.csv('disability_weights_df.csv')

#read in sets for TB durations
setwd(indir_delay_data)
durations_ref_df<-read.csv('accepted_calibration_sets_ref_df.csv')%>%
  select(c('sim_id',
           'omega_',
           'pi_67.', 
           'upsilon_1.'	,
           'upsilon_2.'	,
           'upsilon_3.',
           'upsilon_4.'))

#discount ref
#discount factor
R = 0.03
eval_start_yr <- 2018

disc_factor_ref_df<-data.frame(year = 2018:2027)
disc_factor_ref_df<-disc_factor_ref_df%>%
  mutate(FY = 1/(1+R)^(year - eval_start_yr))

##combine simulation estimates##
setwd(indir_outputs)
all_files_in_outdir <- list.files(pattern="*.csv")

#create reference data frame for each sim id and each program
cost_total_df<-data.frame() #total costs by sim id and program and sensitivity id (0 is average)
cost_itemize_df<-data.frame() #costs by line item (IPT, ART, TB care, O HIV care)
DALY_total_df<-data.frame() #YLL, YLD, DALY
mort_itemize_df<-data.frame() #TB mortality by HIV status and gender and Omort by gender
mort_total_df<-data.frame() #total mortality
inc_total_df<-data.frame() #total inc
n_ipt_df<-data.frame() #nipt

#######read in model outputs and calculate costs, TB cases, deaths and DALYs for each sim####
for (i in 1:length(all_files_in_outdir)){
  
  print(i)
  #pull sim id of current for loop eval#
  sim_id_temp<-str_split(all_files_in_outdir[i], '_')
  sim_id_temp<-str_split(sim_id_temp[[1]][7], '.csv')
  sim_id_temp<-as.integer(sim_id_temp[[1]][1])
  
  setwd(indir_outputs)
  model_outputs_temp_df<-read.csv(all_files_in_outdir[i])%>%
    filter(year >= 2018)%>%
    filter(year <= 2027)%>%
    filter(time != 78) #filter out starting conditions
  
  #create durations data frame to adjust costs#
  durations_ref_df_filter<-durations_ref_df%>%
    filter(sim_id == sim_id_temp)
  
  durations_ref_df_temp<-data.frame(duration_adj_desc = as.character(),
                                duration_adj = as.numeric())
  
  #ocare from year to monthly costs
  durations_ref_df_temp<-rbind(durations_ref_df_temp,
                           data.frame(duration_adj_desc = 'yearly',
                                      duration_adj = 1/12))
  
  #ipt costs from yearly to six months
  ipt_duration_temp<-(durations_ref_df_filter$omega_[1])
  durations_ref_df_temp<-rbind(durations_ref_df_temp,
                           data.frame(duration_adj_desc = 'ipt_duration_adj',
                                      duration_adj = ipt_duration_temp/12))
  
  #active TB from yearly costs to 1/(upsilon*pi_{6,7})
  tb_h1_duration_temp<-(durations_ref_df_filter$upsilon_1.[1]*durations_ref_df_filter$pi_67.[1])
  tb_h2_duration_temp<-(durations_ref_df_filter$upsilon_2.[1]*durations_ref_df_filter$pi_67.[1])
  tb_h3_duration_temp<-(durations_ref_df_filter$upsilon_3.[1]*durations_ref_df_filter$pi_67.[1])
  tb_h4_duration_temp<-(durations_ref_df_filter$upsilon_4.[1]*durations_ref_df_filter$pi_67.[1])
  
  durations_ref_df_temp<-rbind(durations_ref_df_temp,
                           data.frame(duration_adj_desc = 'tb_active_duration_adj_h1',
                                      duration_adj = tb_h1_duration_temp))
  durations_ref_df_temp<-rbind(durations_ref_df_temp,
                           data.frame(duration_adj_desc = 'tb_active_duration_adj_h2',
                                      duration_adj = tb_h2_duration_temp))
  durations_ref_df_temp<-rbind(durations_ref_df_temp,
                           data.frame(duration_adj_desc = 'tb_active_duration_adj_h3',
                                      duration_adj = tb_h3_duration_temp))
  durations_ref_df_temp<-rbind(durations_ref_df_temp,
                           data.frame(duration_adj_desc = 'tb_active_duration_adj_h4',
                                      duration_adj = tb_h4_duration_temp))
  
  
  #calculate the number of individuals in each cost category from model outputs#
  
  #ocare h
  model_outputs_temp_df$n_other_h2<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[TB_SET, DR_SET, 2, G_SET])])
  model_outputs_temp_df$n_other_h3<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[TB_SET, DR_SET, 3, G_SET])])
  model_outputs_temp_df$n_other_h4<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[TB_SET, DR_SET, 4, G_SET])])
  
  #n on IPT
  model_outputs_temp_df$n_ipt<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[c(2,5), DR_SET, 4, G_SET])])
  
  #TB care by DR status and HIV status (due to varying durations)
  #DS
  model_outputs_temp_df$n_tb_r1_h1<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 1, 1, G_SET])])
  model_outputs_temp_df$n_tb_r1_h2<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 1, 2, G_SET])])
  model_outputs_temp_df$n_tb_r1_h3<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 1, 3, G_SET])])
  model_outputs_temp_df$n_tb_r1_h4<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 1, 4, G_SET])])
  
  #MDR
  model_outputs_temp_df$n_tb_r2_h1<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 2, 1, G_SET])])
  model_outputs_temp_df$n_tb_r2_h2<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 2, 2, G_SET])])
  model_outputs_temp_df$n_tb_r2_h3<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 2, 3, G_SET])])
  model_outputs_temp_df$n_tb_r2_h4<-rowSums(model_outputs_temp_df[,4+c(N_t_r_h_g_ref[6, 2, 4, G_SET])])
  
  ###subset model outputs used to generate costs #####
  cost_model_temp_df<-model_outputs_temp_df%>%
    select(c('sim_id', 'year', 
             'time', 'program_id', 
             'n_other_h2',
             'n_other_h3',
             'n_other_h4',
             'n_ipt',
             'n_tb_r1_h1',
             'n_tb_r1_h2',
             'n_tb_r1_h3',
             'n_tb_r1_h4',
             'n_tb_r2_h1',
             'n_tb_r2_h2',
             'n_tb_r2_h3',
             'n_tb_r2_h4'))
  
  
  #join model outputs for cost model with duration, discount factor, and input cost data#
  cost_model_temp_df_melt <- melt(cost_model_temp_df, id = c("sim_id", "year", "time", "program_id"))
  cost_model_temp_df_melt<-cost_model_temp_df_melt%>%
    mutate(variable = if_else(variable == 'n_other_h4', 
                              if_else(program_id == 1, 
                                      'n_other_h4_facility', 'n_other_h4_community'), 
                              as.character(variable)))

  #average
  cost_model_temp_df_join_avg<-cost_model_temp_df_melt%>%
    full_join(input_costs_ref_df, by = c('variable'))%>%
    left_join(durations_ref_df_temp, by = c('duration_adj_desc'))%>%
    left_join(disc_factor_ref_df, by = c('year'))%>%
    mutate(undisc_val_cost_duration_adj = cost*duration_adj*value,
           disc_val_cost_duration_adj = cost*duration_adj*value*FY)
  
  #itemize average cost by program, line item and year #
  #cost_itemize_df_temp_by_year_by_sim<-cost_model_temp_df_join%>%
  #  group_by(sim_id, program_id, cost_desc, year)%>%
  #  summarise(undisc_cost_itemized = sum(undisc_val_cost_duration_adj),
  #            disc_cost_itemized = sum(disc_val_cost_duration_adj))
  
  #cost_itemize_df_temp<-cost_itemize_df_temp_by_year_by_sim%>%
  #  group_by(sim_id, program_id, cost_desc)%>%
  #  summarise(undisc_cost_itemized = sum(undisc_cost_itemized),
  #            disc_cost_itemized = sum(disc_cost_itemized))
  
  #summarise total costs by program ##
  cost_total_df_temp_avg<-cost_model_temp_df_join_avg%>%
    group_by(sim_id, program_id)%>%
    summarise(undisc_cost_total = sum(undisc_val_cost_duration_adj),
              disc_cost_total = sum(disc_val_cost_duration_adj))
  
  cost_total_df_temp_avg$sensitivity_id<-0
  
  if (i == 1){
    cost_total_df<-cost_total_df_temp_avg
    #cost_itemize_df<-cost_itemize_df_temp
  } else {
    cost_total_df<-rbind(cost_total_df, cost_total_df_temp_avg)
    #cost_itemize_df<-rbind(cost_itemize_df, cost_itemize_df_temp)
  }
  
  
  #run cost sensitivity scenarios
  for (s in unique(cost_sensitivity_scenarios_df$sensitivity_id)){
    cost_sensitivity_scenarios_df_current<-cost_sensitivity_scenarios_df%>%
      filter(sensitivity_id == s)
    cost_desc_sensitivity_temp<-as.character(cost_sensitivity_scenarios_df_current$cost_desc[1])
    sensitivity_desc_temp<-as.character(cost_sensitivity_scenarios_df_current$sensitivity_desc[1])
    cost_adj_temp<-if_else(sensitivity_desc_temp == "high", 1.25, 0.75)
    
    input_costs_ref_df_sensitivity<-input_costs_ref_df%>%
      mutate(cost_adj = if_else(cost_desc == cost_desc_sensitivity_temp,
                                cost_adj_temp, 1))%>%
      mutate(cost = cost*cost_adj)%>%
      select(-c('cost_adj'))
    
    cost_model_temp_df_join_sensitivity_temp<-cost_model_temp_df_melt%>%
      full_join(input_costs_ref_df_sensitivity, by = c('variable'))%>%
      left_join(durations_ref_df_temp, by = c('duration_adj_desc'))%>%
      left_join(disc_factor_ref_df, by = c('year'))%>%
      mutate(undisc_val_cost_duration_adj = cost*duration_adj*value,
             disc_val_cost_duration_adj = cost*duration_adj*value*FY)
    
    #summarise total costs by program #
    cost_total_df_temp_sensitivity<-cost_model_temp_df_join_sensitivity_temp%>%
      group_by(sim_id, program_id)%>%
      summarise(undisc_cost_total = sum(undisc_val_cost_duration_adj),
                disc_cost_total = sum(disc_val_cost_duration_adj))
    
    cost_total_df_temp_sensitivity$sensitivity_id <- s
    
    cost_total_df<-rbind(cost_total_df, cost_total_df_temp_sensitivity)
    
  }
  
  ###pull n_ipt overtime###
  n_ipt_temp_df<-model_outputs_temp_df%>%
    select(c('sim_id', 'year', 
             'time', 'program_id', 
             'n_ipt'))
  
  if (i == 1){
    n_ipt_df<-n_ipt_temp_df
  } else {
    n_ipt_df<-rbind(n_ipt_df, n_ipt_temp_df)
  }
  
  ####DALYs###
  setwd(indir_outputs)
  model_outputs_temp_df<-read.csv(all_files_in_outdir[i])%>%
    filter(year >= 2018)%>%
    filter(year <= 2027)%>%
    filter(time != 78) #filter out starting conditions
  
  #calculate years lived with disability#
  YLD_temp_df_melt<-melt(model_outputs_temp_df, id = c("year", "sim_id", "time", "program_id"))
  YLD_temp_df_melt<-rename(YLD_temp_df_melt, compartment_id = variable)
  
  YLD_temp_df_calc<-YLD_temp_df_melt%>%
    left_join(disability_weights_ref_df, by = c('compartment_id'))%>%
    left_join(disc_factor_ref_df, by = c('year'))%>%
    filter(!is.na(TB_compartment))%>%
    mutate(undisc_disability = value * (1/12) * weight_value,
           disc_disability =  value * (1/12) * weight_value * FY)
  
  YLD_tot_df_temp<-YLD_temp_df_calc%>%
    group_by(sim_id, program_id)%>%
    summarise(undisc_YLD_total = sum(undisc_disability),
              disc_YLD_total = sum(disc_disability))
  
  #calculate years of life lost#
  model_outputs_temp_df<-read.csv(all_files_in_outdir[i])%>%
    filter(time <= 10) #filter out warm up period 
  
  YLL_df_temp<-model_outputs_temp_df%>%
    mutate(Tb_mort_neg_female_t_per_100k_ppl = Tb_mort_neg_female_cumulative-lag(Tb_mort_neg_female_cumulative),
           Tb_mort_neg_male_t_per_100k_ppl = Tb_mort_neg_male_cumulative-lag(Tb_mort_neg_male_cumulative),
           Tb_mort_pos_female_t_per_100k_ppl = Tb_mort_pos_female_cumulative-lag(Tb_mort_pos_female_cumulative),
           Tb_mort_pos_male_t_per_100k_ppl = Tb_mort_pos_male_cumulative-lag(Tb_mort_pos_male_cumulative),
           O_mort_male_t_per_100k_ppl = O_mort_male_cumulative-lag(O_mort_male_cumulative),
           O_mort_female_t_per_100k_ppl = O_mort_female_cumulative-lag(O_mort_female_cumulative))%>%
    filter(time > 0)%>% #filter out starting conditions after calculating january 2018 mortality
    select(c("sim_id", 'year', "time", "program_id",
             "Tb_mort_neg_male_t_per_100k_ppl", #
             "Tb_mort_neg_female_t_per_100k_ppl",#
             "Tb_mort_pos_male_t_per_100k_ppl",#
             "Tb_mort_pos_female_t_per_100k_ppl",#
             "O_mort_male_t_per_100k_ppl", 
             "O_mort_female_t_per_100k_ppl"))
  
  
  YLL_df_temp_melt<-melt(YLL_df_temp, 
                         id = c("year", "sim_id", "time", "program_id"))
  
  YLL_df_temp_summarised<-YLL_df_temp_melt%>%
    group_by(sim_id, year, time, program_id)%>%
    summarise(mortality_total = sum(value))%>%
    ungroup()%>%
    mutate(time = time + 2018 - 0.083333)%>% ###to convert time to 2018.083333 = January 2018, etc.
    mutate(year = floor(time))%>% #match year with converted time for discounting.
    mutate(YLL_total_time = 2028 - time)%>%
    mutate(YLL_time_n_mort = YLL_total_time*mortality_total)%>%
    left_join(disc_factor_ref_df, by = c('year'))%>%
    mutate(undisc_YLL = YLL_time_n_mort,
           disc_YLL = YLL_time_n_mort*FY)%>%
    group_by(sim_id, program_id)%>%
    summarise(undisc_YLL_total = sum(undisc_YLL),
              disc_YLL_total  = sum(disc_YLL))
  
  DALY_total_df_temp<-YLL_df_temp_summarised%>%
    left_join(YLD_tot_df_temp, by = c("sim_id", "program_id"))%>%
    mutate(undisc_DALY = undisc_YLL_total + undisc_YLD_total,
           disc_DALY = disc_YLL_total + disc_YLD_total)
 
  if (i == 1){
    DALY_total_df<-DALY_total_df_temp
  } else {
    DALY_total_df<-rbind(DALY_total_df, DALY_total_df_temp)
  }
  
  #get mortality by group
  mort_itemize_df_temp<-YLL_df_temp_melt%>%
    group_by(sim_id, program_id, variable)%>%
    summarise(total = sum(value))
  
  if (i == 1){
    mort_itemize_df<-mort_itemize_df_temp
  } else {
    mort_itemize_df<-rbind(mort_itemize_df, mort_itemize_df_temp)
  }
  
  #get total mortality
  mort_total_df_temp<-YLL_df_temp_melt%>% #YLL df has mort values.
    filter(!grepl('O_mort', variable))%>%
    mutate(time = time + 2018 - 0.083333)%>% ###to convert time to 2018.083333 = January 2018, etc.
    mutate(year = floor(time))%>% #match year with converted time for discounting.
    left_join(disc_factor_ref_df, by = "year")%>%
    mutate(disc_mort_value = value * FY)%>%
    group_by(sim_id, program_id)%>%
    summarise(undisc_mort = sum(value),
              disc_mort = sum(disc_mort_value))
  
  if (i == 1){
    mort_total_df<-mort_total_df_temp
  } else {
    mort_total_df<-rbind(mort_total_df, mort_total_df_temp)
  }
  
  #pull incidence
  inc_total_df_temp<-model_outputs_temp_df%>%
    select(c("year", "time", "sim_id", "program_id",
             "Tb_inc_neg_male_cumulative", 
             "Tb_inc_neg_female_cumulative", 
             "Tb_inc_pos_male_cumulative",    
             "Tb_inc_pos_female_cumulative"))%>%
    mutate(Tb_inc_neg_female_t_per_100k_ppl = Tb_inc_neg_female_cumulative-lag(Tb_inc_neg_female_cumulative),
           Tb_inc_neg_male_t_per_100k_ppl = Tb_inc_neg_male_cumulative-lag(Tb_inc_neg_male_cumulative),
           Tb_inc_pos_female_t_per_100k_ppl = Tb_inc_pos_female_cumulative-lag(Tb_inc_pos_female_cumulative),
           Tb_inc_pos_male_t_per_100k_ppl = Tb_inc_pos_male_cumulative-lag(Tb_inc_pos_male_cumulative))%>%
    filter(time > 0)%>% #filter out starting conditions after calculating january 2018 mortality
    mutate(Tb_inc_total = Tb_inc_neg_female_t_per_100k_ppl+
             Tb_inc_neg_male_t_per_100k_ppl+
             Tb_inc_pos_female_t_per_100k_ppl+
             Tb_inc_pos_male_t_per_100k_ppl)%>%
    select(c("sim_id", 'year', "time", "program_id",
             "Tb_inc_total"))%>%
    mutate(time = time + 2018 - 0.083333)%>% ###to convert time to 2018.083333 = January 2018, etc.
    mutate(year = floor(time))%>% #match year with converted time for discounting.
    left_join(disc_factor_ref_df, by = "year")%>%
    mutate(disc_inc = Tb_inc_total * FY)%>%
    group_by(sim_id, program_id)%>%
    summarise(undisc_inc = sum(Tb_inc_total),
              disc_inc = sum(disc_inc))
  
  if (i == 1){
    inc_total_df<-inc_total_df_temp
  } else {
    inc_total_df<-rbind(inc_total_df, inc_total_df_temp)
  }
  
  
}

######SUMMARISED OVER ALL SIMULATIONS######

######program, incremental health gains, costs, and ICERs######

incremental_ICER_calcs_df<-cost_total_df%>%
  left_join(DALY_total_df, by = c('sim_id', 
                                  'program_id'))%>%
  left_join(mort_total_df, by = c("sim_id", "program_id"))%>%
  left_join(inc_total_df, by = c("sim_id", "program_id"))

incremental_ICER_calcs_df_melt<-melt(incremental_ICER_calcs_df, 
                             id = c("sim_id", "program_id", "sensitivity_id"))%>%
  filter(!grepl("YLL", variable))%>%
  filter(!grepl("YLD", variable))

incremental_ICER_calcs_df_melt$program_id<-paste0('p', incremental_ICER_calcs_df_melt$program_id) 
incremental_ICER_calcs_df_melt$program_variable <- paste0(incremental_ICER_calcs_df_melt$program_id, 
                                              "_", incremental_ICER_calcs_df_melt$variable)

incremental_ICER_calcs_df_melt<-incremental_ICER_calcs_df_melt%>%
  select(-c('program_id', 'variable'))

incremental_ICER_calcs_df_cast<-dcast(incremental_ICER_calcs_df_melt, sim_id + sensitivity_id ~ program_variable)


#######incremental health gains, costs, and ICERs EXPECTED COST #######
incremental_ICER_calcs_df_exp_cost<-incremental_ICER_calcs_df_cast%>%
  filter(sensitivity_id == 0)%>%
  select(-c('sensitivity_id'))%>%
  mutate(p2_p1_disc_incremental_cost = p2_disc_cost_total-p1_disc_cost_total,
         p3_p1_disc_incremental_cost = p3_disc_cost_total-p1_disc_cost_total,
         p3_p2_disc_incremental_cost = p3_disc_cost_total-p2_disc_cost_total,
         p2_p1_disc_incremental_DALY = p1_disc_DALY-p2_disc_DALY,
         p3_p1_disc_incremental_DALY = p1_disc_DALY-p3_disc_DALY,
         p3_p2_disc_incremental_DALY = p2_disc_DALY-p3_disc_DALY,
         p2_p1_undisc_incremental_cost = p2_undisc_cost_total-p1_undisc_cost_total,
         p3_p1_undisc_incremental_cost = p3_undisc_cost_total-p1_undisc_cost_total,
         p3_p2_undisc_incremental_cost = p3_undisc_cost_total-p2_undisc_cost_total,
         p2_p1_undisc_incremental_DALY = p1_undisc_DALY-p2_undisc_DALY,
         p3_p1_undisc_incremental_DALY = p1_undisc_DALY-p3_undisc_DALY,
         p3_p2_undisc_incremental_DALY = p2_undisc_DALY-p3_undisc_DALY)%>%
  mutate(p2_p1_disc_ICER = (p2_p1_disc_incremental_cost/(p2_p1_disc_incremental_DALY)),
         p3_p1_disc_ICER = (p3_p1_disc_incremental_cost/(p3_p1_disc_incremental_DALY)),
         p3_p2_disc_ICER = (p3_p2_disc_incremental_cost/(p3_p2_disc_incremental_DALY)),
         p2_p1_undisc_ICER = (p2_p1_undisc_incremental_cost/(p2_p1_undisc_incremental_DALY)),
         p3_p1_undisc_ICER = (p3_p1_undisc_incremental_cost/(p3_p1_undisc_incremental_DALY)),
         p3_p2_undisc_ICER = (p3_p2_undisc_incremental_cost/(p3_p2_undisc_incremental_DALY)))%>%
  mutate(p3_p1_disc_inc = p1_disc_inc - p3_disc_inc,
         p2_p1_disc_inc = p1_disc_inc - p2_disc_inc,
         p3_p2_disc_inc = p2_disc_inc - p3_disc_inc,
         p3_p1_undisc_inc = p1_undisc_inc - p3_undisc_inc,
         p2_p1_undisc_inc = p1_undisc_inc - p2_undisc_inc,
         p3_p2_undisc_inc = p2_undisc_inc - p3_undisc_inc)%>%
  mutate(p3_p1_disc_inc_per_dollar = p3_p1_disc_incremental_cost/p3_p1_disc_inc,
         p2_p1_disc_inc_per_dollar = p2_p1_disc_incremental_cost/p2_p1_disc_inc,
         p3_p2_disc_inc_per_dollar = p3_p2_disc_incremental_cost/p3_p2_disc_inc,
         p3_p1_undisc_inc_per_dollar = p3_p1_undisc_incremental_cost/p3_p1_undisc_inc,
         p2_p1_undisc_inc_per_dollar = p2_p1_undisc_incremental_cost/p2_p1_undisc_inc,
         p3_p2_undisc_inc_per_dollar = p3_p2_undisc_incremental_cost/p3_p2_undisc_inc)%>%
mutate(p3_p1_disc_mort = p1_disc_mort - p3_disc_mort,
       p2_p1_disc_mort = p1_disc_mort - p2_disc_mort,
       p3_p2_disc_mort = p2_disc_mort - p3_disc_mort,
       p3_p1_undisc_mort = p1_undisc_mort - p3_undisc_mort,
       p2_p1_undisc_mort = p1_undisc_mort - p2_undisc_mort,
       p3_p2_undisc_mort = p2_undisc_mort - p3_undisc_mort)%>%
  mutate(p3_p1_disc_mort_per_dollar = p3_p1_disc_incremental_cost/p3_p1_disc_mort,
         p2_p1_disc_mort_per_dollar = p2_p1_disc_incremental_cost/p2_p1_disc_mort,
         p3_p2_disc_mort_per_dollar = p3_p2_disc_incremental_cost/p3_p2_disc_mort,
         p3_p1_undisc_mort_per_dollar = p3_p1_undisc_incremental_cost/p3_p1_undisc_mort,
         p2_p1_undisc_mort_per_dollar = p2_p1_undisc_incremental_cost/p2_p1_undisc_mort,
         p3_p2_undisc_mort_per_dollar = p3_p2_undisc_incremental_cost/p3_p2_undisc_mort)
  


incremental_ICER_calcs_df_exp_cost_melt<-melt(incremental_ICER_calcs_df_exp_cost,
                                              id = c("sim_id"))

incremental_ICER_calcs_df_exp_cost_mean_max_min<-incremental_ICER_calcs_df_exp_cost_melt%>%
  group_by(variable)%>%
  summarise(mean = mean(value),
            min = min(value),
            max = max(value))


######ICERs SENSITIVITY COST #######
#comparing programs 1 and 3 only
ICERS_sensitivity_cost_df<-incremental_ICER_calcs_df_cast%>%
  select(c("sim_id", "sensitivity_id",
           "p1_disc_DALY", "p3_disc_DALY",
           "p1_disc_cost_total", "p3_disc_cost_total"))%>%
  mutate(incremental_disc_DALY = p1_disc_DALY - p3_disc_DALY,
         incremental_disc_cost = p3_disc_cost_total - p1_disc_cost_total)%>%
  mutate(p3_p1_undisc_ICER = incremental_disc_cost/incremental_disc_DALY)%>%
  group_by(sensitivity_id)%>%
  summarise(est = mean(p3_p1_undisc_ICER))

ICERS_sensitivity_cost_df<-ICERS_sensitivity_cost_df%>%
  left_join(cost_sensitivity_scenarios_df, by = c("sensitivity_id"))

avg_ICER<-filter(ICERS_sensitivity_cost_df, sensitivity_id == 0)

ICERS_sensitivity_cost_df2<-filter(ICERS_sensitivity_cost_df, 
                                   sensitivity_id != 0)

test<-ICERS_sensitivity_cost_df2%>%
  select(c('cost_desc', 'sensitivity_desc', 'est'))

test2<-dcast(test,
             cost_desc ~ sensitivity_desc)

test2$UI<-abs(test2$high-test2$low)


test2<-test2%>%
  mutate(UI_rank = rank(UI))%>%
  arrange(UI_rank)

width <- 0.95
type<-test2$cost_desc%>% unlist() %>% levels()


tornado_plot_df<-test2%>%
  mutate(xmin=as.numeric(UI_rank)-width/2,
         xmax=as.numeric(UI_rank)+width/2)%>%
  select(-c('UI', 'UI_rank'))

tornado_plot_df_melt<-melt(tornado_plot_df, id = c("cost_desc",
                                                   "xmin", "xmax"))

tornado_plot_df_melt$avg_ICER <-avg_ICER$est

tornado_plot_df2<-tornado_plot_df_melt%>%
  ungroup()%>%
  mutate(ymin = pmin(value, avg_ICER),
         ymax = pmax(value, avg_ICER))
  
ggplot() + 
  geom_rect(data = tornado_plot_df2,
            aes(ymax=ymax, 
                ymin=ymin, 
                xmax=xmax, 
                xmin=xmin, 
                fill=variable))+
  theme_bw() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = avg_ICER$est)+
  geom_hline(yintercept = 750, linetype=2, color = "grey")+
  scale_x_continuous(breaks = c(1:length(tornado_plot_df2$cost_desc)), 
                     labels = tornado_plot_df2$cost_desc) +
  coord_flip()

#nipt
n_ipt_df_by_community_v_facility<-n_ipt_df%>%
  mutate(program_group = if_else(program_id == 3, "community_based_ipt", "facility_based_ipt"))%>%
  mutate(value_avg = if_else(program_id == 3, n_ipt, n_ipt/2))%>%
  group_by(sim_id, program_group)%>%
  summarise(n_ipt_total = sum(value_avg))%>%
  group_by(program_group)%>%
  summarise(mean = round(mean(n_ipt_total)),
            min = round(min(n_ipt_total)),
            max = round(max(n_ipt_total)))

write.csv(n_ipt_df,
          'n_ipt_df.csv',
          row.names = FALSE)

##if want to print active TB durations
active_tb_durations_summary_df<-durations_ref_df%>%
  mutate(hiv1 = 1/(upsilon_1.*2),
         hiv2 = 1/(upsilon_2.*2),
         hiv3 = 1/(upsilon_3.*2),
         hiv4 = 1/(upsilon_4.*2))%>%
  select(-c('omega_',
            'pi_67.', 
            'upsilon_1.'	,
            'upsilon_2.'	,
            'upsilon_3.',
            'upsilon_4.'))

active_tb_durations_summary_df_melt<-melt(active_tb_durations_summary_df, id = "sim_id")
active_tb_durations_summary_df2<-active_tb_durations_summary_df_melt%>%
  group_by(variable)%>%
  summarise(mean = round(mean(value)*12, 1),
            max = round(max(value)*12, 1),
            min = round(min(value)*12, 1))

write.csv(active_tb_durations_summary_df2,
          'active_tb_calibrated_durations.csv',
          row.names = FALSE)

png(file="ICER_threshold.png")#,
    #width=600, height=500)
print(p)
dev.off()

write.csv(perc_cost_effective,
          'perc_cost_effective.csv',
          row.names = FALSE)
