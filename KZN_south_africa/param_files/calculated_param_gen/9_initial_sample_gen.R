#draws initial samples using latin hypercube

#clean workspace
rm(list = ls())
gc()

library(lhs)

set.seed(1)
n_samples <- 10000
n_params <- 37

#update parameter file
outdir <- paste0(here(),'/param_files/')

raw_lhs_df <- as.data.frame(randomLHS(n_samples, n_params))
calib_param_names<-c('beta_1', #effective contact rate by g
                     'beta_2',
                     'phi_2', #relative transmisibility by h
                     'phi_3',
                     'phi_4',
                     'varepsilon', #fraction MDR
                     'iota_2', #reduced risk infection on IPT by dr
                     'xi', #increased risk of infection after active TB
                     'zeta', #reduced risk of reinfection after LTBI
                     'pi_34', #recent infection period
                     'pi_36', #base rates of TB prog
                     'pi_46', #base prog to active TB
                     'pi_56',
                     'pi_86',
                     'theta_2', #relative risk of prog by HIV
                     'theta_3',
                     'theta_4',
                     'upsilon_1', #recovery rate by HIV
                     'upsilon_2',
                     'upsilon_3',
                     'upsilon_4',
                     'pi_76', #TB relapse
                     'prop_eligible_males', #prop CD4 more eligible
                     'prop_eligible_females',
                     'eta_12_1', #incidence rate
                     'eta_12_2',
                     'eta_23_1', #prog of CD4 counts
                     'eta_23_2',
                     'mort_base_1',
                     'mort_base_2',
                     'increased_mort_risk_hiv_only_2',
                     'increased_mort_risk_hiv_only_3',
                     'increased_mort_risk_hiv_only_4',
                     'increased_mort_risk_TB_hiv_1',
                     'increased_mort_risk_TB_hiv_2',
                     'increased_mort_risk_TB_hiv_3',
                     'increased_mort_risk_TB_hiv_4'
                     )

colnames(raw_lhs_df)<-calib_param_names
raw_lhs_df$sim_id<-1:nrow(raw_lhs_df)

setwd(outdir)
write.csv(raw_lhs_df, 'sim_calibration_ref_df.csv', row.names = FALSE)
