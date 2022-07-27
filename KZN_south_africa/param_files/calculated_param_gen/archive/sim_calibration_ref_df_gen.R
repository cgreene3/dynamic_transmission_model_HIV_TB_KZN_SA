#create sim id ref df based on set of calibrated values

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('readxl', 'here', 'dplyr', 'reshape2', 'ggplot2', 'stringr', 'normalr'), require, character.only=T)

#update parameter file
outdir <- paste0(here(),'/param_files/')

#######Create calibration ref df#########

#mort values to test
#hiv only
h2_test<-c(8,12)
h3_test<-c(35,45)
h4_test<-c(1.2,1.5)

#tb/hiv coinfection
tb_h1_test<-c(14,17)
tb_h2_test<-c(22,30)
tb_h3_test<-c(40,60)
tb_h4_test<-c(17,20)

#betas to test
beta_1_test<-seq(from = 6, to = 12, by = 1)
beta_2_test<-seq(from = 6, to = 12, by = 1)

#create a sim id for each combination
sim_calib_id<-c()
b1_contact<-c()
b2_contact<-c()
hiv2_mort<-c()
hiv3_mort<-c()
hiv4_mort<-c()
tb_hiv1_mort<-c()
tb_hiv2_mort<-c()
tb_hiv3_mort<-c()
tb_hiv4_mort<-c()

sim_id_temp<-1

for (b1 in beta_1_test){
  for (b2 in beta_2_test){
    for (hiv2 in h2_test){
      for (hiv3 in h3_test){
        for (hiv4 in h4_test){
          for (tb_hiv1 in tb_h1_test){
            for (tb_hiv2 in tb_h2_test){
              for (tb_hiv3 in tb_h3_test){
                for (tb_hiv4 in tb_h4_test){
                  
                  sim_calib_id<-c(sim_calib_id, sim_id_temp)
                  b1_contact<-c(b1_contact, b1)
                  b2_contact<-c(b2_contact, b2)
                  hiv2_mort<-c(hiv2_mort, hiv2)
                  hiv3_mort<-c(hiv3_mort, hiv3)
                  hiv4_mort<-c(hiv4_mort, hiv4)
                  tb_hiv1_mort<-c(tb_hiv1_mort, tb_hiv1)
                  tb_hiv2_mort<-c(tb_hiv2_mort, tb_hiv2)
                  tb_hiv3_mort<-c(tb_hiv3_mort, tb_hiv3)
                  tb_hiv4_mort<-c(tb_hiv4_mort, tb_hiv4)
                  sim_id_temp<-sim_id_temp+1
                }
              }
            }
          }
        }
      }
    }
  }
}

sim_calibration_ref_df<-data.frame(sim_calib_id, b1_contact, b2_contact,
                                   hiv2_mort, hiv3_mort, hiv4_mort,
                                   tb_hiv1_mort, tb_hiv2_mort,
                                   tb_hiv3_mort, tb_hiv4_mort)

setwd(outdir)
write.csv(sim_calibration_ref_df, 'sim_calibration_ref_df.csv', row.names = FALSE)
