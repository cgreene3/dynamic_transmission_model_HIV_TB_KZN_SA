#Cost calculations for community IPT analysis
#November 2021
#Updated in December 2022 to reflect monthly costs instead of annual.
#Note that I need to specify which version of the input costs to use by the date of the csv.

rm(list=ls())


library(dplyr)
library(reshape2)

home_dir <- ("H:/Research Projects/Ross K01 transmission model/")
in_dir <- (paste0(home_dir, "Data/"))
out_dir <- (paste0(home_dir, "Analysis/"))

#Read datasets
cost_date<-'20221202'
cc<-read.csv(file=paste0(in_dir, "Costs/cost_components_mo_", cost_date, ".csv"))
df<-read.csv(file=paste0(in_dir, "GBD disability weights/disability_weights_matched.csv"))

df$cost<-0

#Cost models for HIV. Each HIV compartment includes outpatient and inpatient HIV-related (non-TB) costs.
#HIV4 compartment includes policy_1 scenario for facility-based care and policy_2 and 3 scenarios for community-based care.
hiv2<-sum(cc$mo_cost[which(cc$var_name=="outpt_hiv")], cc$mo_cost[which(cc$var_name=="inpt_hiv_cd4_200+")])
hiv3<-sum(cc$mo_cost[which(cc$var_name=="outpt_hiv")], cc$mo_cost[which(cc$var_name=="inpt_hiv_cd4_0200")])
hiv4_prog1<-sum(cc$mo_cost[which(cc$var_name=="fac_ART")], cc$mo_cost[which(cc$var_name=="inpt_hiv_art")])
hiv4_prog2<-sum(cc$mo_cost[which(cc$var_name=="com_art")], cc$mo_cost[which(cc$var_name=="inpt_hiv_art")])

#Cost models for IPT. DILI stands for drug-induced liver injury, which is the primary toxicity from IPT.
dili_cost<-sum(cc$mo_cost[which(cc$var_name=="ipt_provider_dili")], cc$mo_cost[which(cc$var_name=="ipt_labs_dili")])
dili_prob<-cc$mo_cost[which(cc$var_name=="ipt_prob_dili")]
on_ipt<-sum(cc$mo_cost[which(cc$var_name=="ipt_med")], cc$mo_cost[which(cc$var_name=="ipt_provider")], dili_cost*dili_prob)

#Cost models for treatment of active TB (DS and MDR)
#The total cost of TB treatment for HIV neg and HIV pos is the same, but the variables differ because the duration of active TB differs by HIV status. The cost of TB treatment is divided over different durations of active TB.
tb6_ds_hiv_neg<-cc$mo_cost[which(cc$var_name=="tb_ds_hiv_neg")]
tb6_ds_hiv_pos<-cc$mo_cost[which(cc$var_name=="tb_ds_hiv_pos")]
tb6_mdr_hiv_neg<-cc$mo_cost[which(cc$var_name=="tb_mdr_hiv_neg")]
tb6_mdr_hiv_pos<-cc$mo_cost[which(cc$var_name=="tb_mdr_hiv_pos")]

#Matching cost models to annual compartments 
#Need to have a separate dataset for care programs 2 and 3 that incorporates the costs of community-based ART

#Start with HIV costs
df$cost[df$HIV==2]<-hiv2
df$cost[df$HIV==3]<-hiv3
df$cost[df$HIV==4]<-hiv4_prog1

#Add IPT costs to TB compartments 2 and 5 with corresponding HIV costs
df$cost[df$HIV==2 & (df$TB==2|df$TB==5)]<-hiv2+on_ipt
df$cost[df$HIV==3 & (df$TB==2|df$TB==5)]<-hiv3+on_ipt
df$cost[df$HIV==4 & (df$TB==2|df$TB==5)]<-hiv4_prog1+on_ipt

#Code TB costs by HIV and MDR status and replace earlier values for people with TB=6
#DS drug-susceptible active TB
df$cost[df$TB==6 & df$DR==1 & df$HIV==1] <-tb6_ds_hiv_neg
df$cost[df$TB==6 & df$DR==1 & df$HIV==2] <-tb6_ds_hiv_pos+hiv2
df$cost[df$TB==6 & df$DR==1 & df$HIV==3] <-tb6_ds_hiv_pos+hiv3
df$cost[df$TB==6 & df$DR==1 & df$HIV==4] <-tb6_ds_hiv_pos+hiv4_prog1
#MDR multi-drug resistant active TB
df$cost[df$TB==6 & df$DR==2 & df$HIV==1] <-tb6_mdr_hiv_neg
df$cost[df$TB==6 & df$DR==2 & df$HIV==2] <-tb6_mdr_hiv_pos+hiv2
df$cost[df$TB==6 & df$DR==2 & df$HIV==3] <-tb6_mdr_hiv_pos+hiv3
df$cost[df$TB==6 & df$DR==2 & df$HIV==4] <-tb6_mdr_hiv_pos+hiv4_prog1

df$program_id=1

#Make a version of this dataset with costs for community-based care for programs 2 and 3
df_program2<-df
df_program2$cost[df_program2$HIV==4]<-hiv4_prog2
df_program2$cost[df_program2$HIV==4 & (df_program2$TB==2|df_program2$TB==5)]<-hiv4_prog2+on_ipt
df_program2$cost[df_program2$TB==6 & df_program2$DR==1 & df_program2$HIV==4] <-tb6_ds_hiv_pos+hiv4_prog2
df_program2$cost[df_program2$TB==6 & df_program2$DR==2 & df_program2$HIV==4] <-tb6_mdr_hiv_pos+hiv4_prog2


#Match the variables needed to join in the loop in 02_cost eff calculator (variable and program_id)
df_program2$program_id=2
df_program3<-df_program2
df_program3$program_id=3

df_all<-rbind(df, df_program2)
df_all<-rbind(df_all, df_program3)

#Export costs
write.csv(df_all, file=paste0(out_dir,"costs.csv"))


#using tidyverse conditional recoding https://rpubs.com/prlicari13/541675
#Add library(tidyverse) if trying to make this work
#Keeping this block of tidyverse code because it's pretty and may be useful in a limited example.
#The problem here was trying to do sequential blocks it would overwrite the earlier cost variable with NA.
#df<- df %>%
#  mutate(cost = case_when(
#    HIV==2 ~ hiv2,
#    HIV==3 ~ hiv3,
#    HIV==4 ~ hiv4_pol1
#  ))


  