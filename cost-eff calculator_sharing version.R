#Skeleton code of cost calculation for community IPT analysis
#July 2021

rm(list=ls())


library(dplyr)
library(reshape2)

home_dir <- ("Y:/Research Projects/Ross K01 transmission model/")
in_dir <- (paste0(home_dir, "Data/"))
out_dir <- (paste0(home_dir, "Analysis/"))

#Import model results
#--Input file will be population in each compartment 2018-2027 under each policy (x20? rows for different parameter combinations)
outputs <- read.csv(file=paste0(home_dir, "output/20210721/past_2018/out_df_sim_id_5656_to_2028.csv"))

#Temporary fix to add policy to this data frame. Anticipate deleting these next three lines when results have code in them
outputs$policy <- 1
outputs2<- read.csv(file=paste0(home_dir, "output/20210721/past_2018/out_df_sim_id_5656_to_2028.csv"))
outputs2$policy <-2
outputs<-rbind(outputs, outputs2)

##Separating analysis datasets--------------------------------------------------------------------------------------
#Select one row per year (midpoint ends in .5) for costing analysis years 2018-2027
mids<-c(28.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5)
output_mids<-filter(outputs, time %in% mids)

#Split off mortality and incidence data since those are cumulative and analysis only uses the last row
death_vars<-c("year", "sim_id", "time", "policy", "TB_mort_HIV_neg_male", "TB_mort_HIV_neg_female", "TB_mort_HIV_pos_male", "TB_mort_HIV_pos_female")
deaths<-outputs[death_vars]
deaths_end<-deaths[ which(deaths$time==38),] #Selects the last row for the final cumulative values
inc_vars<-c("year", "sim_id", "time", "policy", "TB_incidence_HIV_neg_male", "TB_incidence_HIV_neg_female", "TB_incidence_HIV_pos_male", "TB_incidence_HIV_pos_female")
incidence<-outputs[inc_vars]
incidence_end<-incidence[ which(incidence$time==38),] #Selects the last row for the final cumulative values

##Total deaths averted and cases averted per policy
#Holding space to do this calculation before reshaping the datasets

##YLLs calculations--------------------------------------------------------------------------------------------------
#Import mortality age fractions for YLLs
tb_male_mort_fractions<-read.csv(file=paste0(out_dir, "Mort_age_fractions_TB_male.csv"))
tb_female_mort_fractions<-read.csv(file=paste0(out_dir, "Mort_age_fractions_TB_female.csv"))
tbhiv_male_mort_fractions<-read.csv(file=paste0(out_dir, "Mort_age_fractions_TBHIV_male.csv"))
tbhiv_female_mort_fractions<-read.csv(file=paste0(out_dir, "Mort_age_fractions_TBHIV_female.csv"))

#Simplify the age fraction dataframes to vectors to use them more easily in next calcs.
tb_male_mort_vec<-tb_male_mort_fractions$mean
tb_female_mort_vec<-tb_female_mort_fractions$mean
tbhiv_male_mort_vec<-tbhiv_male_mort_fractions$mean
tb_hiv_female_mort_vec<-tbhiv_female_mort_fractions$mean

#Reshape deaths dataset to multiply by mort
deaths_end_long<-melt(deaths_end,  id.vars=c("year", "sim_id", "time", "policy")) #Will need to think more about sim_id

gbd_age_groups<-c("age_8", "age_9","age_10", "age_11", "age_12","age_13", "age_14", "age_15","age_16")
deaths_end_long[gbd_age_groups]<-NA

#Multiply each mort fraction by the endpoint deaths - Stuck here now

deaths_end_long <- deaths_end_long[deaths_end_long$variable=='TB_mort_HIV_neg_male',]
#multiply deaths_end_long$value by the vector tb_male_mort_vec to get the values for columns 7:15


#This way worked with a single policy row, but doesn't hold up to multiple policies or sim_ids
tb_num_mort_male<-deaths_end$TB_mort_HIV_neg_male*tb_male_mort_fractions$mean
tb_num_mort_female<-deaths_end$TB_mort_HIV_neg_female*tb_female_mort_fractions$mean
tbhiv_num_mort_male<-deaths_end$TB_mort_HIV_pos_male*tbhiv_male_mort_fractions$mean
tbhiv_num_mort_female<-deaths_end$TB_mort_HIV_pos_female*tbhiv_female_mort_fractions$mean

df_male<-rbind(tb_num_mort_male, tbhiv_num_mort_male) #Next could mult this by life table values for YLL male
df_female<-rbind(tb_num_mort_female, tbhiv_num_mort_female)



#Import the cost worksheet
costs<-read.csv(paste0(out_dir, "ART_IPT_model_cost_inputs.csv"))
#--Data management to sum the costs that apply to each compartment and policy



#Multiply the per-person cost of each compartment-year-policy combination by the population of each compartment for each parameter set

#Sum the costs over years and compartments to generate a set of costs per policy (same number of estimates as parameter sets)

#Calculate the mean and CIs of the cost estimates per policy

#Import the mean number of TB cases and deaths averted (or calculate them here?) for the top parameter sets

#Divide the mean cost estimate by number of TB cases averted to calculate cost per case averted

#Divide the mean cost estimate by number of deaths averted to calculate cost per death averted



#Calculate DALYs averted using steps outlined by Marx et al 
#--Apply YLL and YLD to TB and HIV deaths in each policy
#--Sum YLL and YLD for DALYS from deaths
#--Apply YLD to HIV and TB compartments in each policy
#--Sum with DALYs from deaths to get full DALYS
#--Divide total cost by DALYs to get costs per DALY averted


#Estimate YLL using distribution of ages at death from HIV and TB among adults in South Africa, by year
#Calculate the average age of death (by year, gender, TB, TB-HIV) in South Africa OR calculate the proportion of deaths by age group in a year?

#Extract the life expectancy by year - how to figure this out across age groups?
#http://ghdx.healthdata.org/record/ihme-data/gbd-2019-life-tables-1950-2019