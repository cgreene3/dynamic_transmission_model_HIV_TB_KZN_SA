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
#Select one row per year (midpoint ends in .5) for costing analysis years 2018-2027 #Consider improving with modulo %% or seq
mids<-seq(from=28.5, to = 37.5, by=1)
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

#Reshape deaths dataset to multiply by mort
deaths_end_long<-melt(deaths_end,  id.vars=c("year", "sim_id", "time", "policy")) #Will need to think more about sim_id

#Add the matching variables to the mortality fraction dataframes, join them to the age distributions, and multiply by age fractions
tb_male_mort_fractions$variable<-"TB_mort_HIV_neg_male"
tb_female_mort_fractions$variable<-"TB_mort_HIV_neg_female"
tbhiv_male_mort_fractions$variable<-"TB_mort_HIV_pos_male"
tbhiv_female_mort_fractions$variable<-"TB_mort_HIV_pos_female"
df<-rbind(tb_male_mort_fractions, tb_female_mort_fractions, tbhiv_male_mort_fractions, tbhiv_female_mort_fractions)
df2<-df %>% left_join(deaths_end_long, by=c("variable"))
df2$fractional_deaths<-df2$mean*df2$value

#Next step to bring in life tables to finish YLLs calculation



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