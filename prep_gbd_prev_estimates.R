#This file takes gbd prevalence estimates from a csv generated using the results tool 
#http://ghdx.healthdata.org/gbd-results-tool and aggregates them to use in model calibration.
#Data are already
#July 28, 2020

rm(list = ls())

#Define input directory
#indir<-"~/ws/GitHub/ws.epi_model_HIV_TB/epi_model_HIV_TB/param_files/" #chelsea
indir<-("C:/Users/jross/repos/epi_model_HIV_TB/param_files/")

#Define output directory - Set to outdir loc
outdir<-("~/GitHub/epi_model_HIV_TB/test_outputs/Jun25") #chelsea
#outdir<- #jen

df<-read.csv(paste0(indir, "GBD_prev_1990_2017_jul28.csv"))

#Section 1 - Assembling calibration dataset. For each year (1990 to 2017) generate TB prevalence rate (per 100,000 population) by HIV+, HIV-, and male/female

#Group 1 - Active TB among HIV-negative males. One value per year (1990 to 2017). Use the number in "val"
    #For each year, sum across these three estimates: sex_id==1 & (cause_id==934|cause_id==946|cause_id==947)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 2  - Active TB prevalence among HIV-negative females
    #For each year, sum across these three estimates: sex_id==2 & (cause_id==934|cause_id==946|cause_id==947)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 3 - Active TB HIV+ males
    #For each year, sum across these three estimates: sex_id==1 & (cause_id==948|cause_id==949|cause_id==950)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 4 - Active TB HIV+ females
    #For each year, sum across these three estimates: sex_id==2 & (cause_id==948|cause_id==949|cause_id==950)
    #Use metric_id==3 to get rates and measure_id==5 for prevalence



#Section 2 - Line graph of values in groups 1-4 with years on x axis



#Section 3 (later priority) - LTBI prevalence estimates over time

#Group 5 - LTBI prevalence over time for males
    #Use sex_id==1 and cause_id==954
    #Use metric_id==3 to get rates and measure_id==5 for prevalence

#Group 6 - LTBI prevalence over time for females
    #Use sex_id==2 and cause_id==954
    #Use metric_id==3 to get rates and measure_id==5 for prevalence


 
#Section 4 - Line graph of values in groups 5&6 with years on x-axis. Values should be substantially higher than groups 1-4.

