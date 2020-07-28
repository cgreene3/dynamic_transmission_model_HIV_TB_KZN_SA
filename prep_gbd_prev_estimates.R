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

#For each year (1990 to 2017) generate TB prevalence rate (per 100,000 population) by HIV+, HIV-, and male/female

#Skeleton code that needs to be converted to apply functions. More to come here.
#Active TB HIV negative male

#Active TB prevalence HIV negative female

#Active TB HIV+ male

#Active TB HIV+ female
