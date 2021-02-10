#Estimating background mortality rate from GBD
#Uses (1) ART coverage
#(2) proportion of pop in HIV +, CD4 < 200, & CD4 >= 200,
#(3) prog from CD4 >= 200 to <200 = 0.129533679 per year

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'readxl'), require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/param_files')
setwd(indir)

calib_years<-1990:2017
genders <-c('Males', 'Females')
eta_23 <- 0.129533679

hiv_transmission_df<-read_excel('hiv_input_gen_data.xlsx', sheet = 'Sheet1')

hiv_transmission_df$year<-as.integer(hiv_transmission_df$year)

eta_34<-rep(0, times = nrow(hiv_transmission_df))
eta_24<-rep(0, times = nrow(hiv_transmission_df))

counter = 1

for (r in 1:nrow(hiv_transmission_df)){
  
  #when ART only available to cd4 less
  if (hiv_transmission_df$year[counter] >= 2005){
    if (hiv_transmission_df$year[counter]<2015){
      
      n4_prop_next_yr <- hiv_transmission_df$n4_prop[counter+1]
      n4_prop_current_yr <- hiv_transmission_df$n4_prop[counter]
      n3_prop_current_yr <- hiv_transmission_df$n3_prop[counter]
      
      eta_34[counter]<-(n4_prop_next_yr-n4_prop_current_yr)/n3_prop_current_yr
    }
  }
  
  if (hiv_transmission_df$year[counter]>= 2015){
    n4_prop_next_yr <- hiv_transmission_df$n4_prop[counter+1]
    n4_prop_current_yr <- hiv_transmission_df$n4_prop[counter]
    n3_prop_current_yr <- hiv_transmission_df$n3_prop[counter]
    n2_prop_current_yr <-hiv_transmission_df$n2_prop[counter]
    
    
    art_inititation<-(n4_prop_next_yr-n4_prop_current_yr)/(n3_prop_current_yr+n2_prop_current_yr)
    eta_24[counter]<-art_inititation
    eta_34[counter]<-art_inititation
    
    
    #eta_24[counter]<-eta_34[counter]
    
    #eta_24
    #n2_prop_next_yr<- hiv_transmission_df$n2_prop[counter+1]
    #n2_prop_current_yr<- hiv_transmission_df$n2_prop[counter]
    #n1_prop_current_yr<- hiv_transmission_df$n1_prop[counter]
    #eta_12<-hiv_transmission_df$hiv_incidence[counter]
    
    #numerator<-(n2_prop_next_yr-
    #              n2_prop_current_yr-
    #              (eta_12*n1_prop_current_yr)+
    #              (eta_23*n2_prop_current_yr)
    #)
    
    #eta_24_current_yr<-numerator/(-n2_prop_current_yr)
    #eta_24[counter]<-eta_24_current_yr
    
    #eta_34
   # n4_prop_next_yr<- hiv_transmission_df$n4_prop[counter+1]
  #  n4_prop_current_yr<- hiv_transmission_df$n4_prop[counter]
  #  n3_prop_current_yr<- hiv_transmission_df$n3_prop[counter]
    
    
   # delta_n4<-n4_prop_next_yr-n4_prop_current_yr
  #  n2_move<-(eta_24_current_yr*n2_prop_current_yr)
    
   # eta_34_current_yr<-(delta_n4-n2_move)/(n3_prop_current_yr)
    #eta_34[counter]<-eta_34_current_yr
    
  }
  
  counter <- counter + 1
} 

hiv_transmission_df$eta_24<-eta_24
hiv_transmission_df$eta_34<-eta_34

hiv_transmission_df<-hiv_transmission_df%>%
  filter(year<=2017)

write.csv(hiv_transmission_df, 'hiv_transmission_df.csv')

