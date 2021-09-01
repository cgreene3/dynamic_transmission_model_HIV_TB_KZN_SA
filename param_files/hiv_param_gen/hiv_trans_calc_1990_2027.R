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
indir <- paste0(here(),'/param_files/hiv_param_gen')
outdir <- paste0(here(),'/param_files')
setwd(indir)

calib_years<-1990:2027
genders <-c('Males', 'Females')
eta_23 <- 0.129533679

hiv_transmission_df<-read_excel('hiv_input_gen_data_1980_2028.xlsx', sheet = 'Sheet1')

hiv_transmission_df$year<-as.integer(hiv_transmission_df$year)

eta_34<-rep(0, times = nrow(hiv_transmission_df))

eta_24<-rep(0, times = nrow(hiv_transmission_df))

for (yr in calib_years){
  for (g in genders){
    
    print(yr)
    
    #before 2004 no ART is available
    if(yr >= 2004){
      
      current_yr_age_gender_eval_temp<-paste0(yr, '_', g)
      nxt_yr_age_gender_eval_temp<-paste0(yr+1, '_', g)
      
      current_yr_row_temp<-which(grepl(current_yr_age_gender_eval_temp, hiv_transmission_df$year_gender))
      nxt_yr_row_temp<-which(grepl(nxt_yr_age_gender_eval_temp, hiv_transmission_df$year_gender))
      
      n4_prop_next_yr <- hiv_transmission_df$art_coverage[nxt_yr_row_temp]
      n4_prop_current_yr <- hiv_transmission_df$art_coverage[current_yr_row_temp]
      n3_prop_current_yr <- hiv_transmission_df$PLHIV_3_CD4_200_less[current_yr_row_temp]
      n2_prop_eligible_current_yr <- hiv_transmission_df$n2_prop_eligible[current_yr_row_temp]
      n2_prop_current_yr <- hiv_transmission_df$PLHIV_2_CD4_200_more[current_yr_row_temp]*
        n2_prop_eligible_current_yr
      
      #to account for changes from mortality and aging out
      hiv_prev_change<-hiv_transmission_df$hiv_prevalence[nxt_yr_row_temp]-
        hiv_transmission_df$hiv_prevalence[current_yr_row_temp]
      
      #relevant_change<-((n3_prop_current_yr+n2_prop_current_yr)*(hiv_prev_change))
      
      
      art_inititation<-(n4_prop_next_yr-n4_prop_current_yr-hiv_prev_change)/(n3_prop_current_yr+n2_prop_current_yr)
      eta_24[current_yr_row_temp]<-art_inititation*n2_prop_eligible_current_yr
      eta_34[current_yr_row_temp]<-art_inititation
      
      #added to v5 to test impacts of increasing ART initiation rates
      if(g == 'Females'){
        if(yr >= 2016){
          eta_24[current_yr_row_temp]<-art_inititation*n2_prop_eligible_current_yr*2
          eta_34[current_yr_row_temp]<-art_inititation*2
        }
      }
      
    }
  }
    
}

hiv_transmission_df$eta_24<-eta_24
hiv_transmission_df$eta_34<-eta_34

#hiv_transmission_df<-hiv_transmission_df%>%
#  filter(year<=2017)


###create plots####
hiv_transition_rate_plot_df<-hiv_transmission_df%>%
  select(c('year', 'gender', 'eta_24', 'eta_34'))%>%
  melt(id.vars = c('year', 'gender'))%>%
  mutate(transition_from = if_else(variable == 'eta_24',
                                   'Transition from HIV+, CD4 >=200', 
                                   'Transition from HIV+, CD4 <200'),
         transition_categories = paste0(gender,": ",transition_from))%>%
  rename(initiation_rate = value)

#png('hiv_transitions_march_1.png', width=450,height=350,res=100)
print(ggplot(hiv_transition_rate_plot_df, 
             aes(x = year, 
                 y = initiation_rate,
                 colour = transition_categories))+
        geom_point() + 
        geom_line(aes(color = transition_categories)))
#   +
#        labs(title = 'test',
#             y = 'Mortality rate')+
#        ylim(0,.5))
#dev.off()

setwd(outdir)
write.csv(hiv_transmission_df, 'hiv_transmission_df_v5.csv')
