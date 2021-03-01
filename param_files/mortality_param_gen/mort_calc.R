#Estimating background mortality rate from GBD
#Uses GBD 2019 - Number of Deaths
#SA population estimates

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2', 'stringr'), require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/mortality_param_gen')
outdir <- paste0(here(),'/param_files')
setwd(indir)



#read in pop files

#put GBD estimates into seperate folder so can read in all
setwd(paste0(indir,'/GBD_pop_estimates')) 
pop_files = list.files(pattern="*.CSV")

pop_df<-data.frame()

lapply(pop_files, function(file){
  temp_df <-read.csv(file)%>%
    filter(location_id == 196,
           age_group_id >= 8, #between 15
           age_group_id <=15, #less than 60
           sex_id != 3) #only females and males
  pop_df <<- rbind(pop_df, temp_df)
})

pop_df<-pop_df%>%
  group_by(sex_id, sex_name, year_id)%>%
  summarise(expected_total_pop = sum(val))%>%
  rename(year = year_id)

#prevalence and mortality data frames consider the following disease states
#HIV/AIDS - causeid 298
#DS TB - cause id 934
#MDR TB - cause ids 946, 947

#read in prevalence data
setwd(indir)
prev_df<-read.csv('IHME_prev_1990_2017_Feb28.csv')
prev_df<-prev_df%>%
  group_by(sex_id, year)%>%
  summarise(expected_disease_pop = sum(val))
  

non_disease_mort_df<-pop_df%>%
  left_join(prev_df, by = c('sex_id', 'year'))%>%
  mutate(expected_non_disease_pop = expected_total_pop-expected_disease_pop)

#read in mortality data
mort_df<-read.csv('IHME_deaths_1990_2017_Feb28.csv')%>%
  mutate(disease_cat = if_else(cause_name == 'All causes', 'base', 'disease'))%>%
  group_by(sex_id, year, disease_cat)%>%
  summarise(expected_disease_mort = sum(val))%>%
  ungroup()%>%
  mutate(base_mort_calcs = if_else(disease_cat == 'base', 
                                   expected_disease_mort, 
                                   -expected_disease_mort))%>%
  group_by(sex_id, year)%>%
  summarise(expected_non_disease_mort_total = sum(base_mort_calcs))

non_disease_mort_df<-non_disease_mort_df%>%
  left_join(mort_df, by = c('sex_id', 'year'))%>%
  mutate(expected_non_disease_mort_rate = expected_non_disease_mort_total/expected_non_disease_pop)



non_disease_mort_graph<-ggplot(non_disease_mort_df, aes(x = year, 
                                y = expected_non_disease_mort_rate,
                                group = sex_name)) + 
  geom_point() + 
  geom_line(aes(color = sex_name))+
  labs(title = "Non-disease specific mortality\novertime by gender",
       y = 'Non-disease specific mortality rate',
       colour = "Gender")+
  ylim(0,.01)

png('non_disease_mort_graph.png', width=450,height=350,res=100)
print(non_disease_mort_graph)
dev.off()


#write data to read into epi model

#calculate base death rates by gender
df_input_data <- non_disease_mort_df

TB_SET <- 1:8
HIV_SET <- 1:4
#calculate disease specific mort
df_input_data<-do.call("rbind", replicate(length(TB_SET)*length(HIV_SET), 
                                          df_input_data, 
                                          simplify = FALSE))

df_input_data<-df_input_data%>%
  mutate(TB_compartment = 0,
         HIV_compartment = 0)%>%
  rename(G_compartment = sex_id)

df_input_data$year<-as.integer(df_input_data$year)
df_input_data<-df_input_data%>%
    arrange(year)%>%
    arrange(G_compartment)

counter <-1

for (yr in unique(df_input_data$year)){
    for (g in 1:2){
        for (t in TB_SET){
            for (h in HIV_SET){
                df_input_data$TB_compartment[counter]<-t
                df_input_data$HIV_compartment[counter]<-h
                counter<-counter+1
            }
        }
    }
}

#adjusted mortality
HIV_2_increase <- 5
HIV_3_increase <- 10
HIV_4_increase <- 1.2

TB_Active_HIV_1<-20 #20x total increase in mort rate from non-disease mort
TB_Active_HIV_2<-5 #50x total increase in mort rate from non-disease mort
TB_Active_HIV_3<-10 #100x total increase in mort rate from non-disease mort
TB_Active_HIV_4<-25 #30x total increase in mort rate from non-disease mort

df_input_data<-df_input_data%>%
    mutate(hiv_adj = if_else(HIV_compartment == 1, 1,
                           if_else(HIV_compartment == 2, HIV_2_increase,
                                   if_else(HIV_compartment == 3, HIV_3_increase,
                                           HIV_4_increase))))%>%
    mutate(tb_adj = if_else(TB_compartment!=6, 1,
                            if_else(HIV_compartment == 1, TB_Active_HIV_1,
                                    if_else(HIV_compartment == 2, TB_Active_HIV_2,
                                            if_else(HIV_compartment == 3, TB_Active_HIV_3,
                                                    TB_Active_HIV_4)))))%>%
    mutate(mort_rate = expected_non_disease_mort_rate*hiv_adj*tb_adj)

mort_rate_plots_df<-df_input_data%>%
  mutate(cause = if_else((HIV_compartment==1 & TB_compartment!=6), 
                         'TB/HIV Negative',
                         if_else((HIV_compartment==1 & TB_compartment==6),
                                 'Active TB, HIV Negative',
                         if_else((HIV_compartment!=1 & TB_compartment!=6),
                         'No TB, HIV Positive',
                         'Active TB, HIV Positive'))))%>%
  group_by(year, cause, sex_name)%>%
  summarise(mort_rate = mean(mort_rate))%>%
  rename(Gender = sex_name)


for (cause1 in unique(mort_rate_plots_df$cause)){
  plot_file_name<-str_replace_all(string=cause1, pattern=" ", repl="_")
  plot_file_name<-str_replace_all(string=plot_file_name, pattern=",", repl="")
  plot_file_name<-str_replace_all(string=plot_file_name, pattern="/", repl="_")
  png(paste0(plot_file_name,'.png'), width=450,height=350,res=100)
  print(ggplot(mort_rate_plots_df%>%filter(cause == cause1), aes(x = year, 
                                 y = mort_rate,
                                 group = Gender)) + 
    geom_point() + 
    geom_line(aes(color = Gender))+
    labs(title = paste0("Mortality overtime for ", cause1),
         y = 'Mortality rate')+
    ylim(0,.5))
  dev.off()
}



setwd(outdir)
write.csv(df_input_data, 'mort_df.csv')

