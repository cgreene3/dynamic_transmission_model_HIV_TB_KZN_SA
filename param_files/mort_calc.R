#Estimating background mortality rate from GBD
#Uses GBD 2019

#Still need to check denominators

#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 'ggplot2'), require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files')
outdir <- paste0(here(),'/param_files')
setwd(indir)

all_cause <- read.csv('gbd_2019_all-cause_mort_rates_south_africa_1990-2019_20210126.csv')
tb <- read.csv('gbd_2019_tb_mort_rates_south_africa_1990-2019_20210126.csv')
hiv <- read.csv('gbd_2019_hiv_mort_rates_south_africa_1990-2019_20210126.csv')

#prep data for merge with dplyr renaming
tb <- tb %>% rename(tb_val=val)
tb <- tb %>% select(measure_id, measure_name, location_name, sex_name, age_name, metric_id, metric_name, year, tb_val)

hiv <- hiv %>% rename(hiv_val=val)
hiv <- hiv %>% select(measure_id, measure_name, location_name, sex_name, age_name, metric_id, metric_name, year, hiv_val)

all_cause <- all_cause %>% rename(all_val=val)
all_cause <- all_cause %>% select(measure_id, measure_name, location_name, sex_name, age_name, metric_id, metric_name, year, all_val)

#merge with dplyr syntax
df1 <- left_join(all_cause, tb)
df <- left_join(df1, hiv)
rm(df1)
df$base_val <- df$all_val-df$tb_val-df$hiv_val

#reshape data wide to long for plotting using reshape2 package
df_long <- melt(df, id.vars=c('measure_id', 'measure_name', 'location_name', 'sex_name', 'age_name', 'metric_id', 'metric_name', 'year'),
                measure.vars = c('all_val', 'tb_val', 'hiv_val', 'base_val'), variable.name = 'cause', value.name = 'rate')


#plot
p1 <- ggplot(data=df_long, aes(x=year, y=rate, group=cause))+
    geom_line(aes(color=cause))+
    facet_wrap(~sex_name,  ncol=3)+
    labs(title="All-cause and cause-specific mortality rates among adults 15-49 in South Africa",x="Year", y = "Mortality rate per 100,000 population")+
    theme_minimal()
    
p1

#write data to read into epi model

#calculate base death rates by gender
df_input_data <- df_long%>%
    filter(sex_name != 'Both',
           cause == 'base_val')%>%
    mutate(G_compartment = if_else(sex_name == 'Male', 1, 2),
           base_value = rate/100000)%>%
    select(c('year', 'base_value', 'G_compartment'))

TB_SET <- 1:8
HIV_SET <- 1:4
#calculate disease specific mort
df_input_data<-do.call("rbind", replicate(length(TB_SET)*length(HIV_SET), df_input_data, simplify = FALSE))
df_input_data$TB_compartment <-rep(0, times = length(df_input_data))
df_input_data$HIV_compartment <- rep(0, times = length(df_input_data))
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

TB_Active_HIV_1<-20
TB_Active_HIV_2<-50
TB_Active_HIV_3<-100
TB_Active_HIV_4<-30

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
    mutate(total_mort = base_value*hiv_adj*tb_adj)%>%
    select(c('year', 'TB_compartment',  
             'HIV_compartment', 'G_compartment',
             'total_mort'))

setwd(outdir)
write.csv(df_input_data, 'mort_df.csv')
