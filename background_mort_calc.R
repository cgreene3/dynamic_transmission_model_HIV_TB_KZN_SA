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

