####calculate HIV incidence rates####
###determine if there is a clear trend before 2000 and after 2017###


#clean workspace
rm(list = ls())
gc()

#load packages
sapply(c('here', 'dplyr', 'reshape2', 
         'ggplot2', 'readxl', 'RColorBrewer'), 
       require, character.only=T)

#Need to set project (upper R corner of screen) to epi_model_HIV_TB for here to work
indir <- paste0(here(),'/param_files/hiv_incidence_gen')
outdir <- paste0(here(),'/param_files/hiv_incidence_gen/Nov_17')
setwd(indir)

#read data
F15_regions_ref_df<-read.csv('uganda_subcounties.csv')
incidence_rates_df<-read.csv('IHME_SSA_HIV_INC_MORT_2000_2018_INCIDENCE_ADMIN_2_Y2021M01D26.CSV')%>%
  filter(ADM0_NAME == 'Uganda')

#to make sure left and right columns match
F15_regions_ref_df$DName2016<-toupper(F15_regions_ref_df$DName2016)
F15_regions_ref_df$CName2016<-toupper(F15_regions_ref_df$CName2016)
incidence_rates_df$ADM1_NAME<-toupper(incidence_rates_df$ADM1_NAME)

#Some times the incidence rates admin 1 name is matched with D name others its matched with C name

#add D names are in incidence rates
F15_regions_DNAME_df<-F15_regions_ref_df%>%
  left_join(incidence_rates_df, by = c("DName2016" = "ADM1_NAME"))%>%
  mutate(rename(., 'ADM1_NAME' = 'DName2016'))%>%
  select(c('F15Regions', 'ADM1_NAME', 'ADM1_NAME', 'ADM2_NAME', 'year', 'metric', 'mean'))
         
#these are the ones matched with CName (more than F15)
F15_regions_CNAME_df<-F15_regions_ref_df%>%
  filter(CName2016%in% c('BUGWERI',
                         'KAPELEBYONG',
                         'KASSANDA',
                         'KIKUUBE',
                         'KWANIA',
                         'NABILATUK'))%>%
  left_join(incidence_rates_df, by = c("CName2016" = "ADM1_NAME"))%>%
  mutate(rename(., 'ADM1_NAME' = 'CName2016'))%>%
  select(c('F15Regions', 'ADM1_NAME', 'ADM2_NAME', 'year', 'metric', 'mean'))

#remove duplicates

#combine data frames
F15_regions_master_df<-rbind(F15_regions_DNAME_df, F15_regions_CNAME_df)

#remove duplicated rows
F15_regions_master_df<-F15_regions_master_df%>%
  distinct()

#grouped by F15 regions
F15_regions_master_df<-dcast(F15_regions_master_df, 
            F15Regions + ADM1_NAME + ADM2_NAME + year ~ metric,
            fun.aggregate = sum)%>%
  mutate(Population = (Count/Rate)*100000)

F15_regions_summarised_df<-F15_regions_master_df%>%
  group_by(F15Regions, year)%>%
  summarise(total_count = sum(Count),
            total_pop = sum(Population))%>%
  mutate(rate = total_count/total_pop,
         rate_per_100k = rate*100000)

setwd(outdir)
write.csv(F15_regions_summarised_df, 
          'F15_regions_incidence_rates_2000_2018.csv',
          row.names = FALSE)

#create graph
mycolors = c(brewer.pal(name="Blues", n = 8)[6:8],
             brewer.pal(name="Reds", n = 8)[6:8],
             brewer.pal(name="Greens", n = 8)[6:8],
             brewer.pal(name="Purples", n = 8)[6:8],
             brewer.pal(name="Greys", n = 8)[6:8])

#jpeg("F15_regions_incidence_rates_2000_2018.jpeg", dpi=700)

ggsave(filename = "F15_regions_incidence_rates_2000_2018.png", 
       width = 10, height = 10, 
       device='png', dpi=700)

ggplot(F15_regions_summarised_df, 
       aes(x = year, y = rate_per_100k, 
           colour = F15Regions)) +
  geom_line()+
  scale_color_manual(values = mycolors)+
  ylim(0,max(F15_regions_summarised_df$rate_per_100k)+100)

dev.off()

ggsave(filename = "F15_regions_new_infections_2000_2018.png", 
       width = 10, height = 10, 
       device='png', dpi=700)

ggplot(F15_regions_summarised_df, 
       aes(x = year, y = total_count, 
           colour = F15Regions)) +
  geom_line()+
  scale_color_manual(values = mycolors)+
  ylim(0,max(F15_regions_summarised_df$total_count)+100)

dev.off()


uganda_summarised_df<-F15_regions_master_df%>%
  group_by(year)%>%
  summarise(total_count = sum(Count),
            total_pop = sum(Population))%>%
  mutate(rate = total_count/total_pop,
         rate_per_100k = rate*100000)

ggsave(filename = "uganda_new_infections_2000_2018.png", 
       width = 10, height = 10, 
       device='png', dpi=700)

ggplot(uganda_summarised_df, 
       aes(x = year, y = total_count)) +
  geom_line()+
  ylim(0,max(uganda_summarised_df$total_count)+100)

dev.off()
