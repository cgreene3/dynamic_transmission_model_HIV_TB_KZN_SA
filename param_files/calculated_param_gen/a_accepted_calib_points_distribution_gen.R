#clean workspace
rm(list = ls())
gc()

#generates histograms of accepted parameter values
#used to assess calibrated parameter value ranges fit

library(readxl)
library(here)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)

#update parameter file
outdir_sample <- paste0(here(),'/param_files/input_parameters')
indir_params <- paste0(here(),'/param_files/calculated_param_gen/input_data')
outdir_graphs <-paste0(here(),'/param_files/distribution_of_accepted_points_graphs')
indir_calib_analysis<-paste0(here(), '/calibration_analysis/')

setwd(indir_calib_analysis)
best_results_df_plot_dist<-read.csv("almost_best_calibration_sets_ref_df.csv")#%>%
  #filter(total_in_confidence == 20) #change if want to see dist of parameter sets
#with 18, 19, or 20 in calib targets

setwd(indir_params)
model_params_df<-read_excel('KZN_SA_model_parameters.xlsx')%>%
  select(c('model_matched_param', 'min', 'max', 'value', 'calibrated'))

model_params_df$model_matched_param<-str_replace_all(model_params_df$model_matched_param, c("," = "."))

n_calib_params<-sum(model_params_df$calibrated == 'yes')
calib_param_names<-model_params_df%>%
  filter(calibrated == 'yes')%>%
  select(c('model_matched_param'))
calib_param_names<-as.vector(t(calib_param_names))

#raw_lhs_df <- as.data.frame(randomLHS(n_samples, n_calib_params))
#colnames(raw_lhs_df)<-calib_param_names

for (mp in model_params_df$model_matched_param){
  if(model_params_df[model_params_df$model_matched_param == mp, 
                     c('calibrated')] == "yes"){
    
    max_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                         c('max')])
    min_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                         c('min')])
    
    #fit to beta distribution
    selected_vals_temp<-unlist(best_results_df_plot_dist%>%select(mp))
    
    bin_width<-(max_temp-min_temp)/10
    bins<-seq(min_temp, max_temp, by = bin_width)
    bins_ref_df<-data.frame(bins)
    bins_ref_df<-bins_ref_df%>%
      mutate(bins_array_temp = as.numeric(row.names(bins_ref_df)),
             bins_avg = bins+(bin_width/2))
    
    selected_vals_bin<-rep(0, times = length(selected_vals_temp))
    
    #id bin selected vals 
       for (i in 1:length(selected_vals_temp)){
         for (b in 1:length(bins)){
           if (selected_vals_temp[i] >= bins[b]){
             if (selected_vals_temp[i] < bins[b+1]){
               selected_vals_bin[i]<-b
             }
           }
         }
       }
    
    
    selected_vals_temp_df<-data.frame(selected_vals_bin)
    selected_vals_temp_df<-selected_vals_temp_df%>%
      mutate(prop_of_total = 1/nrow(selected_vals_temp_df))%>%
      group_by(selected_vals_bin)%>%
      summarise(prop_of_total_best = sum(prop_of_total))%>%
      left_join(bins_ref_df, by = c('selected_vals_bin'='bins_array_temp'))
    
    setwd(outdir_graphs)
    png(file=paste0(mp, '.jpg'))
    plot(ggplot()+
           geom_bar(data = selected_vals_temp_df,
                    aes(y = prop_of_total_best, x = bins_avg),
                    stat = "identity")+
           xlim(min_temp, max_temp)+
           ggtitle(mp)+
           xlab('Value')+
           ylab('Density'))
    dev.off()
    
    
    #hist(selected_vals_temp,
    #     main = mp,
    #     #xlim = c(min_temp, max_temp),
    #     freq = FALSE,
    #     xlab = "value",
    #     breaks = 10)
    #dev.off()
  #   
  #   #range01 <- function(x){(x-min_temp)/(max_temp-min_temp)} #so can fit to beta dist
  #   #selected_vals_01<-range01(selected_vals_temp)
  #   
  #   #shapes<-ExtDist::eBeta(selected_vals_01, 
  #   #                       method = "MOM")
  #   
  # 
  #   #sample_01_beta_temp<-qbeta(unlist(raw_lhs_df%>%select(mp)),
  #   #                           shape1 = shapes$shape1, 
  #   #                           shape2 = shapes$shape2)
  #   
  #   #sample_df_temp<-as.data.frame(qunif(sample_01_beta_temp,
  #   #                                    min_temp, max_temp))
  #   #colnames(sample_df_temp) = 'sample'
  #   
  #   bin_width<-(max_temp-min_temp)/10
  #   bins<-seq(min_temp, max_temp, by = bin_width)
  #   bins_ref_df<-data.frame(bins)
  #   bins_ref_df<-bins_ref_df%>%
  #     mutate(bins_array_temp = as.numeric(row.names(bins_ref_df)),
  #            bins_avg = bins+(bin_width/2))
  #   
  #   
  #   selected_vals_bin<-rep(0, times = length(selected_vals_temp))
  #   #fitted_vals_bin<-rep(0, times = nrow(sample_df_temp))
  #   
  #   #id bin selected vals 
  #   for (i in 1:length(selected_vals_temp)){
  #     for (b in 1:length(bins)){
  #       if (selected_vals_temp[i] >= bins[b]){
  #         if (selected_vals_temp[i] < bins[b+1]){
  #           selected_vals_bin[i]<-b
  #         }
  #       }
  #     }
  #   }
  #   
  #   selected_vals_temp_df<-data.frame(selected_vals_bin)
  #   selected_vals_temp_df<-selected_vals_temp_df%>%
  #     mutate(prop_of_total = 1/nrow(selected_vals_temp_df))%>%
  #     group_by(selected_vals_bin)%>%
  #     summarise(prop_of_total_best = sum(prop_of_total))%>%
  #     left_join(bins_ref_df, by = c('selected_vals_bin'='bins_array_temp'))%>%
  #     mutate(variable = 'initial best sample values')
  #   
  #   #id bin fitted vals 
  #   # for (i in 1:length(sample_df_temp$sample)){
  #   #   for (b in 1:length(bins)){
  #   #     if (sample_df_temp$sample[i] >= bins[b]){
  #   #       if (sample_df_temp$sample[i] < bins[b+1]){
  #   #         fitted_vals_bin[i]<-b
  #   #       }
  #   #     }
  #   #   }
  #   # }
  #   
  #   #fitted_vals_temp_df<-data.frame(fitted_vals_bin)
  #   #fitted_vals_temp_df<-fitted_vals_temp_df%>%
  #   #  mutate(prop_of_total = 1/nrow(fitted_vals_temp_df))%>%
  #   #  group_by(fitted_vals_bin)%>%
  #   #  summarise(prop_of_total_best = sum(prop_of_total))%>%
  #   #  left_join(bins_ref_df, by = c('fitted_vals_bin'='bins_array_temp'))%>%
  #   #  mutate(variable = 'revised fitted sample values')
  #   
  #   #fitted_selected_vals_temp_df<-rbind(fitted_vals_temp_df%>%
  #   #                                      select(c('bins_avg', 'variable', 'prop_of_total_best')),
  #   #                                    selected_vals_temp_df%>%select(c('bins_avg', 'variable', 'prop_of_total_best')))
  #   
  #   selected_vals_temp_df<-selected_vals_temp_df%>%
  #     select(c('bins_avg', 'prop_of_total_best'))
  #   
  #   setwd(outdir_graphs)
  #   png(file=paste0(mp, '.jpg'))
  #   plot(ggplot()+
  #          geom_bar(data = selected_vals_temp_df,
  #                   aes(y = prop_of_total_best, x = bins_avg),
  #                   stat = "identity")+
  #          xlim(min_temp, max_temp)
  #   plot(ggplot(data=selected_vals_temp_df, 
  #               aes(x=bins_avg, y=prop_of_total_best)) +
  #          geom_bar(stat="identity")
  #        
  #        
  #     ggplot()+
  #          geom_bar(data = selected_vals_temp_df%>%
  #                     select(c('bins_avg', 'prop_of_total_best')),
  #                   aes(y=prop_of_total_best,
  #                     x = bins_avg),
  #                   #position="dodge", 
  #                   stat="identity")+
  #          ggtitle(mp)+
  #          xlab('Value')+
  #          ylab('Density')
  #     )
  #   dev.off()
  #   
  #   
  #   
  #   #colnames(sample_df_temp)<-mp
  #   #sample_df<-cbind(sample_df, sample_df_temp)
  # #}else {
  # #  value_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
  # #                                         c('value')])
  # #  sample_df_temp<-as.data.frame(rep(value_temp, times = n_samples))
  # #  colnames(sample_df_temp)<-mp
  # #  sample_df<-cbind(sample_df, sample_df_temp)
  }
}


##if generating another sample
#setwd(outdir_sample)
#write.csv(sample_df, 'calibration_sets_df_beta.csv', row.names = FALSE)