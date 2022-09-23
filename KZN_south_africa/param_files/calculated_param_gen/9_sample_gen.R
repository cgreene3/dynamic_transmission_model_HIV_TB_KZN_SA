#clean workspace
rm(list = ls())
gc()


library(readxl)
library(here)
library(dplyr)
library(ExtDist)
library(lhs)
library(ggplot2)
library(stringr)
library(reshape2)

n_samples <- 100000
run_id <- 2 #specify if first, second, third sample, etc.
sample_df <- data.frame(matrix(nrow = n_samples, ncol = 0))

set.seed(as.integer(run_id)) #use a different seed for each sample

#update parameter file
outdir <- paste0(here(),'/param_files/input_parameters')
indir_params <- paste0(here(),'/param_files/calculated_param_gen/input_data')
indir_accepted_params<- paste0(here(), '/calibration_analysis/run_', run_id-1)
outdir_beta_dist_fit<-paste0(here(),'/param_files/distribution_of_accepted_points_graphs')

setwd(indir_params)
model_params_df<-read_excel('KZN_SA_model_parameters.xlsx')%>%
  select(c('model_matched_param', 'min', 'max', 'value', 'calibrated'))

model_params_df$model_matched_param<-str_replace_all(model_params_df$model_matched_param, c("," = "."))

if(run_id>1){
setwd(indir_accepted_params)
best_results_beta_dist<-read.csv('best_results_beta_dist.csv')%>%
  filter(total_in_confidence == 20)
}
n_calib_params<-sum(model_params_df$calibrated == 'yes')
calib_param_names<-model_params_df%>%
  filter(calibrated == 'yes')%>%
  select(c('model_matched_param'))
calib_param_names<-as.vector(t(calib_param_names))

raw_lhs_df <- as.data.frame(randomLHS(n_samples, n_calib_params))
colnames(raw_lhs_df)<-calib_param_names

#find_mode <- function(x) {
#  u <- unique(x)
#  tab <- tabulate(match(x, u))
#  u[tab == max(tab)]
#}

if(run_id == 1){
  
  for (mp in model_params_df$model_matched_param){
    if(model_params_df[model_params_df$model_matched_param == mp, 
                       c('calibrated')] == "yes"){
      
      max_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                c('max')])
      min_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                c('min')])
      
      sample_df_temp<-as.data.frame(qunif(unlist(raw_lhs_df%>%select(mp)),
                      min_temp, max_temp))
      colnames(sample_df_temp)<-mp
      sample_df<-cbind(sample_df, sample_df_temp)
    } else {
      value_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                           c('value')])
      sample_df_temp<-as.data.frame(rep(value_temp, times = n_samples))
      colnames(sample_df_temp)<-mp
      sample_df<-cbind(sample_df, sample_df_temp)
      
    }
  }
}

if(run_id > 1){
  for (mp in model_params_df$model_matched_param){
    if(model_params_df[model_params_df$model_matched_param == mp, 
                     c('calibrated')] == "yes"){
      
      max_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                           c('max')])
      min_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                           c('min')])
      
      #fit to beta distribution
      selected_vals_temp<-unlist(best_results_beta_dist%>%select(mp))
      
      range01 <- function(x){(x-min_temp)/(max_temp-min_temp)} #so can fit to beta dist
      selected_vals_01<-range01(selected_vals_temp)
      
      shapes<-ExtDist::eBeta(selected_vals_01, 
                             method = "MOM")
      
      #hist(selected_vals_01)
      sample_01_beta_temp<-qbeta(unlist(raw_lhs_df%>%select(mp)),
                                        shape1 = shapes$shape1, 
                                        shape2 = shapes$shape2)
      
      sample_df_temp<-as.data.frame(qunif(sample_01_beta_temp,
                                          min_temp, max_temp))
      colnames(sample_df_temp) = 'sample'
      
      bin_width<-(max_temp-min_temp)/10
      bin_width_fitted_values<-(max_temp-min_temp)/10
      bins<-seq(min_temp, max_temp, by = bin_width)
      bins_fitted_values<-seq(min_temp, max_temp, by = bin_width)
      bins_ref_df<-data.frame(bins)
      bins_ref_df<-bins_ref_df%>%
        mutate(bins_array_temp = as.numeric(row.names(bins_ref_df)),
               bins_avg = bins+(bin_width/2))


      selected_vals_bin<-rep(0, times = length(selected_vals_temp))
      fitted_vals_bin<-rep(0, times = nrow(sample_df_temp))

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
        left_join(bins_ref_df, by = c('selected_vals_bin'='bins_array_temp'))%>%
        mutate(variable = 'initial best sample values')
      
      #id bin fitted vals 
      for (i in 1:length(sample_df_temp$sample)){
        for (b in 1:length(bins)){
          if (sample_df_temp$sample[i] >= bins[b]){
            if (sample_df_temp$sample[i] < bins[b+1]){
              fitted_vals_bin[i]<-b
            }
          }
        }
      }
      
      fitted_vals_temp_df<-data.frame(fitted_vals_bin)
      fitted_vals_temp_df<-fitted_vals_temp_df%>%
        mutate(prop_of_total = 1/nrow(fitted_vals_temp_df))%>%
        group_by(fitted_vals_bin)%>%
        summarise(prop_of_total_best = sum(prop_of_total))%>%
        left_join(bins_ref_df, by = c('fitted_vals_bin'='bins_array_temp'))%>%
        mutate(variable = 'revised fitted sample values')
      
      #min_df<-as.data.frame(t(c(0,0,3,min_temp)))
      #max_df<-as.data.frame(t(c(6,0,3,max_temp)))
      #colnames(min_df)<-colnames(fitted_vals_temp_df)
      #colnames(max_df)<-colnames(fitted_vals_temp_df)
      
      #fitted_vals_temp_df<-rbind(fitted_vals_temp_df,min_df,max_df)
      
      fitted_selected_vals_temp_df<-rbind(fitted_vals_temp_df%>%
        select(c('bins_avg', 'variable', 'prop_of_total_best')),
        selected_vals_temp_df%>%select(c('bins_avg', 'variable', 'prop_of_total_best')))
      
      setwd(outdir_beta_dist_fit)
      png(file=paste0(mp, '.jpg'))
      plot(ggplot()+
             geom_bar(data = fitted_selected_vals_temp_df%>%filter(variable != 'revised fitted sample values'),
                      aes(#fill=variable,
                        y=prop_of_total_best,
                        x = bins_avg),
                      position="dodge", stat="identity")+
             ggtitle(mp)+
             xlab('Value')+
             ylab('Density'))
      # plot(ggplot()+
      #   geom_bar(data = selected_vals_temp_df,
      #            aes(#fill=variable,
      #              y=prop_of_total_best,
      #              x = bins_avg),
      #            position="dodge", stat="identity")+
      #   ggtitle(mp)+
      #   xlab('Value')+
      #   ylab('Density')+
      #   geom_line(data = fitted_vals_temp_df, 
      #             aes(y = prop_of_total_best,
      #                 x = bins_avg), colour = 'red'))#+
      dev.off()
      
      
      
    colnames(sample_df_temp)<-mp
    sample_df<-cbind(sample_df, sample_df_temp)
  }else {
    value_temp<-as.numeric(model_params_df[model_params_df$model_matched_param == mp, 
                                           c('value')])
    sample_df_temp<-as.data.frame(rep(value_temp, times = n_samples))
    colnames(sample_df_temp)<-mp
    sample_df<-cbind(sample_df, sample_df_temp)
  }
  }
}

sample_df$sim_id<-seq(1:nrow(sample_df))

setwd(outdir)
write.csv(sample_df, paste0('sample_df_', run_id, '.csv'), row.names = FALSE)



#######





# png(file=paste0(mp, '.jpg'))
# t_hist<-hist(selected_vals_temp, prob = TRUE, main = mp)
# t_hist$density<-t_hist$density/sum(t_hist$density)
# plot(t_hist, freq = FALSE, main = mp)
# t_density<-density(sample_df_temp$sample)
# if(max(t_density$y) > 1){
#   t_density$y <-(t_density$y/100)
# }
# lines(t_density, col = 'red')
# dev.off()
#overlay sampled points
#density function
#prop in each bucket.
# 
# ggplot(data = hist_df, aes(x = selected_vals_temp))+
#      geom_histogram(aes(y = stat(count/sum(count))),
#                     binwidth = bin_width)

#geom_bar(data = vals_temp_df_melt%>%filter(variable=='prop_of_total_best'),
#         aes(#fill=variable,
#             y=value,
#             x = bins_avg),
#         position="dodge", stat="identity")+
#geom_line(data = vals_temp_df_melt%>%filter(variable=='prop_of_total_fitted'),
#         aes(#fill=variable,
#           y=value,
#           x = bins_avg),stat="identity")+
#ggtitle(mp) +
#xlab("Value") + ylab("Density")
# #x<-rbeta(1000, shape1 = shapes$shape1, shape2 = shapes$shape2)
# #density_df<-as.data.frame(qunif(x, min_temp, max_temp))
# #colnames(density_df) = 'prob'
#
# # selected_vals_df<-as.data.frame(selected_vals_temp)
# # 
# # 
# # 
# # ggplot()+
# #   geom_histogram(data = selected_vals_df,
# #                  mapping = aes(x = selected_vals_temp,
# #                                y = stat(count/sum(count))),
# #                  binwidth = bin_width)+
# #   geom_density(data = sample_df_temp, 
# #                mapping = aes(x = sample))
# # 
# # ggplot()+
# #   #geom_density(sample_df_temp, mapping=aes(x = sample))+
# #   geom_histogram(best_results_beta_dist%>%select(mp), 
# #                  mapping = aes(x = mp,
# #                                y = stat(count/sum(count))))
# #   ggplot(data = best_results_beta_dist, mapping=aes(x = mp))+
# #   geom_histogram(aes(y = stat(count/sum(count))),
# #                                   binwidth = bin_width)
# #   #geom_histogram(data = selected_vals_temp_df)
# #   geom_density(, aes(sample))
# # 

#   #+
#geom_density(data = density_df, aes(prob))
# 
# 
# 
# 
# 
# 
# 
# hist_df<-as.data.frame(cbind(selected_vals_temp,
#                bins_array_temp))
# hist_df<-hist_df%>%
#   mutate(prop_of_total = 1/nrow(hist_df))%>%
#   group_by(bins_array_temp)%>%
#   summarise(prop_of_total_ = sum(prop_of_total))
#   
#   count(bins_array_temp)%>%
#   ungrou
#   mutate(prop_of_total = n/sum(n))
#   group_by(bins_array_temp)%>%
#   summarise(prop_of_total_selected_vals = nrow(bins_array_temp))
# 
# 
# 
# t<-ggplot(data = hist_df, aes(x = selected_vals_temp))+
#   geom_histogram(aes(y = stat(count/sum(count))),
#                  binwidth = bin_width)
# 
# 
# 
#   geom_line()
# 
# 
# +
#   stat_function(fun = qbeta, args = list(shape1 = shapes$shape1, 
#                                          shape2 = shapes$shape2))
# 
# t<-t+lines(density(sample_df_temp$sample), 
#            col = "red") 
# 
# ggplot(sample_df_temp, aes(x=weight)) + 
#   geom_density()
# 
# hist_df<-hist_df%>%
#   
# 
# ggplot(data = hist_df, aes(bins_array_temp))+
#   geom_bar()
# 
# hist_df<-hist_df%>%
#   mutate(Density = 1/nrow(hist_df))%>%
#   group_by(bins_array_temp)%>%
#   summarise(Density = sum(Density))
# #   left_join()
# 
# setwd(outdir_beta_dist_fit)
# 
# hist_df<-as.data.frame(selected_vals_temp)
# 
# 
# 
# 
# ggplot(data = )
#   geom_density(data = sample_df_temp$`qunif(sample_01_beta_temp, min_temp, max_temp)`)
#                
# 
# 
# png(file=paste0(mp, '.jpg'))
# hist_info<-hist(selected_vals_temp, 
#                 plot = FALSE)
# hist_info$density <- hist_info$counts /    # Compute density values
#   sum(hist_info$counts)*100
# plot(hist_info, freq = FALSE)              # Plot histogram with percentages
# lines(density(sample_df_temp$sample), 
#       col = "red")   
# dev.off()