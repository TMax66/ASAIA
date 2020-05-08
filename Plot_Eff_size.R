
require(ggplot2)
require(gridExtra)
require(tidyr)
require(ggbeeswarm)

source("geom_flat_violin.R")

#Function that resamples a vector (with replacement) and calculates the median value
boot_mean = function(x) {
  mean(sample(x, replace = TRUE))
}
#Number of bootstrap samples
nsteps=1000
#Confidence level
Confidence_Percentage = 95
Confidence_level = Confidence_Percentage/100
alpha=1-Confidence_level
lower_percentile=(1-Confidence_level)/2
upper_percentile=1-((1-Confidence_level)/2)
i=0


df_raw_data <- read.csv("prova2.csv")#<<-partendo da file in formato long
df_raw_data<-df_raw_data[,-1]
names(df_raw_data)<-c("Condition","value")

prova<-function(data, control){
  df_summary <- data.frame(Condition=levels(factor(data$Condition)), 
                           n=tapply(data$value, data$Condition, length),
                           mean=tapply(data$value,data$Condition, mean))

  df_summary <- data.frame(Condition=levels(factor(data$Condition)), 
                           n=tapply(data$value,data$Condition, length),
                           mean=tapply( data$value,  data$Condition, mean))
  
  df_all_new_mean <- data.frame(Condition=levels(factor( data$Condition)), 
                                resampled_mean=tapply( data$value,
                                                       data$Condition, boot_mean), id=i)
  df_all_new_mean[order(df_all_new_mean$Condition),]
  
  for (i in 1:nsteps) {
    df_boostrapped_mean <- data.frame(Condition=levels(factor(data$Condition)),
                                      resampled_mean=tapply(data$value, 
                                                            data$Condition, boot_mean), id=i)
   df_all_new_mean <- bind_rows(df_all_new_mean, df_boostrapped_mean)
  }
  
  df_summary$ci_mean_hi <- tapply(df_all_new_mean$resampled_mean, 
                                  df_all_new_mean$Condition, quantile, probs=upper_percentile)
  
  df_summary$ci_mean_lo <- tapply(df_all_new_mean$resampled_mean, df_all_new_mean$Condition, 
                                  quantile, probs=lower_percentile)
  
  
  df_all_new_mean[order(df_all_new_mean$Condition),]
  
  df_spread <- spread(df_all_new_mean, key=Condition, value=resampled_mean)
  
  df_spread_differences <- df_spread[,2:ncol(df_spread)] - df_spread[,control]
  
  df_differences <- gather(df_spread_differences, Condition, diff_mean)
   
  df_diff_summary <- data.frame(Condition=levels(factor(df_raw_data$Condition)), 
                                diff_mean=tapply(df_differences$diff_mean, df_differences$Condition, mean))
  df_diff_summary$ci_mean_hi <- tapply(df_differences$diff_mean, df_differences$Condition, quantile, probs=upper_percentile)
  df_diff_summary$ci_mean_lo <- tapply(df_differences$diff_mean, df_differences$Condition, quantile, probs=lower_percentile)
  
  return(list(df_differences, df_diff_summary))
  

  
  
  
}



