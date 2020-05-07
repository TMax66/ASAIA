
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

prova<-function(data){
  df_summary <- data.frame(Condition=levels(factor(data$Condition)), 
                           n=tapply(data$value, data$Condition, length),
                           mean=tapply(data$value,data$Condition, mean))
  return(df_summary)

}


df_summary <- data.frame(Condition=levels(factor(df_raw_data$Condition)), 
                         n=tapply(df_raw_data$value, df_raw_data$Condition, length),
                         mean=tapply(df_raw_data$value, df_raw_data$Condition, mean))
#Generate a dataframe that collects all bootstrapped median values
df_all_new_mean <- data.frame(Condition=levels(factor(df_raw_data$Condition)), 
                              resampled_mean=tapply(df_raw_data$value,
                                                    df_raw_data$Condition, boot_mean), id=i)