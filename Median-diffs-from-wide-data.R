# Script that:
# Reads csv file (with header)
# Calculates median and 95%CI by a percentile bootstrap
# Calculates the difference between medians and the 95%CI, reflecting the magnitude and uncertainty of the difference
# Plots: (A) the data, (B) bootstrapped medians, (C) differences between bootstrapped medians and (D) the median and 95% CI
# Created by: 
# Joachim Goedhart, @joachimgoedhart, 2018


#Requires the packages ggplot2, gridExtra, tidyr
require(ggplot2)
require(gridExtra)
require(tidyr)
require(ggbeeswarm)

#Requires a new geom for ggplot2 that plots half of a violin plot (for the distribution of bootstrapped differences)

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

#Read the (non-tidy) data in csv format
df_raw_wide <- read.csv("prova.csv")

df_raw_wide<-df_raw_wide[,-1]
#Convert the wide data to a long format (and remove NA)
df_raw_data <-  na.omit(gather(df_raw_wide, Condition, value))


#Generate a dataframe that keeps a summary of the data
df_summary <- data.frame(Condition=levels(factor(df_raw_data$Condition)), 
          n=tapply(df_raw_data$value, df_raw_data$Condition, length),
          mean=tapply(df_raw_data$value, df_raw_data$Condition, mean))

#Generate a dataframe that collects all bootstrapped median values
df_all_new_mean <- data.frame(Condition=levels(factor(df_raw_data$Condition)), 
                      resampled_mean=tapply(df_raw_data$value,
                      df_raw_data$Condition, boot_mean), id=i)


#####################################################################################
######################## Calculate Median and 95%CI by bootstrap ###################

#Perform the resampling nsteps number of times (typically 1,000-10,000x)
for (i in 1:nsteps) {

	#Caclulate the median from a boostrapped sample (resampled_median) and add to the dataframe
	df_boostrapped_mean <- data.frame(Condition=levels(factor(df_raw_data$Condition)),
	                       resampled_mean=tapply(df_raw_data$value, 
	                       df_raw_data$Condition, boot_mean), id=i)
	
	#Add the new median to a datafram that collects all the resampled median values
	df_all_new_mean <- bind_rows(df_all_new_mean, df_boostrapped_mean)
}

#Calculate the confidence interval of the boostrapped medians, based on percentiles, and add to the dataframe that summarizes the data
df_summary$ci_mean_hi <- tapply(df_all_new_mean$resampled_mean, 
              df_all_new_mean$Condition, quantile, probs=upper_percentile)

df_summary$ci_mean_lo <- tapply(df_all_new_mean$resampled_mean, df_all_new_mean$Condition, 
              quantile, probs=lower_percentile)

#####################################################################################
############ Generate dataframe with differences between the Median and 95%CI  ######

#Order by Condition
df_all_new_mean[order(df_all_new_mean$Condition),]

#Convert the boostrapped dataset from long to wide format
df_spread <- spread(df_all_new_mean, key=Condition, value=resampled_mean)





#Subtract the Column with header "wt" from the other columns and move these 'differences' into a new dataframe
#df_spread_differences <- df_spread[,2:ncol(df_spread)] - df_spread[,"wt"]
df_spread_differences <- df_spread[,2:ncol(df_spread)] - df_spread[,"AsaiapHM4"]

#Convert the dataframe with differences between medians into a long format
df_differences <- gather(df_spread_differences, Condition, diff_mean)

#Calculate the summary of the differences and put these values in a dataframe 'df_diff_summary'
df_diff_summary <- data.frame(Condition=levels(factor(df_raw_data$Condition)), 
                    diff_mean=tapply(df_differences$diff_mean, df_differences$Condition, mean))

#Determine the CI of the differences, based on percentiles, and add to the dataframe that summarizes the differences
df_diff_summary$ci_mean_hi <- tapply(df_differences$diff_mean, df_differences$Condition, quantile, probs=upper_percentile)
df_diff_summary$ci_mean_lo <- tapply(df_differences$diff_mean, df_differences$Condition, quantile, probs=lower_percentile)


# #Sort dataframe Conditions according to median value
# df_summary$Condition <- reorder (df_summary$Condition, df_summary$median)
# df_raw_data$Condition <- reorder (df_raw_data$Condition, df_raw_data$value, median)
# 
# #Order according to original order, i.e. as read from csv (this is achieved by levels=df_raw_data$Condition)
# df_summary$Condition <- factor(df_summary$Condition, levels=df_raw_data$Condition)
# df_raw_data$Condition<- factor(df_raw_data$Condition, levels=df_raw_data$Condition)
# df_diff_summary$Condition <- factor(df_diff_summary$Condition, levels=df_raw_data$Condition)
# df_differences$Condition  <- factor(df_differences$Condition, levels=df_raw_data$Condition)
# df_all_new_medians$Condition <- factor(df_all_new_medians$Condition, levels=df_raw_data$Condition)

##################### Generate plots for output ###################

######### PLOT with raw data + median and Confidence Interval ########
meanplot <- ggplot(df_summary, aes(x = Condition, y = mean),width=0.02)+
  geom_boxplot(data=df_raw_data,aes(x=Condition, y=value),width=0.2, alpha=0.3, color=c("blue", "red"))+
  geom_quasirandom(data = df_raw_data, aes(x=Condition, y=value),
                   varwidth = TRUE, cex=2, alpha=0.4)+theme_light(base_size = 12) + 
  theme(panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill=NA)) + 
  theme(axis.text.y = element_text(size=8))+theme(legend.position="none")+
  scale_x_discrete(labels=c(expression(Asaia^{italic(pHM4)}),expression(Asaia^{italic(WSP)})))+
  theme(panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill=NA)) + 
  labs(x="", y="Observed raw data of IL6")

effectplot <- ggplot(df_differences, aes(x = Condition, y = diff_mean))+
  geom_flat_violin(aes(x=Condition),position = position_nudge(x = 0, y = 0),
                   color=NA,scale = "width", alpha =0.7,data = df_differences)+ 
  geom_linerange(aes(y= diff_mean, ymin = ci_mean_lo, ymax = ci_mean_hi), 
                 color="black", size =1, data = df_diff_summary)+
  geom_point(aes(x=df_diff_summary[2,1],y = df_diff_summary[2,2]), shape = 20,color = "black",fill=NA,size = 3,data = df_diff_summary)+
  geom_hline(yintercept = 0, col = "black", size = .5) +
  theme_light(base_size = 12)+
  scale_x_discrete("", labels=c("", expression(Asaia^{italic(WSP)} - Asaia^{italic(phM4)})))+
  theme(panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill=NA)) + 
  theme(legend.position="none")+labs(y="mean of differences")+scale_y_continuous(limits=c(-1000,4000))



(meanplot)/(effectplot)


#  geom_jitter(data = df_raw_data, aes(x=Condition, y=value), position=position_jitter(0.3), cex=1, alpha=0.4)+
  #geom_errorbar(data=df_summary, aes(x=Condition, ymin=mean, ymax=mean), 
 #color="black",width=.8, size=1, alpha=1)+
# geom_linerange(aes(ymin = ci_median_lo, ymax = ci_median_hi), color="black", size =1, data = df_summary)+
# geom_point(aes(y = median), shape = 21,color = "black",fill=NA,size = 3, data = df_summary)+
  #ggtitle(paste("Data and median"))  +
  # ylab("") + xlab("")+
  
  #coord_flip()+
  #NULL


# ######### PLOT with all medians from bootstrapping ########
# Bootstrapplot <- ggplot(df_summary, aes(x=Condition, color=Condition)) +
#   geom_quasirandom(data = df_all_new_medians, aes(x=Condition, y=resampled_median), 
#                    varwidth = TRUE, cex=1, alpha=0.1)+
#   
# #  geom_errorbar(data=df_summary, aes(x=Condition, ymin=median, ymax=median), color="black",width=.8, size=1, alpha=1)+
#   
#   geom_point(aes(y = mean), shape = 21,color = "black",fill=NA,size = 4, data = df_summary)+
#   ggtitle(paste("Bootstrapped Medians"))  +
#   ylab("Area [um^2]")+theme_light(base_size = 14) + 
#   theme(panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill=NA)) + theme(axis.text.y = element_text(size=12))+theme(legend.position="none")+
#   #Set the axis limits here - autoscale if not set
#   #ylim(-0.02,0.3)+
#   theme(axis.text.y = element_text(size=0))+
#   labs(x = NULL)+
#   coord_flip()+
#   NULL
# 
# ######### PLOT with all differences from median ########
# DiffBootplot <- ggplot(df_differences, aes(x=Condition, color=Condition)) +
#   geom_quasirandom(data = df_differences, aes(x=Condition, y=diff_median), varwidth = TRUE, cex=1, alpha=0.1)+
# # geom_errorbar(data=df_summary, aes(x=Condition, ymin=median, ymax=median), color="black",width=.8, size=1, alpha=1)+
#   geom_point(aes(y = diff_median), shape = 21,color = "black",fill=NA,size = 4,data = df_diff_summary)+
#   geom_hline(yintercept = 0, col = "black", size = .5) +
#   ggtitle(paste("Differences between medians"))  +
#   ylab("Difference [um^2]")+theme_light(base_size = 14) + 
#   theme(panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black", fill=NA)) + theme(axis.text.y = element_text(size=12))+theme(legend.position="none")+
#   #Set the axis limits here - autoscale if not set
#   #ylim(-0.02,0.3)+
#   theme(axis.text.y = element_text(size=0))+
#   labs(x = NULL)+
#   coord_flip()+
#   NULL


########### PLOT with distribution of differences, median difference and Confidence Interval ########
  
# g <- arrangeGrob(meanplot, #Bootstrapplot, DiffBootplot, 
#                  effectplot, ncol=1, nrow=2)
# 
# plot(g)
# 
# g+rotate()
# 
# 
# #Uncomment to save the plot
# # ggsave(file="Output-differences.pdf", g)
