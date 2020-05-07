library(readxl)
library(tidyverse)
library(dabestr)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(BayesPostEst)
library(patchwork)
library(gridExtra)
library(MKinfer)
library(ggpubr)
library(pscl)
rm(list=ls())
options(scipen = .999)
dati <- read_excel("Dati AsaiaWSP-Leishmania_MODIFICATO.xlsx")
#nleish<-read_excel("nparinmacro.xlsx", sheet = "Foglio3")
#rateinf<-read_excel("nminf.xlsx")
#cito<-read_excel("cito.xlsx")
#Cytokine####

########TABELLA DESCRITTIVA###############



# 
# tab2<-dati %>% 
#   pivot_longer(cols=4:12, names_to = "Parametro", values_to = "value") %>% 
#   group_by(gruppo,time,Parametro) %>% 
#   summarise(n=n(), mean=round(mean(value, na.rm=TRUE),2), sd=round(sd(value,na.rm=TRUE),2),
#             qnt_25  = round(quantile(value, probs= 0.25,na.rm = TRUE),2),
#             median  = round(quantile(value, probs= 0.50,na.rm = TRUE),2),
#             qnt_75  = round(quantile(value, probs= 0.75,na.rm = TRUE),2)) %>% 
#   data.frame()
#   
#   
#   
# write.table(tab2, file="descrittiva.csv")
#   
# 
# 
# dati %>% 
#   pivot_longer(cols=4:12, names_to = "Parametro", values_to = "value") %>% 
#   filter(Parametro=="iNOS") %>% 
#   group_by(gruppo,time,Parametro) %>%
#   filter(gruppo=="Leishmania") %>% 
#   summarise(  mean= mean(value, na.rm=TRUE), sd=sd(value,na.rm=TRUE),
#             qnt_25  =  quantile(value, probs= 0.25,na.rm = TRUE),
#             median  =  quantile(value, probs= 0.50,na.rm = TRUE),
#             qnt_75  =  quantile(value, probs= 0.75,na.rm = TRUE)) %>% 
#   data.frame()
# 
# options(scipen = 999)
# tab <- read_excel("tab.xlsx", sheet = "Foglio2", col_types = c("numeric", "numeric"))
# 
# x<-data.frame("var"=round((tab$difference/tab$`control_group (mean)`)*100,2))
# write.table(x, file="var.csv")






  #IL6#######
IL6<-dati[,2:4]


 
######bootstrap t-test IL6  24h####  
x=IL6 %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 



boot.t.test(IL6 ~ gruppo, data=x)

xx=IL6 %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 
xx%>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(IL6),
            sd=sd(IL6)) %>% 
  
  as.data.frame()

boot.t.test(iL6 ~ gruppo, data=xx)

z=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Med","AsaiapHM4")) %>% 
  mutate(gruppo=factor(gruppo, levels=c("AsaiapHM4","Med"))) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6)
  


 boot.t.test(IL6 ~ gruppo, data=z)

z2=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z2)

z3=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z3)

z4=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z4)


#####grafici IL6 24h#####
il6A24h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=factor(gruppo)) %>% 
  dabest(gruppo, IL6, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )


#ggplot(il6A24h, aes(x=gruppo, y=IL6))


plot(il6A24h,rawplot.type = "sinaplot", group.summaries = "mean_sd", float.contrast = FALSE,
            rawplot.ylabel = "Observed data of  IL6 pg/ml ",
            effsize.ylabel = "Unpaired mean difference",
            axes.title.fontsize = 12,
            rawplot.ylim = c(0,8000),
            effsize.ylim = c(-1000,3500))


#####new graph


  
px<-il6A24h[["data"]]
px<-px[,-2]

z<-split(px, px$gruppo)

p<-cbind(z$AsaiapHM4[,2] ,z$AsaiaWSP[,2])
names(p)<-c("AsaiapHM4","AsaiaWSP")
write.csv(p, file="prova.csv")

write.csv(px, file="prova2.csv")



il6B24h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  # mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, IL6, 
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )




plot(il6B24h,rawplot.type = "sinaplot", group.summaries = "mean_sd", float.contrast = FALSE,
            rawplot.ylabel = "Observed data of IL6 pg/ml ",
            effsize.ylabel = "Unpaired mean difference",
            axes.title.fontsize = 12,
         rawplot.ylim = c(0,8000),
         effsize.ylim = c(-2000,3500))
dev.off()












il6C<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  #mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, IL6,
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE, 
  )

p3<-plot(il6C,rawplot.type = "sinaplot", group.summaries = "mean_sd", float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)




###PLOT 1 IL6 24h####

(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(IL6) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
# between the sample mean Log10(IL6) in AsaiaWSP (3.75 pg/ml) and AsaiapHM4 (3.60 pg/ml) was 0.15 pg/ml, with a
# bootstraped 95% confidence interval from 0.07 to 0.22 pg/ml; the bootstraped 
# Welch two-sample t-test P value was 0.0008.
# 
# In panel B is reported the unpaired mean comparison of Log10(IL6) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(IL6) in AsaiaWSP L (3.36 pg/ml) and AsaiapHM4 (3.22 pg/ml) was 0.14 pg/ml, with a
# bootstraped 95% confidence interval from -0.20 to 0.50 pg/ml; the bootstraped 
# Welch two-sample t-test P value was 0.44. ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))


###PLOT 2 IL6 multiple contrast 24h####
(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(IL6) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(IL6) in AsaiapHM4 (3.60 pg/ml) and Med (0.99 pg/ml) was 2.61 pg/ml, with a
bootstraped 95% confidence interval from 2.25 to 3.01 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

2) The difference between the sample mean Log10(IL6) in AsaiaWSP (3.74 pg/ml) and Med (0.99 pg/ml) was 2.75 pg/ml, with a
bootstraped 95% confidence interval from 2.42 to 3.16 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

3) The difference between the sample mean Log10(IL6) in Leishmania (0.79 pg/ml) and Med (0.99 pg/ml) was -0.21 pg/ml, with a
bootstraped 95% confidence interval from -0.61 to 0.24 pg/ml; the bootstraped Welch two-sample t-test P value was 0.41. 

4) The difference between the sample mean Log10(IL6) in LPS (4.35 pg/ml) and Med (0.99 pg/ml) was 3.36 pg/ml, with a
bootstraped 95% confidence interval from 2.58 to 4.33 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))


##################################################################################################################
##################################################################################################################

######bootstrap t-test IL6  48h####  
x=IL6 %>% 
  #select(-time) %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

x %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(IL6),
            sd=sd(IL6)) %>% 
  as.data.frame()

boot.t.test(IL6 ~ gruppo, data=x)

xx=IL6 %>% 
  #select(-time) %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

xx %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(IL6),
            sd=sd(IL6)) %>% 
  as.data.frame()

boot.t.test(IL6 ~ gruppo, data=xx)

z=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z)

z2=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z2)

z3=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z3)

z4=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(IL6 ~ gruppo, data=z4)

#####grafici IL6 48h#####
il6A48h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, IL6, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(il6A48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  IL6 pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,8000),
         effsize.ylim = c(-2000,5500)
)

il6B48h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, IL6, 
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(il6B48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of IL6 pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,8000),
         effsize.ylim = c(-2500,5500))

il6C48h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  # mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, IL6, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(il6C,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)

###PLOT 1 IL6 48h####

(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(IL6) between Asaia WSP vs AsaiaphM4 as control group at 48h: The difference 
# between the sample mean Log10(IL6) in AsaiaWSP (3.58 pg/ml) and AsaiapHM4 (3.20 pg/ml) was 0.38 pg/ml, with a
# bootstraped 95% confidence interval from 0.19 to 0.58 pg/ml; the bootstraped 
# Welch two-sample t-test P value was 0.0090.
# 
# In panel B is reported the unpaired mean comparison of Log10(IL6) Asaia WSP L vs AsaiaphM4 L as control group at 48h: The difference 
# between the sample mean Log10(IL6) in AsaiaWSP L (3.37 pg/ml) and AsaiapHM4 (3.09 pg/ml) was 0.29 pg/ml, with a
# bootstraped 95% confidence interval from -0.18 to 0.77 pg/ml; the bootstraped 
# Welch two-sample t-test P value was 0.25. ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))
  ###PLOT 2 IL6 48h####

(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(IL6) between different Treatment Group vs Med as control group at 48h: 
         
1) The difference between the sample mean Log10(IL6) in AsaiapHM4 (3.20 pg/ml) and Med (0.73 pg/ml) was 2.48 pg/ml, with a
bootstraped 95% confidence interval from 2.00 to 3.09 pg/ml; the bootstraped Welch two-sample t-test P value was 0.0026. 

2) The difference between the sample mean Log10(IL6) in AsaiaWSP (3.58 pg/ml) and Med (0.73 pg/ml) was 2.86 pg/ml, with a
bootstraped 95% confidence interval from 2.39 to 3.47 pg/ml; the bootstraped Welch two-sample t-test P value was 0.0008. 

3) The difference between the sample mean Log10(IL6) in Leishmania (0.80 pg/ml) and Med (0.72 pg/ml) was 0.07 pg/ml, with a
bootstraped 95% confidence interval from -0.43 to 0.72 pg/ml; the bootstraped Welch two-sample t-test P value was 0.79. 

4) The difference between the sample mean Log10(IL6) in LPS (4.33 pg/ml) and Med (0.73 pg/ml) was 3.61 pg/ml, with a
bootstraped 95% confidence interval from 2.59 to 4.86 pg/ml; the bootstraped Welch two-sample t-test P value was 0.0002.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

#####TNF#######

TNF<-dati[,c(2:3, 5)]


######bootstrap t-test TNF 24h####  
x=TNF %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

x %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(TNF),
            sd=sd(TNF)) %>% 
  as.data.frame()



boot.t.test(TNF ~ gruppo, data=x)

xx=TNF %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

xx %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(TNF),
            sd=sd(TNF))%>% 
  as.data.frame()

boot.t.test(TNF ~ gruppo, data=xx)


z=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo, data=z)

z2=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo, data=z2)

z3=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
 # mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo,data=z3)

z4=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo, data=z4)



####grafici TNF 24h####
 A24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, TNF, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of TNF pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,3500)#,
         # effsize.ylim = c(-0.3,0.4)
)

B24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, TNF,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  TNF pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4000),
         effsize.ylim = c(-1000,2000)
)



tnf24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  #mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, TNF, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)



#####PLOT 1 TNF 24h#####

(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(TNF) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
# between the sample mean Log10(TNF) in AsaiaWSP (2.95 pg/ml) and AsaiapHM4 (2.87 pg/ml) was 0.08 pg/ml, with a
# bootstraped 95% confidence interval from -0.13 to 0.28 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.41.
# 
# In panel B is reported the unpaired mean comparison of Log10(TNF) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(TNF) in AsaiaWSP L (2.41 pg/ml) and AsaiapHM4 (2.29 pg/ml) was 0.11 pg/ml, with a
# bootstraped 95% confidence interval from -0.60 to 0.87 pg/ml; the bootstraped Welch two-sample t-test P value was 0.76. ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))

#### PLOT 2 TNF 24h multipe contrast####

(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(TNF) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(TNF) in AsaiapHM4 (2.88 pg/ml) and Med (0.76 pg/ml) was 2.11 pg/ml, with a
bootstraped 95% confidence interval from 1.84 to 2.41 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

2) The difference between the sample mean Log10(TNF) in AsaiaWSP (2.95 pg/ml) and Med (0.76 pg/ml) was 2.20 pg/ml, with a
bootstraped 95% confidence interval from 1.89 to 2.50 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

3) The difference between the sample mean Log10(TNF) in Leishmania (0.81 pg/ml) and Med (0.76 pg/ml) was 0.05 pg/ml, with a
bootstraped 95% confidence interval from -0.39 to 0.41 pg/ml; the bootstraped Welch two-sample t-test P value was 0.81. 

4) The difference between the sample mean Log10(TNF) in LPS (3.43 pg/ml) and Med (0.76 pg/ml) was 2.66 pg/ml, with a
bootstraped 95% confidence interval from 2.40 to 2.95 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))



#####################################################################################################################
####################################################################################################################
#################################################################################################################

####bootstrap TNF 48h#######

x=TNF %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 
x %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(TNF),
            sd=sd(TNF)) %>% 
  as.data.frame()


boot.t.test(TNF ~ gruppo, data=x)

xx=TNF %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

xx %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(TNF),
            sd=sd(TNF)) %>% 
  as.data.frame()


boot.t.test(tnf ~ gruppo, data=xx)

z=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo, data=z)

z2=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo, data=z2)

z3=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo,data=z3)

z4=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  drop_na(TNF) 

boot.t.test(TNF ~ gruppo, data=z4)

#####grafici TNF 48h#####

A24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, TNF, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of TNF pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,450),
         effsize.ylim = c(-150,150)
)

B24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, TNF,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of TNF pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,450),
         effsize.ylim = c(-150,300)
)



tnf48h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  #mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, TNF, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)


#####Plot 1 TNF 48h#####
(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(TNF) between Asaia WSP vs AsaiaphM4 as control group at 48h: The difference 
# between the sample mean Log10(TNF) in AsaiaWSP (2.17 pg/ml) and AsaiapHM4 (2.13 pg/ml) was 0.05 pg/ml, with a
# bootstraped 95% confidence interval from -0.14 to 0.22 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.68.
# 
# In panel B is reported the unpaired mean comparison of Log10(TNF) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(TNF) in AsaiaWSP L (2.23 pg/ml) and AsaiapHM4 (2.13 pg/ml) was 0.10 pg/ml, with a
# bootstraped 95% confidence interval from -0.13 to 0.38 pg/ml; the bootstraped Welch two-sample t-test P value was 0.45. ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))








#####Plot 2 TNF 48h multipe contrast####
(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(TNF) between different Treatment Group vs Med as control group at 48h: 
         
1) The difference between the sample mean Log10(TNF) in AsaiapHM4 (2.13 pg/ml) and Med (0.89 pg/ml) was 1.23 pg/ml, with a
bootstraped 95% confidence interval from 0.46 to 2.13 pg/ml; the bootstraped Welch two-sample t-test P value was 0.040. 

2) The difference between the sample mean Log10(TNF) in AsaiaWSP (2.17 pg/ml) and Med (0.89 pg/ml) was 1.28 pg/ml, with a
bootstraped 95% confidence interval from 0.49 to 2.14 pg/ml; the bootstraped Welch two-sample t-test P value was 0.040. 

3) The difference between the sample mean Log10(TNF) in Leishmania (1.17 pg/ml) and Med (0.89 pg/ml) was 0.28 pg/ml, with a
bootstraped 95% confidence interval from -0.64 to 1.27 pg/ml; the bootstraped Welch two-sample t-test P value was 0.053. 

4) The difference between the sample mean Log10(TNF) in LPS (3.02 pg/ml) and Med (0.89 pg/ml) was 2.13 pg/ml, with a
bootstraped 95% confidence interval from 1.29 to 3.11 pg/ml; the bootstraped Welch two-sample t-test P value was 0.020.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))

#################################################################################################################
#################################################################################################################

#####IL1beta##################
IL1b<-dati[,c(2:3,6)]
#####bootstrap IL1beta 24h#####
x=IL1b %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  drop_na(IL1beta) 
x %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(IL1beta),
            sd=sd(IL1beta)) %>% 
  as.data.frame()

boot.t.test(IL1beta ~ gruppo, data=x)

xx=IL1b %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  drop_na(IL1beta) 

xx %>% 
  group_by(time,gruppo ) %>% 
  summarize(Mean=mean(IL1beta),
            sd=sd(IL1beta)) %>% 
  as.data.frame()

boot.t.test(IL1beta ~ gruppo, data=xx)

z=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  drop_na(IL1beta) 

boot.t.test(IL1beta ~ gruppo,data=z)

z2=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  drop_na(IL1beta) 

boot.t.test(IL1beta ~ gruppo, data=z2)

z3=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  drop_na(IL1beta) 

boot.t.test(IL1beta ~ gruppo,data=z3)

z4=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  drop_na(IL1beta) 

boot.t.test(IL1beta ~ gruppo, data=z4)


#####grafici IL1beta 24h######
A24h<-IL1b %>%
  filter(!is.na(IL1beta)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  dabest(gruppo, IL1beta, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  IL1 beta pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,450)
         # effsize.ylim = c(-0.3,0.4)
)

B24h<-IL1b %>%
  filter(!is.na(IL1beta)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il1b=log10(IL1beta)) %>% 
  dabest(gruppo, IL1beta,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of IL1 beta pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,450)
         # effsize.ylim = c(-0.1,0.1)
)



ilb24h<-IL1b %>%
  filter(!is.na(IL1beta)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  #mutate(il1b=log10(IL1beta)) %>% 
  dabest(gruppo, IL1beta, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL1 beta) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)





#####PLOT 1 IL1 beta 24h#####
(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(IL1 beta) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
# between the sample mean Log10(IL1 beta) in AsaiaWSP (2.55 pg/ml) and AsaiapHM4 (2.47 pg/ml) was 0.08 pg/ml, with a
# bootstraped 95% confidence interval from -0.004 to 0.16 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.080.
# 
# In panel B is reported the unpaired mean comparison of Log10(IL1 beta) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(IL1 beta) in AsaiaWSP L (2.48 pg/ml) and AsaiapHM4 (2.47 pg/ml) was 0.18 pg/ml, with a
# bootstraped 95% confidence interval from -0.02 to 0.06 pg/ml; the bootstraped Welch two-sample t-test P value was 0.373. ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))


#####PLOT 2 IL1 beta 24h#######
(p3)+plot_annotation(caption=
 
"In this figure is reported the unpaired mean comparison of Log10(IL1 beta) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(IL1 beta) in AsaiapHM4 (2.47 pg/ml) and Med (1.14 pg/ml) was 1.32 pg/ml, with a
bootstraped 95% confidence interval from 1.16 to 1.50 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

2) The difference between the sample mean Log10(IL1 beta) in AsaiaWSP (2.55 pg/ml) and Med (1.14 pg/ml) was 1.41 pg/ml, with a
bootstraped 95% confidence interval from 1.25 to 1.57 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

3) The difference between the sample mean Log10(IL1 beta) in Leishmania (1.68 pg/ml) and Med (1.14  pg/ml) was 0.54 pg/ml, with a
bootstraped 95% confidence interval from 0.37 to 0.71 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

4) The difference between the sample mean Log10(IL1 beta) in LPS (2.08 pg/ml) and Med (1.14  pg/ml) was 0.93 pg/ml, with a
bootstraped 95% confidence interval from 0.77 to 1.10 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))





####NO dati a 48h per IL1 beta
#########################################################################################################################
#########################################################################################################################
####IL12p40#####

IL12p40<-dati[,c(2:3,7)]

#####bootstrap IL12p40 24h#####
x=IL12p40 %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  drop_na(IL12p40) 

boot.t.test(IL12p40 ~ gruppo, data=x)

xx=IL12p40 %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  drop_na(IL12p40) 

boot.t.test(IL12p40  ~ gruppo, data=xx)

z=IL12p40 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  drop_na(IL12p40) 

boot.t.test(IL12p40  ~ gruppo,data=z)

z2=IL12p40 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  drop_na(IL12p40) 

boot.t.test(IL12p40  ~ gruppo, data=z2)

z3=IL12p40 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  drop_na(IL12p40) 

boot.t.test(IL12p40  ~ gruppo,data=z3)

z4=IL12p40 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  drop_na(IL12p40) 

boot.t.test(IL12p40  ~ gruppo, data=z4)

#####grafici IL12p40 24h######

A24h<-IL12p40 %>%
  filter(!is.na(IL12p40)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  dabest(gruppo, IL12p40, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of IL12p40 pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12#,
         # rawplot.ylim = c(0,4.5),
         # effsize.ylim = c(-0.35,0.65)
)

B24h<-IL12p40 %>%
  filter(!is.na(IL12p40)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(il12b=log10(IL12p40)) %>% 
  dabest(gruppo, IL12p40,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  IL12p40 pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
          axes.title.fontsize = 12#,
         # rawplot.ylim = c(0,4.5),
         # effsize.ylim = c(-0.35,0.65)
)



il12p40<-IL12p40 %>%
  filter(!is.na(IL12p40)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
 # mutate(il12b=log10(IL12p40)) %>% 
  dabest(gruppo, IL12p40, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL12p40) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)





####PLOT 1 IL12p40 24h####
(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(IL12p40) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
# between the sample mean Log10(IL12p40) in AsaiaWSP (3.76 pg/ml) and AsaiapHM4 (3.63 pg/ml) was 0.134 pg/ml, with a
# bootstraped 95% confidence interval from -0.15 to 0.42 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.400.
# 
# In panel B is reported the unpaired mean comparison of Log10(IL12p40) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(IL12p40) in AsaiaWSP L (3.65 pg/ml) and AsaiapHM4 (3.50 pg/ml) was 0.15 pg/ml, with a
# bootstraped 95% confidence interval from -0.16 to 0.45 pg/ml; the bootstraped Welch two-sample t-test P value was 0.322 ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))

####PLOT 2 IL12p40 24h#####
(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(IL12p40) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(IL12p40) in AsaiapHM4 (3.63 pg/ml) and Med (0.72 pg/ml) was 2.91 pg/ml, with a
bootstraped 95% confidence interval from 2.63 to 3.23 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

2) The difference between the sample mean Log10(IL12p40) in AsaiaWSP (3.74 pg/ml) and Med (0.72 pg/ml) was 3.04 pg/ml, with a
bootstraped 95% confidence interval from 2.77 to 3.39 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

3) The difference between the sample mean Log10(IL12p40) in Leishmania (0.50 pg/ml) and Med (0.72  pg/ml) was -0.22 pg/ml, with a
bootstraped 95% confidence interval from -0.50 to 0.15 pg/ml; the bootstraped Welch two-sample t-test P value was 0.315. 

4) The difference between the sample mean Log10(IL12p40) in LPS (3.46 pg/ml) and Med (0.72  pg/ml) was 2.74 pg/ml, with a
bootstraped 95% confidence interval from 2.37 to 3.16 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))





#########################################################################################################################
#########################################################################################################################

####ROS#####

ROS<-dati[,c(2:3,8)]
#####bootstrap ros 24h#####
x=ROS %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(ros=log10(ROS)) %>% 
  drop_na(ROS) 

boot.t.test(ROS ~ gruppo, data=x)

xx=ROS %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
 # mutate(ros=log10(ROS)) %>% 
  drop_na(ROS) 

boot.t.test(ROS ~ gruppo, data=xx)

z=ROS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(ros=log10(ROS)) %>% 
  drop_na(ROS) 

boot.t.test(ROS ~ gruppo,data=z)

z2=ROS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
 # mutate(ros=log10(ROS)) %>% 
  drop_na(ROS) 

boot.t.test(ROS ~ gruppo, data=z2)

z3=ROS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(ros=log10(ROS)) %>% 
  drop_na(ROS) 

boot.t.test(ROS ~ gruppo, data=z3)

z4=ROS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(ros=log10(ROS)) %>% 
  drop_na(ROS) 

boot.t.test(ROS ~ gruppo, data=z4)




  
#####grafici ros 24h######
A24h<-ROS %>%
  filter(!is.na(ROS)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(ros=log10(ROS)) %>% 
  dabest(gruppo, ROS, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  ROS pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         # rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-2000,12500)
)

B24h<-ROS %>%
  filter(!is.na(ROS)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(ros=log10(ROS)) %>% 
  dabest(gruppo, ROS, 
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of ROS pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12#,
         # rawplot.ylim = c(0,4.5),
         # effsize.ylim = c(-0.35,0.65)
)



ros<-ROS %>%
  filter(!is.na(ROS)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(ros=log10(ROS)) %>% 
  dabest(gruppo, ROS,  
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(ROS) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)


####PLOT 1 ros 24h####
(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(ROS) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
# between the sample mean Log10(ROS) in AsaiaWSP (4.00 pg/ml) and AsaiapHM4 (3.54 pg/ml) was 0.46 pg/ml, with a
# bootstraped 95% confidence interval from 0.20 to 0.71 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.0004.
# 
# In panel B is reported the unpaired mean comparison of Log10(ROS) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(ROS) in AsaiaWSP L (3.87 pg/ml) and AsaiapHM4 (3.81 pg/ml) was 0.07 pg/ml, with a
# bootstraped 95% confidence interval from -0.14 to 0.35 pg/ml; the bootstraped Welch two-sample t-test P value was 0.598 ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))


####PLOT 2 ros 24h#####
(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(ROS) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(ROS) in AsaiapHM4 (3.54 pg/ml) and Med (3.51 pg/ml) was 0.02 pg/ml, with a
bootstraped 95% confidence interval from -0.28 to 0.26 pg/ml; the bootstraped Welch two-sample t-test P value was 0.897. 

2) The difference between the sample mean Log10(ROS) in AsaiaWSP (4.00 pg/ml) and Med (3.51 pg/ml) was 0.49 pg/ml, with a
bootstraped 95% confidence interval from 0.30 to 0.66 pg/ml; the bootstraped Welch two-sample t-test P value was 0.0002. 

3) The difference between the sample mean Log10(ROS) in Leishmania (3.79 pg/ml) and Med (3.51  pg/ml) was 0.28 pg/ml, with a
bootstraped 95% confidence interval from 0.11 to 0.47 pg/ml; the bootstraped Welch two-sample t-test P value was 0.008. 

4) The difference between the sample mean Log10(ROS) in LPS (4.82 pg/ml) and Med (3.51  pg/ml) was 1.31 pg/ml, with a
bootstraped 95% confidence interval from 0.96 to 1.59 pg/ml; the bootstraped Welch two-sample t-test P value was 0.0002.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))




####################################################################################################################
####################################################################################################################

####Nitriti############
Nitriti<-dati[,c(2:3,9)]
#####bootstrap NO 24h#####
x=Nitriti %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=x)

xx=Nitriti %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=xx)

z=Nitriti %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
 # mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo,data=z)

z2=Nitriti %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=z2)

z3=Nitriti %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=z3)

z4=Nitriti %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=z4)



#####grafici NO 24h######
A24h<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  dabest(gruppo, Nitriti, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of NO μM ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
          rawplot.ylim = c(0,120)
         # effsize.ylim = c(-0.35,0.85)
)

B24h<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  dabest(gruppo, Nitriti,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  NO μM ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,140)#,
         # effsize.ylim = c(-0.35,0.65)
)



No24h<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  dabest(gruppo, Nitriti,  
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(NO) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)



####PLOT 1 NO 24h####
(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(NO) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
# between the sample mean Log10(NO) in AsaiaWSP (1.83 pg/ml) and AsaiapHM4 (1.82 pg/ml) was 0.01 pg/ml, with a
# bootstraped 95% confidence interval from -0.13 to 0.14 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.894
# 
# In panel B is reported the unpaired mean comparison of Log10(ROS) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
# between the sample mean Log10(NO) in AsaiaWSP L (1.87 pg/ml) and AsaiapHM4 (1.88 pg/ml) was -0.01 pg/ml, with a
# bootstraped 95% confidence interval from -0.15 to 0.12 pg/ml; the bootstraped Welch two-sample t-test P value was 0.598 ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))


####PLOT 2 NO 24h####
(p3)+plot_annotation(caption=
"In this figure is reported the unpaired mean comparison of Log10(NO) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(NO) in AsaiapHM4 (1.82 pg/ml) and Med (1.75 pg/ml) was 0.07 pg/ml, with a
bootstraped 95% confidence interval from -0.09 to 0.24 pg/ml; the bootstraped Welch two-sample t-test P value was 0.455. 

2) The difference between the sample mean Log10(NO) in AsaiaWSP (1.83 pg/ml) and Med (1.75 pg/ml) was 0.08 pg/ml, with a
bootstraped 95% confidence interval from -0.09 to 0.25 pg/ml; the bootstraped Welch two-sample t-test P value was 0.376. 

3) The difference between the sample mean Log10(NO) in Leishmania (1.78 pg/ml) and Med (1.75 pg/ml) was 0.02 pg/ml, with a
bootstraped 95% confidence interval from -0.13 to 0.19 pg/ml; the bootstraped Welch two-sample t-test P value was 0.787. 

4) The difference between the sample mean Log10(NO) in LPS (1.90 pg/ml) and Med (1.75  pg/ml) was 0.15 pg/ml, with a
bootstraped 95% confidence interval from -0.04 to 0.36 pg/ml; the bootstraped Welch two-sample t-test P value was 0.199.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))


######################################################################
#####bootstrap NO 48h#####
x=Nitriti %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=x)

xx=Nitriti %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=xx)

z=Nitriti %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo,data=z)

z2=Nitriti %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=z2)

z3=Nitriti %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=z3)

z4=Nitriti %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  drop_na(Nitriti) 

boot.t.test(Nitriti ~ gruppo, data=z4)




  #####grafici NO 48h######
A48h<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  dabest(gruppo, Nitriti, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  NO μM ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,160),
         effsize.ylim = c(-90,90)
)

B48h<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  dabest(gruppo, Nitriti,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of  NO μM",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,160),
         effsize.ylim = c(-90,90)
)



No48h<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(nit=log10(Nitriti)) %>% 
  dabest(gruppo, Nitriti,  
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(NO) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)

####PLOT 1 NO 48h####
(p1|p2)+plot_annotation(tag_levels = 'A', 
#                         caption = "
# 
# In panel A is reported the unpaired mean comparison of Log10(NO) between Asaia WSP vs AsaiaphM4 as control group at 48h: The difference 
# between the sample mean Log10(NO) in AsaiaWSP (2.11 pg/ml) and AsaiapHM4 (1.89 pg/ml) was 0.23 pg/ml, with a
# bootstraped 95% confidence interval from 0.13 to 0.32 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.0004
# 
# In panel B is reported the unpaired mean comparison of Log10(ROS) Asaia WSP L vs AsaiaphM4 L as control group at 48h: The difference 
# between the sample mean Log10(NO) in AsaiaWSP L (2.10 pg/ml) and AsaiapHM4 (2.04 pg/ml) was 0.06 pg/ml, with a
# bootstraped 95% confidence interval from -0.05 to 0.16 pg/ml; the bootstraped Welch two-sample t-test P value was 0.318 ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))
####PLOT 2 NO 48h####
(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(NO) between different Treatment Group vs Med as control group at 48h: 
         
1) The difference between the sample mean Log10(NO) in AsaiapHM4 (1.89 pg/ml) and Med (1.88 pg/ml) was 0.01 pg/ml, with a
bootstraped 95% confidence interval from -0.10 to 0.29 pg/ml; the bootstraped Welch two-sample t-test P value was 0.914. 

2) The difference between the sample mean Log10(NO) in AsaiaWSP (2.12 pg/ml) and Med (1.88 pg/ml) was 0.23 pg/ml, with a
bootstraped 95% confidence interval from 0.15 to 0.29 pg/ml; the bootstraped Welch two-sample t-test P value was 0.003. 

3) The difference between the sample mean Log10(NO) in Leishmania (1.88 pg/ml) and Med (1.88 pg/ml) was 0.00 pg/ml, with a
bootstraped 95% confidence interval from -0.06 to 0.06 pg/ml; the bootstraped Welch two-sample t-test P value was 0.94. 

4) The difference between the sample mean Log10(NO) in LPS (2.19 pg/ml) and Med (1.88  pg/ml) was 0.15 pg/ml, with a
bootstraped 95% confidence interval from -0.04 to 0.36 pg/ml; the bootstraped Welch two-sample t-test P value was 0.035.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))
###################################################################################################################
###################################################################################################################
#####iNOS<-dati[,c(2:3,10)]
#####bootstrap NO 24h#####
x=iNOS %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  drop_na(inos) 

boot.t.test(inos ~ gruppo, data=x)

xx=iNOS %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  drop_na(inos) 

boot.t.test(inos ~ gruppo, data=xx)

z=iNOS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  drop_na(inos) 

boot.t.test(inos ~ gruppo,data=z)

z2=iNOS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  drop_na(inos) 

boot.t.test(inos ~ gruppo, data=z2)

z3=iNOS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  drop_na(inos) 

boot.t.test(inos ~ gruppo, data=z3)

z4=iNOS %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(inos=log10(iNOS)) %>% 
  drop_na(inos) 

boot.t.test(inos ~ gruppo, data=z4)


A24h<-iNOS %>%
  filter(!is.na(iNOS)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  dabest(gruppo, inos, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(iNOS) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(-1.8,0.3),
         effsize.ylim = c(-0.35,1.6)
)

B24h<-iNOS %>%
  filter(!is.na(iNOS)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  dabest(gruppo, inos,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(iNOS) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         #rawplot.ylim = c(-1,0.2),
         effsize.ylim = c(-0.1,1.2)
)


C24h<-iNOS %>%
  filter(!is.na(iNOS)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(inos=log10(iNOS)) %>% 
  dabest(gruppo, inos, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(NO) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)

#################################################################################################################
################################################################################################################
IL10<-dati[,c(2:3,11)]
#####bootstrap NO 24h#####
x=IL10 %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  drop_na(il10) 

boot.t.test(il10 ~ gruppo, data=x)

xx=IL10 %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  drop_na(il10) 

boot.t.test(il10 ~ gruppo, data=xx)

z=IL10 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  drop_na(il10) 

boot.t.test(il10 ~ gruppo,data=z)

z2=IL10 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  drop_na(il10) 

boot.t.test(il10 ~ gruppo, data=z2)

z3=IL10 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(il10=log10(IL10)) %>% 
  drop_na(il10) 

boot.t.test(il10 ~ gruppo, data=z3)

z4=IL10 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(il10=log10(IL10)) %>% 
  drop_na(il10) 

boot.t.test(il10 ~ gruppo, data=z4)


#####grafici IL10 24h######
A24h<-IL10 %>%
  filter(!is.na(IL10)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  dabest(gruppo, il10, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL10) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
rawplot.ylim = c(0,3.5),
effsize.ylim = c(-0.35,0.2))


B24h<-IL10 %>%
  filter(!is.na(IL10)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  dabest(gruppo, il10, 
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL10) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,3.5),
         effsize.ylim = c(-0.18,0.15)
)


C24h<-IL10 %>%
  filter(!is.na(IL10)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il10=log10(IL10)) %>% 
  dabest(gruppo, il10, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL10) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)



####PLOT 1 IL10 24h####
(p1|p2)+plot_annotation(tag_levels = 'A', 
                        caption = "

In panel A is reported the unpaired mean comparison of Log10(IL10) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
between the sample mean Log10(IL10) in AsaiaWSP (2.73 pg/ml) and AsaiapHM4 (2.84 pg/ml) was 0.11 pg/ml, with a
bootstraped 95% confidence interval from -0.23 to 0.06 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.402

In panel B is reported the unpaired mean comparison of Log10(IL10) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
between the sample mean Log10(IL10) in AsaiaWSP L (2.72 pg/ml) and AsaiapHM4 (2.71 pg/ml) was 0.01 pg/ml, with a
bootstraped 95% confidence interval from -0.06 to 0.10 pg/ml; the bootstraped Welch two-sample t-test P value was 0.870 ",
                        theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0), 
                                                              family = "Comic Sans MS"))
) &
  theme(plot.tag = element_text(size = 13))


####PLOT 2 IL10 24h#####
(p3)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of Log10(IL10) between different Treatment Group vs Med as control group at 24h: 
         
1) The difference between the sample mean Log10(IL10) in AsaiapHM4 (2.84 pg/ml) and Med (0.42 pg/ml) was 2.42 pg/ml, with a
bootstraped 95% confidence interval from 2.34 to 2.53 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16 

2) The difference between the sample mean Log10(IL10) in AsaiaWSP (2.73 pg/ml) and Med (0.42 pg/ml) was 2.31 pg/ml, with a
bootstraped 95% confidence interval from 2.17 to 2.50 pg/ml; the bootstraped Welch two-sample t-test P value was 0.0012. 

3) The difference between the sample mean Log10(IL10) in Leishmania (2.41 pg/ml) and Med (0.42 pg/ml) was 1.99 pg/ml, with a
bootstraped 95% confidence interval from 1.70 to 2.22 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16. 

4) The difference between the sample mean Log10(IL10) in LPS (2.66 pg/ml) and Med (0.42  pg/ml) was 2.24 pg/ml, with a
bootstraped 95% confidence interval from 2.10 to 2.37 pg/ml; the bootstraped Welch two-sample t-test P value was < 2.2e-16.",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))
####################################################################################################################
######################################################################################################################

######TABELLONE####

A<-il6C$result
B<-il6C48h$result
C<-tnf24h$result
D<-tnf48h$result
E<-ilb24h$result
F<-il12p40$result
G<-ros$result
H<-No24h$result
I<-No48h$result

tab<-data.frame(rbind(A,B,C,D,E,F,G,H,I))

 

tabella<-tab %>% 
  select(variable,test_group, control_group,difference,bca_ci_low,bca_ci_high) %>% 
  mutate(difference=round(difference,2),
         bca_ci_low=round(bca_ci_low,2),
         bca_ci_high=round(bca_ci_high,2))

write.table(tabella, file="tab.csv")



###Arginasi<-dati[,c(2:3,12)]
##############################################################################################################
##############################################################################################################

#####Fagocitosi####

fago <- read_excel("Dati fagocitosi.xlsx")
fago<-na.omit(fago)
names(fago)[4]<-"CFUml"

fago$time<-factor(fago$time, levels=c("1h", "2h", "24h"))

fago<-fago %>% 
  filter(time!="24h")
#####bootstrap fago 1h#####
x=fago %>% 
  filter(time=="1h") %>% 
  filter(gruppo %in% c("Asaia pHM4","Asaia WSP")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  #mutate(cfu=log10(CFUml)) %>% 
  drop_na(CFUml) 

boot.t.test(cfu ~ gruppo, data=x)

xx=fago %>% 
  filter(gruppo %in% c("Asaia pHM4","Asaia pHM4 + LPS")) %>% 
  filter(time=="1h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(cfu=log10(CFUml)) %>% 
  drop_na(cfu) 

boot.t.test(cfu ~ gruppo, data=xx)
#####grafici fago  1h######
A1h<-fago %>%
  filter(!is.na(CFUml)) %>% 
  filter(time=="1h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
 # mutate(cfu=log10(CFUml)) %>% 
  dabest(gruppo, CFUml, 
         idx = list(c("Asaia pHM4","Asaia WSP","Asaia pHM4 + LPS")), 
         paired = FALSE
  )

p1<-plot(A1h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of phagocytosis activity Log10(CFU/ml)",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 10#,
         # rawplot.ylim = c(0,3.5),
         # effsize.ylim = c(-0.35,0.2)
         )

####PLOT FAGO 1h#####
(p1)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of phagocytosis activity (PA) in Log10(UFC/ml) 
between different Treatment Group vs Asaia pHM4 as control group at 1h: 
         
1) The difference between the sample mean PA Log10(UFC/ml) in AsaiaWSP(5.98 UFC/ml) and AsaiapHM4 (5.80 UFC/ml) was 0.19 UFC/ml, with a
bootstraped 95% confidence interval from -0.24 to 0.56 UFC/ml; the bootstraped Welch two-sample t-test P value was 0.228 

2) The difference between the sample mean Log10(UFC/ml) in AsaiapHM4+LPS (5.97 UFC/ml) and AsaiapHM4 (5.80 UFC/ml) was 0.18 UFC/ml, with a
bootstraped 95% confidence interval from -0.30 to 0.53 UFC/ml; the bootstraped Welch two-sample t-test P value was 0.466.",

                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))






#####bootstrap fago 2h#####
x=fago %>% 
  filter(time=="2h") %>% 
  filter(gruppo %in% c("Asaia pHM4","Asaia WSP")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(cfu=log10(CFUml)) %>% 
  drop_na(cfu) 

boot.t.test(cfu ~ gruppo, data=x)

xx=fago %>% 
  filter(gruppo %in% c("Asaia pHM4","Asaia pHM4 + LPS")) %>% 
  filter(time=="2h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(cfu=log10(CFUml)) %>% 
  drop_na(cfu) 

boot.t.test(cfu ~ gruppo, data=xx)





#####grafici fago  2h######
A2h<-fago %>%
  filter(!is.na(CFUml)) %>% 
  filter(time=="2h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(cfu=log10(CFUml)) %>% 
  dabest(gruppo, cfu, 
         idx = list(c("Asaia pHM4","Asaia WSP","Asaia pHM4 + LPS")), 
         paired = FALSE
  )

p1<-plot(A2h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of phagocytosis activity Log10(CFU/ml)",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 10#,
         # rawplot.ylim = c(0,3.5),
         # effsize.ylim = c(-0.35,0.2)
)



####PLOT FAGO 2h#####
(p1)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of phagocytosis activity (PA) in Log10(UFC/ml) 
between different Treatment Group vs Asaia pHM4 as control group at 2h: 
         
1) The difference between the sample mean PA Log10(UFC/ml) in AsaiaWSP(5.92 UFC/ml) and AsaiapHM4 (5.90 UFC/ml) was 0.03 UFC/ml, with a
bootstraped 95% confidence interval from -0.45 to 0.46 UFC/ml; the bootstraped Welch two-sample t-test P value was 0.952 

2) The difference between the sample mean PA Log10(UFC/ml) in AsaiapHM4+LPS (5.86 UFC/ml) and AsaiapHM4 (5.90 UFC/ml) was 0.04 UFC/ml, with a
bootstraped 95% confidence interval from -0.53 to 0.32 UFC/ml; the bootstraped Welch two-sample t-test P value was 0.466.",
                     
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))






#####bootstrap fago 24h#####
x=fago %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Asaia pHM4","Asaia WSP")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(cfu=log10(CFUml)) %>% 
  drop_na(cfu) 

boot.t.test(cfu ~ gruppo, data=x)

xx=fago %>% 
  filter(gruppo %in% c("Asaia pHM4","Asaia pHM4 + LPS")) %>% 
  filter(time=="1h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(cfu=log10(CFUml)) %>% 
  drop_na(cfu) 

boot.t.test(cfu ~ gruppo, data=xx)












##### N. leishmanie#####




nleish<-nleish %>%
  select(gruppo,n.leish) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  dabest(gruppo, n.leish, 
         idx = list(c("Leishmania", "AsaiapHM4","AsaiaWSP", "Anfotericina" )), 
         paired = FALSE
  )

px<-plot(nleish,rawplot.type = "sinaplot", group.summaries = NULL,
     rawplot.ylabel = "number of Leishmania/macrophage",
     effsize.ylabel = "Unpaired mean difference",
     axes.title.fontsize = 10#,
     # rawplot.ylim = c(0,3.5),
     #effsize.ylim = c(-0.35,0.2)
)

(px)+plot_annotation(caption=
                       
"In this figure is reported the unpaired mean comparison of number of Leishmania/macrophages between different Treatment Group vs Leishmania as control group at 24h: 
         
1) The difference between the sample mean of parasites in AsaiapHM4 (1.12) and Med (1.44) was -0.32, with a
bootstraped 95% confidence interval from -0.52 to -0.13; the bootstraped Welch two-sample t-test P value was < 2.2e-16 

2) The difference between the sample mean of parasites in AsaiaWSP (0.63) and Med (1.44) was -0.81  , with a
bootstraped 95% confidence interval from -0.97 to -0.64; the bootstraped Welch two-sample t-test P value was 0.002. 

3) The difference between the sample mean of parasites in Leishmania (0.33) and Med (1.44) was -1.11, with a
bootstraped 95% confidence interval from -1.27 to -0.96; the bootstraped Welch two-sample t-test P value was < 2.2e-16. ",
                     theme=theme(plot.caption=element_text(hjust=0, size=12, margin=margin(t=0),family = "Comic Sans MS"))) &
  theme(plot.tag = element_text(size = 13))


x=nleish %>% 
  select(-replica) %>% 
  filter(gruppo %in% c("Leishmania","AsaiapHM4")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  # mutate(iL6=log10(IL6)) %>% 
  drop_na(n.leish) 

boot.t.test(n.leish ~ gruppo, data=x)

y=nleish %>% 
  select(-replica) %>% 
  filter(gruppo %in% c("Leishmania","AsaiaWSP")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  # mutate(iL6=log10(IL6)) %>% 
  drop_na(n.leish) 

boot.t.test(n.leish ~ gruppo, data=y)

z=nleish %>% 
  select(-replica) %>% 
  filter(gruppo %in% c("Leishmania","Anfotericina")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  # mutate(iL6=log10(IL6)) %>% 
  drop_na(n.leish) 

boot.t.test(n.leish ~ gruppo, data=z)

      



data<-list(exp=nleish)

























###Bayesian Negative binomial regression####

nleish<-nleish %>% 
  mutate(gruppo=factor(gruppo, levels=c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina")))


mod<-stan_glm(n.leish~gruppo, data=nleish, 
                family = neg_binomial_2,adapt_delta = .99)




draws <- as.data.frame(mod)
names(draws)<-c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina", "dispersione")

draws %>% 
  mutate("Leishmania"= exp(Leishmania),
         "AsaiaWSP"= exp(AsaiaWSP-Leishmania), 
         "AsaiapHM4"= exp(AsaiapHM4-Leishmania),
         "Anfotericina"= exp(Anfotericina-Leishmania)) %>% 
  pivot_longer(cols = c(1:4), names_to = "Group", values_to = "postBeta") %>% 
  mutate(Group=factor(Group, levels=c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina"))) %>% 
  ggplot(aes(x = postBeta, y = Group)) +
  geom_halfeyeh(color="navy",.width = c(0.95))+
  labs(x="Bayesian Posterior Distribution of number of Leishmania in
       Macrophages with uncertainity estimation (95% credibility intervals) of Bayesian Posterior Median)")+
  annotate(geom = "text", label="1.44(95%CI: 1.30-1.6)",
           x=1.44,
           y=0.8)+
  annotate(geom = "text", label="0.18(95%CI: 0.14-0.24)",
           x=0.18,
           y=1.8)+
  annotate(geom = "text", label="0.10(95%CI: 0.08-0.0.14)",
           x=0.10,
           y=2.8)+
  annotate(geom = "text", label="0.05(95%CI: 0.04-0.07)",
           x=0.05,
           y=3.8)

  


###Zero-inflated neg binom regression#####

nleish<-nleish %>% 
  mutate(gruppo=factor(gruppo, levels=c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina")))


 nleish %>% 
   ggplot(aes(n.leish))+geom_rug(aes(x = n.leish, y = 0), position = position_jitter(height = 0))+
   facet_wrap(~gruppo)+labs(x="Number of Leishmania in Macrophages", y="N.of Macrophages")+geom_histogram(bins=100)

 #### prop. Leishm##############
 
 x<-nleish %>% 
   mutate(prop=ifelse(n.leish==0,0,1)) 
 













######Citofluometria#######

names(cito)[3:10]<-c("CD40p","CD40gm","CD80p","CD80gm","CD86p","CD86gm","MHCp","MHCgm")

dabcito<-cito %>%
  dabest(gruppo, CD40gm, 
         idx = list(c("AsaiapHM4 +L", "Asaia WSP +L",  "Controllo positivo", "Leishmania", "Medium")),
         paired = FALSE
  )
plot(dabcito)

dabcito2<-cito %>%
  dabest(gruppo, CD80gm, 
         idx = list(c("AsaiapHM4 +L", "Asaia WSP +L",  "Controllo positivo", "Leishmania", "Medium")),
         paired = FALSE
  )
plot(dabcito2)


dabcito3<-cito %>%
  dabest(gruppo, CD86gm, 
         idx = list(c("AsaiapHM4 +L", "Asaia WSP +L",  "Controllo positivo", "Leishmania", "Medium")),
         paired = FALSE
  )
plot(dabcito3)


dabcito4<-cito %>%
  drop_na(MHCgm) %>% 
  dabest(gruppo, MHCgm, 
         idx = list(c("AsaiapHM4 +L", "Asaia WSP +L",  "Controllo positivo", "Leishmania", "Medium")),
         paired = FALSE
  )
plot(dabcito4)


#######################generic cumming plot#######
set.seed(54321)

N = 40
c1 <- rnorm(N, mean = 100, sd = 25)
c2 <- rnorm(N, mean = 100, sd = 50)
g1 <- rnorm(N, mean = 120, sd = 25)
g2 <- rnorm(N, mean = 80, sd = 50)
g3 <- rnorm(N, mean = 100, sd = 12)
g4 <- rnorm(N, mean = 100, sd = 50)
gender <- c(rep('Male', N/2), rep('Female', N/2))
id <- 1: N


wide.data <- 
  tibble::tibble(
    Control1 = c1, Control2 = c2,
    Group1 = g1, Group2 = g2, Group3 = g3, Group4 = g4,
    Gender = gender, ID = id)


my.data   <- 
  wide.data %>%
  tidyr::gather(key = Group, value = Measurement, -ID, -Gender)

library(dabestr)
## Loading required package: boot
## Loading required package: magrittr
two.group.unpaired <- 
  my.data %>%
  dabest(Group, Measurement, 
         # The idx below passes "Control" as the control group, 
         # and "Group1" as the test group. The mean difference
         # will be computed as mean(Group1) - mean(Control1).
         idx = c("Control1", "Group1"), 
         paired = FALSE)
p<-plot(two.group.unpaired, float.contrast = FALSE, group.summaries = NULL, 
        rawplot.type ="swarmplot")

  p+annotate("text", x = 2, y = 130, label = "Some text")
