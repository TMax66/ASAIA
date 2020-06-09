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
options(scipen = .999)
dati <- read_excel("Dati AsaiaWSP-Leishmania_MODIFICATO.xlsx")
nleish<-read_excel("nparinmacro.xlsx", sheet = "Foglio3")
rateinf<-read_excel("nminf.xlsx")
cito<-read_excel("cito.xlsx")
#Cytokine####

#IL6#######
IL6<-dati[,2:4]

######bootstrap t-test IL6  24h####  
x=IL6 %>% 
  filter(gruppo %in% c("AsaiapHM4","AsaiaWSP")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=x)

xx=IL6 %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=xx)

z=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z)

z2=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z2)

z3=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z3)

z4=IL6 %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z4)


#####grafici IL6 24h#####
il6A24h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, iL6, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(il6A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
            rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
            effsize.ylabel = "Unpaired mean difference",
            axes.title.fontsize = 12,
            rawplot.ylim = c(0,4.5),
            effsize.ylim = c(-0.1,0.4)
         )

il6B24h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, iL6, 
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(il6B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
            rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
            effsize.ylabel = "Unpaired mean difference",
            axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.5,0.8))

il6C<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, iL6, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(il6C,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)




###PLOT 1 IL6 24h####

(p1|p2)+plot_annotation(tag_levels = 'A', 
                        caption = "

In panel A is reported the unpaired mean comparison of Log10(IL6) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
between the sample mean Log10(IL6) in AsaiaWSP (3.75 pg/ml) and AsaiapHM4 (3.60 pg/ml) was 0.15 pg/ml, with a
bootstraped 95% confidence interval from 0.07 to 0.22 pg/ml; the bootstraped 
Welch two-sample t-test P value was 0.0008.

In panel B is reported the unpaired mean comparison of Log10(IL6) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
between the sample mean Log10(IL6) in AsaiaWSP L (3.36 pg/ml) and AsaiapHM4 (3.22 pg/ml) was 0.14 pg/ml, with a
bootstraped 95% confidence interval from -0.20 to 0.50 pg/ml; the bootstraped 
Welch two-sample t-test P value was 0.44. ",
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
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=x)

xx=IL6 %>% 
  #select(-time) %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=xx)

z=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z)

z2=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z2)

z3=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z3)

z4=IL6 %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  drop_na(IL6) 

boot.t.test(iL6 ~ gruppo, data=z4)

#####grafici IL6 48h#####
il6A48h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, iL6, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(il6A48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.1,0.9)
)

il6B48h<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, iL6, 
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(il6B48h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.5,0.95))

il6C<-IL6 %>%
  filter(!is.na(IL6)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  mutate(iL6=log10(IL6)) %>% 
  dabest(gruppo, iL6, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(il6C,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL6) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)

###PLOT 1 IL6 48h####

(p1|p2)+plot_annotation(tag_levels = 'A', 
                        caption = "

In panel A is reported the unpaired mean comparison of Log10(IL6) between Asaia WSP vs AsaiaphM4 as control group at 48h: The difference 
between the sample mean Log10(IL6) in AsaiaWSP (3.58 pg/ml) and AsaiapHM4 (3.20 pg/ml) was 0.38 pg/ml, with a
bootstraped 95% confidence interval from 0.19 to 0.58 pg/ml; the bootstraped 
Welch two-sample t-test P value was 0.0090.

In panel B is reported the unpaired mean comparison of Log10(IL6) Asaia WSP L vs AsaiaphM4 L as control group at 48h: The difference 
between the sample mean Log10(IL6) in AsaiaWSP L (3.37 pg/ml) and AsaiapHM4 (3.09 pg/ml) was 0.29 pg/ml, with a
bootstraped 95% confidence interval from -0.18 to 0.77 pg/ml; the bootstraped 
Welch two-sample t-test P value was 0.25. ",
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
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=x)

xx=TNF %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=xx)

z=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=z)

z2=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=z2)

z3=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo,data=z3)

z4=TNF %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=z4)



####grafici TNF 24h####
A24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, tnf, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.3,0.4)
)

B24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, tnf,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-1.2,1.2)
)



C24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, tnf, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)



#####PLOT 1 TNF 24h#####

(p1|p2)+plot_annotation(tag_levels = 'A', 
                        caption = "

In panel A is reported the unpaired mean comparison of Log10(TNF) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
between the sample mean Log10(TNF) in AsaiaWSP (2.95 pg/ml) and AsaiapHM4 (2.87 pg/ml) was 0.08 pg/ml, with a
bootstraped 95% confidence interval from -0.13 to 0.28 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.41.

In panel B is reported the unpaired mean comparison of Log10(TNF) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
between the sample mean Log10(TNF) in AsaiaWSP L (2.41 pg/ml) and AsaiapHM4 (2.29 pg/ml) was 0.11 pg/ml, with a
bootstraped 95% confidence interval from -0.60 to 0.87 pg/ml; the bootstraped Welch two-sample t-test P value was 0.76. ",
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
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=x)

xx=TNF %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=xx)

z=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=z)

z2=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=z2)

z3=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo,data=z3)

z4=TNF %>% 
  filter(time=="48h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(tnf=log10(TNF)) %>% 
  drop_na(tnf) 

boot.t.test(tnf ~ gruppo, data=z4)

#####grafici TNF 48h#####

A24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, tnf, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.3,0.4)
)

B24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, tnf,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-1.2,1.2)
)



C24h<-TNF %>%
  filter(!is.na(TNF)) %>% 
  filter(time=="48h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  mutate(tnf=log10(TNF)) %>% 
  dabest(gruppo, tnf, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(TNF) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)


#####Plot 1 TNF 48h#####
(p1|p2)+plot_annotation(tag_levels = 'A', 
                        caption = "

In panel A is reported the unpaired mean comparison of Log10(TNF) between Asaia WSP vs AsaiaphM4 as control group at 48h: The difference 
between the sample mean Log10(TNF) in AsaiaWSP (2.17 pg/ml) and AsaiapHM4 (2.13 pg/ml) was 0.05 pg/ml, with a
bootstraped 95% confidence interval from -0.14 to 0.22 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.68.

In panel B is reported the unpaired mean comparison of Log10(TNF) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
between the sample mean Log10(TNF) in AsaiaWSP L (2.23 pg/ml) and AsaiapHM4 (2.13 pg/ml) was 0.10 pg/ml, with a
bootstraped 95% confidence interval from -0.13 to 0.38 pg/ml; the bootstraped Welch two-sample t-test P value was 0.45. ",
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
  mutate(il1b=log10(IL1beta)) %>% 
  drop_na(il1b) 

boot.t.test(il1b ~ gruppo, data=x)

xx=IL1b %>% 
  filter(gruppo %in% c("Asaia pHM4 L","AsaiaWSP L")) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  drop_na(il1b) 

boot.t.test(il1b ~ gruppo, data=xx)

z=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiapHM4","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  drop_na(il1b) 

boot.t.test(il1b ~ gruppo,data=z)

z2=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("AsaiaWSP","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  drop_na(il1b) 

boot.t.test(il1b ~ gruppo, data=z2)

z3=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("Leishmania","Med")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  drop_na(il1b) 

boot.t.test(il1b ~ gruppo,data=z3)

z4=IL1b %>% 
  filter(time=="24h") %>% 
  filter(gruppo %in% c("LPS","Med")) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  drop_na(il1b) 

boot.t.test(il1b ~ gruppo, data=z4)


#####grafici IL1beta 24h######
A24h<-IL1b %>%
  filter(!is.na(IL1beta)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  dabest(gruppo, il1b, 
         idx = list(c("AsaiapHM4","AsaiaWSP" )), 
         paired = FALSE
  )

p1<-plot(A24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL1 beta) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.3,0.4)
)

B24h<-IL1b %>%
  filter(!is.na(IL1beta)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(il1b=log10(IL1beta)) %>% 
  dabest(gruppo, il1b,
         idx = list(c("Asaia pHM4 L","AsaiaWSP L" )), 
         paired = FALSE
  )

p2<-plot(B24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL1 beta) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12,
         rawplot.ylim = c(0,4.5),
         effsize.ylim = c(-0.1,0.1)
)



C24h<-IL1b %>%
  filter(!is.na(IL1b)) %>% 
  filter(time=="24h") %>% 
  mutate(gruppo=as.factor(gruppo)) %>%
  mutate(il1b=log10(IL1beta)) %>% 
  dabest(gruppo, il1b, 
         idx = list(c("Med","AsaiapHM4","AsaiaWSP",  "Leishmania", "LPS" )),
         paired = FALSE
  )

p3<-plot(C24h,rawplot.type = "sinaplot", group.summaries = NULL, float.contrast = FALSE,
         rawplot.ylabel = "Observed data of Log10(IL1 beta) pg/ml ",
         effsize.ylabel = "Unpaired mean difference",
         axes.title.fontsize = 12)





#####PLOT 1 IL1 beta 24h#####
(p1|p2)+plot_annotation(tag_levels = 'A', 
                        caption = "

In panel A is reported the unpaired mean comparison of Log10(IL1 beta) between Asaia WSP vs AsaiaphM4 as control group at 24h: The difference 
between the sample mean Log10(IL1 beta) in AsaiaWSP (2.55 pg/ml) and AsaiapHM4 (2.47 pg/ml) was 0.08 pg/ml, with a
bootstraped 95% confidence interval from -0.004 to 0.16 pg/ml; the bootstraped  Welch two-sample t-test P value was 0.080.

In panel B is reported the unpaired mean comparison of Log10(IL1 beta) Asaia WSP L vs AsaiaphM4 L as control group at 24h: The difference 
between the sample mean Log10(IL1 beta) in AsaiaWSP L (2.48 pg/ml) and AsaiapHM4 (2.47 pg/ml) was 0.18 pg/ml, with a
bootstraped 95% confidence interval from -0.02 to 0.06 pg/ml; the bootstraped Welch two-sample t-test P value was 0.373. ",
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


#########################################################################################################################
#########################################################################################################################


IL12p40<-dati[,c(2:3,7)]

IL12p40 %>% 
  filter(!is.na(IL12p40)) %>% 
  mutate(IL12p40=as.numeric(IL12p40)) %>% 
  mutate(log10IL12p40=log10(IL12p40))%>% 
  ggplot(aes(y=log10IL12p40, x=gruppo))+
  stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
  coord_flip()+labs(x="")+
  geom_jitter(width = 0.25, size=1.8)+
  facet_grid(~time)


dab<-IL12p40 %>%
  filter(!is.na(IL12p40)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10IL12p40=log10(IL12p40)) %>% 
  dabest(gruppo, log10IL12p40, 
         idx = list(c("AsaiapHM4","AsaiaWSP" ), 
                    c("Asaia pHM4 L","AsaiaWSP L" )),
         paired = FALSE
  )


plot(dab)


dab1<-IL12p40 %>%
  filter(!is.na(IL12p40)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10IL12p40=log10(IL12p40)) %>% 
  dabest(gruppo, log10IL12p40, 
         idx = list(c("AsaiaWSP","Leishmania", "LPS", "Med")),
         paired = FALSE
  )



plot(dab1)






ROS<-dati[,c(2:3,8)]

ROS %>% 
  filter(!is.na(ROS)) %>% 
  mutate(ROS=as.numeric(ROS)) %>% 
  mutate(log10ROS=log10(ROS))%>% 
  ggplot(aes(y=log10ROS, x=gruppo))+
  stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
  coord_flip()+labs(x="")+
  geom_jitter(width = 0.25)+
  facet_grid(~time)

dab<-ROS %>%
  filter(!is.na(ROS)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10ROS=log10(ROS)) %>% 
  dabest(gruppo, log10ROS, 
         idx = list(c("AsaiapHM4","AsaiaWSP" ), 
                    c("Asaia pHM4 L","AsaiaWSP L" )),
         paired = FALSE
  )


plot(dab)


dab1<-ROS %>%
  filter(!is.na(ROS)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10ROS=log10(ROS)) %>% 
  dabest(gruppo, log10ROS, 
         idx = list(c("AsaiaWSP","Leishmania", "LPS", "Med")),
         paired = FALSE
  )



plot(dab1)







Nitriti<-dati[,c(2:3,9)]

Nitriti %>% 
  filter(!is.na(Nitriti)) %>% 
  mutate(Nitriti=as.numeric(Nitriti)) %>% 
  mutate(log10Nitriti=log10(Nitriti))%>% 
  ggplot(aes(y=log10Nitriti, x=gruppo))+
  stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
  coord_flip()+labs(x="")+
  geom_jitter(width = 0.25)+
  facet_grid(~time)

dab<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10Nitriti=log10(Nitriti)) %>% 
  dabest(gruppo, log10Nitriti, 
         idx = list(c("AsaiapHM4","AsaiaWSP" ), 
                    c("Asaia pHM4 L","AsaiaWSP L" )),
         paired = FALSE
  )


plot(dab)


dab1<-Nitriti %>%
  filter(!is.na(Nitriti)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10Nitriti=log10(Nitriti)) %>% 
  dabest(gruppo, log10Nitriti, 
         idx = list(c("AsaiaWSP","Leishmania", "LPS", "Med")),
         paired = FALSE
  )



plot(dab1)






iNOS<-dati[,c(2:3,10)]

iNOS %>% 
  filter(!is.na(iNOS)) %>% 
  mutate(iNOS=as.numeric(iNOS)) %>% 
  mutate(log10iNOS=log10(iNOS))%>% 
  ggplot(aes(y=log10iNOS, x=gruppo))+
  stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
  coord_flip()+labs(x="")+
  geom_jitter(width = 0.25,size=1.8)+
  facet_grid(~time)


dab<-iNOS %>%
  filter(!is.na(iNOS)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10iNOS=log10(iNOS)) %>% 
  dabest(gruppo, log10iNOS, 
         idx = list(c("AsaiapHM4","AsaiaWSP" ), 
                    c("Asaia pHM4 L","AsaiaWSP L" )),
         paired = FALSE
  )


plot(dab)


dab1<-iNOS %>%
  filter(!is.na(iNOS)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10iNOS=log10(iNOS)) %>% 
  dabest(gruppo, log10iNOS, 
         idx = list(c("AsaiaWSP","Leishmania", "LPS", "Med")),
         paired = FALSE
  )



plot(dab1)







IL10<-dati[,c(2:3,11)]

IL10 %>% 
  filter(!is.na(IL10)) %>% 
  mutate(IL10=as.numeric(IL10)) %>% 
  mutate(log10IL10=log10(IL10))%>% 
  ggplot(aes(y=log10IL10, x=gruppo))+
  stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
  coord_flip()+labs(x="")+
  geom_jitter(width = 0.25,size=1.8)+
  facet_grid(~time)

dab<-IL10 %>%
  filter(!is.na(IL10)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10IL10=log10(IL10)) %>% 
  dabest(gruppo, log10IL10, 
         idx = list(c("AsaiapHM4","AsaiaWSP" ), 
                    c("Asaia pHM4 L","AsaiaWSP L" )),
         paired = FALSE
  )


plot(dab)


dab1<-IL10 %>%
  filter(!is.na(IL10)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10IL10=log10(IL10)) %>% 
  dabest(gruppo, log10IL10,
         idx = list(c("AsaiaWSP","Leishmania", "LPS", "Med")),
         paired = FALSE
  )



plot(dab1)


Arginasi<-dati[,c(2:3,12)]

Arginasi %>% 
    filter(!is.na(Arginasi)) %>% 
    mutate(Arginasi=as.numeric(Arginasi)) %>% 
    mutate(log10Arginasi=log10(Arginasi))%>% 
    ggplot(aes(y=log10Arginasi, x=gruppo))+
    stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
    coord_flip()+labs(x="")+
    geom_jitter(width = 0.25,size=1.8)+
    facet_grid(~time)

  
  
  
  
  

#####Fagocitosi####

fago <- read_excel("Dati fagocitosi.xlsx")
fago<-na.omit(fago)
names(fago)[4]<-"CFUml"

fago$time<-factor(fago$time, levels=c("1h", "2h", "24h"))

fago %>% 
    mutate(log10CFU=log10(CFUml)) %>% 
    ggplot(aes(x=gruppo, y=log10CFU))+
    geom_jitter(width = 0.25, size=1.8)+
    stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
    coord_flip()+labs(x="")+
    facet_wrap(~time)
    
  


dabf<-fago %>%
  filter(!is.na(CFUml)) %>% 
  filter(time!="24h") %>% 
  mutate(log10CFUml=log10(CFUml)) %>% 
  dabest(gruppo, log10CFUml, 
         idx = list(c("Asaia pHM4","Asaia WSP", "Asaia pHM4 + LPS")), 
                    
         paired = FALSE
  )


plot(dabf, color.column=time)
    #####Nparassiti####

nleish<-read_excel("~/Scrivania/pippo/nparinmacro.xlsx", sheet = "Foglio3")

nleish %>% 
  filter(!is.na(n.leish)) %>% 
  mutate(lognleish=log10(n.leish+1)) %>% 
  ggplot(aes(x=gruppo, y=lognleish))+
  geom_jitter(width = 0.25, size=1.8)+
  stat_summary(fun.data = "mean_cl_boot", colour = "blue")+
  coord_flip()+labs(x="")+
  facet_grid(~replica)

ggplot(nleish, aes(n.leish)) + geom_histogram() 



dab<-nleish %>%
  filter(!is.na(n.leish)) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  mutate(log10nleish=log10(n.leish+1)) %>% 
  dabest(gruppo, log10nleish, 
         idx = list(c("AsaiapHM4","AsaiaWSP" , "Leishmania", "Anfotericina")),
         paired = FALSE
  )





plot(dab, color.column = as.factor(replica))




##### N. leishmanie#####




nleish<-nleish %>%
  mutate(gruppo=as.factor(gruppo)) %>% 
  dabest(gruppo, n.leish, 
         idx = list(c("Leishmania", "AsaiapHM4","AsaiaWSP", "Anfotericina" )), 
         paired = FALSE
  )

plot(nleish,rawplot.type = "sinaplot")


x=nleish %>% 
  select(-replica) %>% 
  filter(gruppo %in% c("Leishmania","AsaiaWSP")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  # mutate(iL6=log10(IL6)) %>% 
  drop_na(n.leish) 

boot.t.test(n.leish ~ gruppo, data=x)

y=nleish %>% 
  select(-replica) %>% 
  filter(gruppo %in% c("Leishmania","AsaiapHM4")) %>% 
  mutate(gruppo=as.factor(gruppo)) %>% 
  # mutate(iL6=log10(IL6)) %>% 
  drop_na(n.leish) 

boot.t.test(n.leish ~ gruppo, data=y)

y=nleish %>% 
  select(-replica) %>% 
  filter(gruppo %in% c("Leishmania","AsaiapHM4")) %>% 
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
