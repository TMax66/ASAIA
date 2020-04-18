library(readxl)
library(tidyverse)
library(betareg)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

###dati#####
dt<- read_excel("Dati citofluorimetro AGGIORNATO 15.04.xlsx",sheet = "MHC")

dt<-dt[1:18,]
dt<-dt %>% 
select(-1,-5,-8) 

####MHC####
MHC<-dt %>% 
 select(plate=Piastra, well=Pozzetto, group=gruppo, "cellMHC+" =`N° cell MHC +`, totcell=`N° cell tot`) %>% 
 drop_na() %>% 
  mutate(group=factor(group, 
  levels=c("Controllo positivo","Medium","Leishmania","AsaiapHM4 +L","Asaia WSP +L"))) 
###multilevel poisson ###
MHC %>% 
  group_by(group) %>% 
  summarise(prop=mean(`cellMHC+`/totcell))



mhc<-glmer(`cellMHC+`~group+offset(log(totcell))+(1|plate/well),family=poisson , data=MHC)


#tab_model(mhc, show.icc = FALSE, show.reflvl = TRUE, transform = NULL, show)

####costimolatorie##############

dt<- read_excel("Dati citofluorimetro AGGIORNATO 15.04.xlsx",sheet = "CD")
dt<-dt[1:24,]

####CD40#####
CD40<-dt %>% 
  select(plate=Piastra, well=Pozzetto, group=gruppo, "cellCD40+" =`N°cellule CD40+`, totcell=`N° cell tot`) %>% 
  drop_na() %>% 
  mutate(group=factor(group, 
                      levels=c("Controllo positivo","Medium","Leishmania","AsaiapHM4 +L","Asaia WSP +L"))) 

CD40 %>% 
  group_by(group) %>% 
  summarise(prop=mean(`cellCD40+`/totcell))

cd40<-glmer(`cellCD40+`~group+offset(log(totcell))+(1|plate/well),family=poisson , data=CD40)

####CD80##############
CD80<-dt %>% 
  select(plate=Piastra, well=Pozzetto, group=gruppo, "cellCD80+" =`N° cellule CD80+`, totcell=`N° cellule tot...10`) %>% 
  drop_na() %>% 
  mutate(group=factor(group, 
                      levels=c("Controllo positivo","Medium","Leishmania","AsaiapHM4 +L","Asaia WSP +L"))) 

cd80<-glmer(`cellCD80+`~group+offset(log(totcell))+(1|plate/well),family=poisson , data=CD80)

#####CD86####
CD86<-dt %>% 
  select(plate=Piastra, well=Pozzetto, group=gruppo, "cellCD86+" =`N° cell CD86+`, totcell=`N° cellule tot...14`) %>% 
  drop_na() %>% 
  mutate(group=factor(group, 
                      levels=c("Controllo positivo","Medium","Leishmania","AsaiapHM4 +L","Asaia WSP +L"))) 

cd86<-glmer(`cellCD86+`~group+offset(log(totcell))+(1|plate/well),family=poisson , data=CD86)



##########################################################################################
##########################################################################################
tab_model(mhc, cd40, cd80,cd86,transform =NULL, file="cell.html", show.reflvl=TRUE, 
          show.ci = FALSE, show.re.var = FALSE)

