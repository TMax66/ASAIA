###Zero-inflated neg binom regression#####
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
library(sjPlot)
library(ggeffects)
library(sjmisc)

nleish<-read_excel("nparinmacro.xlsx", sheet = "Foglio3")

nleish<-nleish %>% 
mutate("group"=gruppo) %>%
  select(group,n.leish)


nleish<-nleish %>% 
  mutate(group=factor(group, levels=c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina")))

#nleish %>% 
#ggplot(aes(n.leish))+geom_rug(aes(x = n.leish, y = 0), position = position_jitter(height = 0))+
#facet_wrap(~gruppo)+labs(x="Number of Leishmania in Macrophages", y="N.of Macrophages")+geom_histogram(bins=100)


mod <- zeroinfl(n.leish~ group,dist = "negbin",data = nleish)




m0 <- update(mod, . ~ 1)
# pchisq(2 * (logLik(mod) - logLik(m0)), df = 3, lower.tail=FALSE)

p<-plot_model(mod,grid = FALSE,vline.color = "blue", show.intercept = TRUE,
              show.values = TRUE, show.p = FALSE,transform = NULL)


p$conditional+theme_sjplot()+labs(title="Count model")

#p$zero_inflated+theme_sjplot()+labs(title="Zero-inflated model")



count<-p[["conditional"]][["data"]]


# logit2prob <- function(logit){
#   odds <- exp(logit)
#   prob <- odds / (1 + odds)
#   return(prob)
# }


g1<-count %>% 
  mutate(param = factor(c("Control Group","AsaiapHM4","AsaiaWSP","Amphotericin"), 
                        levels=c("Amphotericin","AsaiaWSP","AsaiapHM4","Control Group"))) %>% 
  ggplot(aes(x = estimate, y=param, label=p.label))+
  geom_point(aes(color=group), size=2.8)+ scale_color_manual(values=c("blue", "red"))+
  geom_segment(aes(x = conf.low, xend = conf.high, y=param, yend=param,
                   color=group))+theme_sjplot()+theme(legend.position = "none") +
  vline_0(color="grey3")+labs(title="Count Model", x="Log-Mean", y="")+
  geom_text(vjust=-1)


zero<-p[["zero_inflated"]][["data"]]

g2<-zero %>% 
  mutate(param = factor(c("Control Group","AsaiapHM4","AsaiaWSP","Amphotericin"), 
                        levels=c("Amphotericin","AsaiaWSP","AsaiapHM4","Control Group"))) %>%
  mutate(xx= c("yes", "no","no","no")) %>% 
  ggplot(aes(x = estimate, y=param,label=p.label))+
  geom_point(aes(color=xx), size=2.8)+ scale_color_manual(values=c("blue", "red"))+
  geom_segment(aes(x = conf.low, xend = conf.high, y=param, yend=param,
                   color=xx))+theme_sjplot()+theme(legend.position = "none")+scale_x_continuous( limits=c(-5, 5))+
  vline_0(color="grey3")+labs(title="Zero-inflated model", x="log-odds", y="")+
  geom_text(vjust=-1)


(g1|g2)+plot_annotation(tag_levels = 'A', 
              title = "Zero-inflated Negative Binomial Model of number of Leishmania in macrophages")



tab_model(mod,transform =NULL , string.ci = "95%CI", 
          file="zinb.html", show.loglik=TRUE, 
          show.aic = TRUE, show.r2 = FALSE)