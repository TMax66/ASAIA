###Zero-inflated neg binom regression#####
    library(readxl)
    library(tidyverse)
    library(dabestr)
    library(rstanarm)
    library(bayesplot)
    library(tidybayes)
   # library(BayesPostEst)
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

# x<-nleish %>% 
#   mutate(Zero=ifelse(n.leish==0, 0, 1))

# %>% 
#   ggplot(aes(x=n.leish))+geom_bar()
#  
#   
  
  


nleish<-nleish %>% 
  mutate(group=factor(group, levels=c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina")))

#nleish %>% 
#ggplot(aes(n.leish))+geom_rug(aes(x = n.leish, y = 0), position = position_jitter(height = 0))+
#facet_wrap(~gruppo)+labs(x="Number of Leishmania in Macrophages", y="N.of Macrophages")+geom_histogram(bins=100)


mod <- zeroinfl(n.leish~ group,dist = "negbin",data = nleish)
p<-plot_model(mod,grid = FALSE,vline.color = "blue", show.intercept = TRUE,
              show.values = TRUE, show.p = FALSE,transform = NULL)
count<-p[["conditional"]][["data"]]

lab<-count[,10]

L<-rnorm(1000000,0.4791586, 0.10890152)
AsaiapHM4<-rnorm(1000000, -0.1052806, 0.09441159)
AsaiaWSP<-rnorm(1000000, -0.8806682, 0.11154448)
Anfo<-rnorm(100000, -1.2729636, 0.15977129)

dt<-data.frame(L,AsaiapHM4, AsaiaWSP,Anfo)



dt %>% 
  pivot_longer(1:4,names_to = "group", values_to = "estimate") %>% 
  mutate(group=factor(group, levels=c("Anfo","AsaiaWSP","AsaiapHM4","L"))) %>% 
  mutate(group=recode(group,"L"="Leishmania (L)",
                      "AsaiapHM4"="Asaia PHM4 + L",
                      "AsaiaWSP"="Asaia WSP + L",
                      "Anfo"= "L+ Amphotericina B")) %>% 
  ggplot(aes(x=estimate, y=group, fill=group))+
  geom_halfeyeh(alpha=0.8, .width = c(0.95))+
  scale_fill_manual(values=c("steelblue3","steelblue2", "steelblue1","gray" ))+
  vline_0(color="red",linetype = 2)+labs(title="Count Model", x="Log-Mean", y="")+
  theme_ggeffects()+ 
  theme(legend.position = 'none')+xlim(-2,1.5)+
  annotate(geom = "text", label="0.48(95%CI:0.26 0.69)",
           x=0.48,
           y=3.8,size=3.5)+
    annotate(geom = "text", label="-0.10(95%CI:-0.29  0.79)",
             x=-0.10,
             y=2.8,size=3.5)+
    annotate(geom = "text", label="-0.88(95%CI:-1.09 -0.66)",
             x=-0.88,
             y=1.8,size=3.5)+
    annotate(geom = "text", label="-1.27(95%CI:-1.59 -0.96)",
             x=-1.27,
             y=0.8,size=3.5)


# +
#   annotate(geom = "text", x = 1.1, y = 4.5, label = "Uncertainty of
#            estimated log mean of parasites count in Lieshmania (L) group with 95%CI", 
#            size=3.5, color="blue",  hjust=0)

















#plot_summs(mod, scale = TRUE, plot.distributions = TRUE, inner_ci_level = .9)
# m0 <- update(mod, . ~ 1)
# # pchisq(2 * (logLik(mod) - logLik(m0)), df = 3, lower.tail=FALSE)
# 
#p<-plot_model(mod,grid = FALSE,vline.color = "blue", show.intercept = TRUE,
         #   show.values = TRUE, show.p = FALSE,transform = NULL)
# 
# 
# p$conditional+theme_sjplot()+labs(title="Count model")
# 
# #p$zero_inflated+theme_sjplot()+labs(title="Zero-inflated model")
# 
#count<-p[["conditional"]][["data"]]


  




  # logit2prob <- function(logit){
#   odds <- exp(logit)
#   prob <- odds / (1 + odds)
#   return(prob)
# }


g1<-count %>% 
mutate(param = factor(c("Leishmania(L)","L+AsaiapHM4","L+AsaiaWSP","L+Amphotericin B"), 
                        levels=c("L+Amphotericin B","L+AsaiaWSP","L+AsaiapHM4","Leishmania(L)"))) 


 
#   ggplot(aes(x = estimate, y=param, label=p.label))+
#   geom_point(aes(color=group), size=2.8)+ scale_color_manual(values=c("blue", "red"))+
#   geom_segment(aes(x = conf.low, xend = conf.high, y=param, yend=param,
#                    color=group))+theme_sjplot()+theme(legend.position = "none") +
#   vline_0(color="grey3")+labs(title="Count Model", x="Log-Mean", y="")+
#   stat_function(fun=dnorm)
#   geom_text(vjust=-1)



  
 
  
  
  
  
  
  
  
  

  
  
  
  
  
zero<-p[["zero_inflated"]][["data"]]

g2<-zero %>% 
  mutate(param = factor(c("Leishmania(L)","AsaiapHM4","AsaiaWSP","L+Amphotericin B"), 
                        levels=c("L+Amphotericin B","AsaiaWSP","AsaiapHM4","Leishmania(L)"))) %>%
  mutate(xx= c("yes", "no","no","no")) %>% 
  ggplot(aes(x = estimate, y=param,label=p.label))+
  geom_point(aes(color=xx), size=2.8)+ scale_color_manual(values=c("blue", "red"))+
  geom_segment(aes(x = conf.low, xend = conf.high, y=param, yend=param,
                   color=xx))+theme_sjplot()+theme(legend.position = "none")+scale_x_continuous( limits=c(-5, 5))+
  vline_0(color="grey3")+labs(title="Zero-inflated Model", x="log-Odds", y="")+
  geom_text(vjust=-1)


(g1|g2)+plot_annotation(tag_levels = 'A', 
              title = "Zero-inflated Negative Binomial Model of number of Leishmania in macrophages")



tab_model(mod,transform =NULL , string.ci = "95%CI", 
          file="zinb.html", show.loglik=TRUE, 
          show.aic = TRUE, show.r2 = FALSE)
