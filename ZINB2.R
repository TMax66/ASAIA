###Zero-inflated neg binom regression#####
    library(readxl)
    library(tidyverse)
    library(dabestr)
    library(rstanarm)
    library(bayesplot)
    library(tidybayes)
    #library(BayesPostEst)
    library(patchwork)
    library(gridExtra)
    library(MKinfer)
    library(ggpubr)
    library(pscl)
    library(sjPlot)
    library(ggeffects)
    library(sjmisc)
    
    nleish<-read_excel("nparinmacro2.xlsx")
    
    
    nleish<-nleish %>% 
      select(group,nleish)
   
     nleish<-nleish %>% 
      mutate(group=factor(group, levels=c("Leishmania","AsaiapHM4","AsaiaWSP","Anfotericina")))
    
    mod <- zeroinfl(nleish~ group,dist = "negbin",data = nleish)
    p<-plot_model(mod,grid = FALSE,vline.color = "blue", show.intercept = TRUE,
                  show.values = TRUE, show.p = FALSE,transform = NULL)
    count<-p[["conditional"]][["data"]]
 
    
    L<-rnorm(1000000,0.66432399, 0.07214864)
    AsaiapHM4<-rnorm(1000000, -0.04252611, 0.08137206)
    AsaiaWSP<-rnorm(1000000, -0.14455194 , 0.08268058)
    
    dt<-data.frame(L,AsaiapHM4, AsaiaWSP)
    
    
    
    dt %>% 
      pivot_longer(1:3,names_to = "group", values_to = "estimate") %>% 
      mutate(group=factor(group, levels=c("AsaiaWSP","AsaiapHM4","L"))) %>% 
      mutate(group=recode(group,"L"="Leishmania (L)",
                          "AsaiapHM4"="Asaia PHM4 + L",
                          "AsaiaWSP"="Asaia WSP + L"
                           )) %>% 
      ggplot(aes(x=estimate, y=group, fill=group))+
      geom_halfeyeh(alpha=0.8, .width = c(0.95))+
      scale_fill_manual(values=c("firebrick4","firebrick4","gray" ))+
      vline_0(color="red",linetype = 2)+labs(title="Coefficients plot", x="Count model estimated coefficients (log-mean)", y="")+
      theme_ggeffects()+ 
      theme(legend.position = 'none')+xlim(-1,1)+
      annotate(geom = "text", label="0.66(95%CI:0.52 0.80)",
               x=0.66,
               y=2.8,size=3.5)+
      annotate(geom = "text", label="-0.04(95%CI:-0.20  0.11)",
               x=-0.04,
               y=1.8,size=3.5)+
      annotate(geom = "text", label="-0.14(95%CI:-0.30 -0.02)",
               x=-0.14,
               y=0.8,size=3.5)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # 
# leish<-nleish %>% 
#   select(group,nleish)
# 
# 
# leish<-leish %>% 
#   mutate(group=factor(group, levels=c("Leishmania","AsaiapHM4","AsaiaWSP")))
# 
# #nleish %>% 
# #ggplot(aes(n.leish))+geom_rug(aes(x = n.leish, y = 0), position = position_jitter(height = 0))+
# #facet_wrap(~gruppo)+labs(x="Number of Leishmania in Macrophages", y="N.of Macrophages")+geom_histogram(bins=100)
# 
# 
# mod <- zeroinfl(nleish~ group,dist = "negbin",data = leish)
# 
# 
# 
# 
# # m0 <- update(mod, . ~ 1)
# # # pchisq(2 * (logLik(mod) - logLik(m0)), df = 3, lower.tail=FALSE)
# # 
# p<-plot_model(mod,grid = FALSE,vline.color = "blue", show.intercept = TRUE,
#              show.values = TRUE, show.p = FALSE,transform = NULL)
# # 
# # 
# #p$conditional+theme_sjplot()+labs(title="Count model")
# # 
# #p$zero_inflated+theme_sjplot()+labs(title="Zero-inflated model")
# # 
# 
# 
# count<-p[["conditional"]][["data"]]
# 
# 
# # logit2prob <- function(logit){
# #   odds <- exp(logit)
# #   prob <- odds / (1 + odds)
# #   return(prob)
# # }
# 
# 
# g1<-count %>% 
#   mutate(param = factor(c("Leishmania(L)","AsaiapHM4","AsaiaWSP"), 
#                         levels=c("AsaiaWSP","AsaiapHM4","Leishmania(L)"))) %>% 
#   ggplot(aes(x = estimate, y=param, label=p.label))+
#   geom_point(aes(color=group), size=2.8)+ scale_color_manual(values=c("blue", "red"))+
#   geom_segment(aes(x = conf.low, xend = conf.high, y=param, yend=param,
#                    color=group))+theme_sjplot()+theme(legend.position = "none") +
#   vline_0(color="grey3")+labs(title="Count Model", x="Log-Mean", y="")+
#   geom_text(vjust=-1)
# 
# 
# zero<-p[["zero_inflated"]][["data"]]
# 
# g2<-zero %>% 
#   mutate(param = factor(c("Leishmania(L)","AsaiapHM4","AsaiaWSP"), 
#                         levels=c("AsaiaWSP","AsaiapHM4","Leishmania(L)"))) %>% 
#   mutate(xx= c("yes", "no","no")) %>% 
#   ggplot(aes(x = estimate, y=param,label=p.label))+
#   geom_point( aes(color=group), size=2.8)+ scale_color_manual(values=c("red", "blue"))+
#   geom_segment(aes(x = conf.low, xend = conf.high, y=param, yend=param,
#                    color=group))+theme_sjplot()+theme(legend.position = "none")+scale_x_continuous( limits=c(-5, 5))+
#   vline_0(color="grey3")+labs(title="Zero-inflated Model", x="log-Odds", y="")+
#   geom_text(vjust=-1)
# 
# 
# (g1|g2)+plot_annotation(tag_levels = 'A', 
#               title = "Zero-inflated Negative Binomial Model of number of Leishmania in macrophages")
# 
# 
# 
# tab_model(mod,transform =NULL , string.ci = "95%CI", 
#           file="zinb2.html", show.loglik=TRUE, 
#           show.aic = TRUE, show.r2 = FALSE)
# 
# 
# 
# 
