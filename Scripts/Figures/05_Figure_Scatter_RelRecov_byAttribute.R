#-------------------#
# scatter plot of   #
# recovery 0 y vs   #
#  20 y             #
#  and lambda
#-------------------#

library(bayestestR)
library(tidyverse)
library(corrr)

load(file="Data/RecovRates.RData")
load(file="Data/ForestAttributes_lambda.RData")
Lambda_Attributes<-Attributes
rm(Attributes)

Recovery_Attributes<-recov_b %>% 
  mutate(attribute=ifelse(attribute=="N_vol","N",attribute),
         attribute=ifelse(attribute=="C_vol","C",attribute),
         attribute=ifelse(attribute=="soilbd","BD",attribute),
         attribute=ifelse(attribute=="sla","SLA",attribute),
         attribute=ifelse(attribute=="WD","WD",attribute),
         attribute=ifelse(attribute=="nfix","NF",attribute),
         attribute=ifelse(attribute=="gini","SH",attribute),
         attribute=ifelse(attribute=="maxDBH","DMAX",attribute),
         attribute=ifelse(attribute=="q0","SR",attribute),
         attribute=ifelse(attribute=="q2","SD",attribute),
         attribute=ifelse(attribute=="chao","SC",attribute),
         attribute=ifelse(attribute=="AGB","AGB",attribute)) %>% 
  dplyr::filter(.,!attribute=="P_vol")

Attributes_Diff20_0<-Recovery_Attributes %>% 
  select(., attribute, attribute_type,chronosequence, age,recov_kld ) %>% 
  mutate(attribute=factor(attribute),
         recov_kld=100*recov_kld,
         age= paste("y",age,sep="")) %>%
  reshape2::dcast(., attribute_type+ attribute + chronosequence ~ age, value.var="recov_kld", mean ) 
  

  RR_20y0y<-ggplot(Attributes_Diff20_0, aes(x=y0,y=y20, color=attribute_type,fill=attribute_type, group=attribute))+
          geom_abline(intercept=0, slope=1,color="black",lty=2)+
          geom_point(size=2,shape=21,alpha=0.5, color="transparent")+
          geom_smooth(method="gam",se=FALSE)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
  scale_fill_manual(name="Attribute group", 
                       labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                       values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    
  labs(y="Relative recovery after 20 y (%)",x="Relative recovery after 0 y (%)")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.17),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))
  
  
  RR_40y20y<-ggplot(Attributes_Diff20_0, aes(x=y20,y=y40, color=attribute_type,fill=attribute_type, group=attribute))+
    geom_abline(intercept=0, slope=1,color="black",lty=2)+
    geom_point(size=2,shape=21,alpha=0.5, color="transparent")+
    geom_smooth(method="gam",se=FALSE)+
    scale_color_manual(name="Attribute group", 
                       labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                       values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    scale_fill_manual(name="Attribute group", 
                      labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                      values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    
    labs(y="Relative recovery after 40 y (%)",x="Relative recovery after 20 y (%)")+
    scale_y_continuous(breaks=c(0,20,40,60,80,100))+
    theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                     axis.text.x = element_text(face="bold",size=8),
                     axis.text.y = element_text(face="bold",size=8),
                     legend.position = c(0.90,0.17),
                     legend.title = element_text(size=8,face="bold"),
                     legend.text = element_text(size=7),
                     panel.background = element_rect(fill="transparent", color=NA))
  
  RR_40y0y<-ggplot(Attributes_Diff20_0, aes(x=y0,y=y40, color=attribute_type,fill=attribute_type, group=attribute))+
    geom_abline(intercept=0, slope=1,color="black",lty=2)+
    geom_point(size=2,shape=21,alpha=0.5, color="transparent")+
    geom_smooth(method="gam",se=FALSE)+
    scale_color_manual(name="Attribute group", 
                       labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                       values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    scale_fill_manual(name="Attribute group", 
                      labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                      values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    
    labs(y="Relative recovery after 40 y (%)",x="Relative recovery after 0 y (%)")+
    scale_y_continuous(breaks=c(0,20,40,60,80,100))+
    theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                     axis.text.x = element_text(face="bold",size=8),
                     axis.text.y = element_text(face="bold",size=8),
                     legend.position = c(0.90,0.17),
                     legend.title = element_text(size=8,face="bold"),
                     legend.text = element_text(size=7),
                     panel.background = element_rect(fill="transparent", color=NA))
  
  
  scatter_ratesSI<-cowplot::plot_grid(RR_20y0y+theme( axis.text.x = element_text(size=12),
                                                    axis.text.y = element_text(size=12),
                                                    axis.title.y = element_text(face="bold",size=12),
                                                    axis.title.x = element_text(face="bold",size=12),
                                                    legend.position = "none"),
                                      RR_40y20y+theme(axis.text.y = element_text(size=12),
                                                  axis.text.x =element_text(size=12),
                                                  axis.title.y = element_text(face="bold",size=12),
                                                  axis.title.x = element_text(face="bold",size=12),
                                                  legend.position = "none"),
                                      # RR_40y0y+theme(legend.position = "none",
                                      #                  axis.text.x = element_text(size=12),
                                      #                  axis.text.y = element_text(face="bold",size=12),
                                      #                  axis.title.y = element_text(size=12),
                                      #                axis.title.x = element_text(face="bold",size=12)),
                                    labels=c("A","B","C"),label_size=12,
                                    ncol=1)
  
  legend<-cowplot::get_legend(RR_20y0y+theme( legend.direction = "horizontal",
                                              legend.text=element_text(face="bold",size=12),
                                              legend.position=c(0.5,0.5),
                                              legend.title=element_blank()))
  
  scatter_ratesSI<-cowplot::plot_grid(legend,scatter_ratesSI,  rel_heights = c(0.05,1),ncol=1)
  
  ggsave(scatter_ratesSI,filename = file.path("Results", "MultiD_RelRecoveryRates_scatter0y20y40y.png"), 
         width    = 25, 
         height   = 37, 
         units    = "cm")
  
  
  # 0 y recovery and lambda
  
  
  Recovery_Zero<-Recovery_Attributes %>% 
                 dplyr::filter(age==0) %>% 
                 group_by(attribute_type, attribute) %>% 
                 summarize(kld_med=median(recov_kld)) %>% 
                 mutate(kld_med=kld_med*100)
  
  Lambdaa<- Lambda_Attributes %>% 
    mutate(attribute=ifelse(attribute=="SoilBulkDensity","BD",attribute),
           attribute=ifelse(attribute=="soilN","N",attribute),
           attribute=ifelse(attribute=="soilC","C",attribute),
           attribute=ifelse(attribute=="sla","SLA",attribute),
           attribute=ifelse(attribute=="wd","WD",attribute),
           attribute=ifelse(attribute=="nfix","NF",attribute),
           attribute=ifelse(attribute=="Gini","SH",attribute),
           attribute=ifelse(attribute=="maxDBH","DMAX",attribute),
           attribute=ifelse(attribute=="q0","SR",attribute),
           attribute=ifelse(attribute=="q2","SD",attribute),
           attribute=ifelse(attribute=="chao","SC",attribute),
           attribute=ifelse(attribute=="AGB","AGB",attribute), 
           attribute_type=ifelse(attribute_type=="diversity","Diversity",attribute_type),
           attribute_type=ifelse(attribute_type=="structure","Structure",attribute_type),
           attribute_type=ifelse(attribute_type=="trait","Function",attribute_type),
           attribute_type=ifelse(attribute_type=="soil","Soil",attribute_type)) %>% 
           dplyr::filter(.,!attribute=="P_vol") %>% 
           group_by(attribute_type,attribute) %>% 
          summarize(med_lambda= median(value))
  
  
  Attributes_lambda_0<- left_join(Recovery_Zero,Lambdaa, by= c("attribute","attribute_type"))
  
  
  RR_lambda_0y<-ggplot(Attributes_lambda_0, aes(x=med_lambda,y=kld_med))+
    geom_smooth(color="black",method="lm",se=FALSE)+
    geom_point(aes(x=med_lambda,y=kld_med, fill=attribute_type),size=2,shape=21,color="transparent")+
    
    scale_color_manual(name="Attribute group", 
                       labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                       values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    scale_fill_manual(name="Attribute group", 
                      labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                      values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    
    labs(x=expression(bold("Recovery rate"~(lambda))),y="Relative recovery after 0 y (%)")+
    scale_y_continuous(breaks=c(0,20,40,60,80,100))+
    theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                     axis.text.x = element_text(face="bold",size=8),
                     axis.text.y = element_text(face="bold",size=8),
                     legend.position = c(0.90,0.17),
                     legend.title = element_text(size=8,face="bold"),
                     legend.text = element_text(size=7),
                     panel.background = element_rect(fill="transparent", color=NA))
  
  
  ggsave(RR_lambda_0y,filename = file.path("Results", "MultiD_Lambda_0yRecovery.png"), 
         width    = 15, 
         height   = 15, 
         units    = "cm")
  
  
  RR_lambda_0yNoLine<-ggplot(Attributes_lambda_0, aes(x=med_lambda,y=kld_med))+
    #geom_smooth(color="black",method="lm",se=FALSE)+
    geom_point(aes(x=med_lambda,y=kld_med, fill=attribute_type),size=2,shape=21,color="transparent")+
    
    scale_color_manual(name="Attribute group", 
                       labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                       values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    scale_fill_manual(name="Attribute group", 
                      labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                      values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
    
    labs(x=expression(bold("Recovery rate"~(lambda))),y="Relative recovery after 0 y (%)")+
    scale_y_continuous(breaks=c(0,20,40,60,80,100))+
    theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                     axis.text.x = element_text(face="bold",size=8),
                     axis.text.y = element_text(face="bold",size=8),
                     legend.position = c(0.90,0.17),
                     legend.title = element_text(size=8,face="bold"),
                     legend.text = element_text(size=7),
                     panel.background = element_rect(fill="transparent", color=NA))
  
  
  ggsave(RR_lambda_0yNoLine,filename = file.path("Results", "MultiD_Lambda_0yRecoveryNoLine.png"), 
         width    = 15, 
         height   = 15, 
         units    = "cm")
  


