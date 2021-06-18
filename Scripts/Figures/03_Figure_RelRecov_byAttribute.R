#----------------#
# Summarize     #
# recovery rates #
# and lambda     #
#----------------#
# Main text: Recovery after 20
# SI: Recovery after 40

library(bayestestR)
library(tidyverse)

# Files too big for Github
# load(file="Data/RecovCurvEst.RData")
# recov_curve<-recov_b
# rm(recov_b)
# 
# load(file="Data/RecovRates.RData")
# 
# load(file="Data/ForestAttributes_lambda.RData")
# Lambda_Attributes<-Attributes
# rm(Attributes)

#------------------#
# Recovery Rates   #
#------------------#

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

#-------------------------#
# Recov. rates / 0 y      #
#-------------------------#

Attributes_0y<-Recovery_Attributes %>% 
  dplyr::filter(.,age=="0") %>%
  mutate(attribute=factor(attribute),
         recov_kld=100*recov_kld) %>% 
  group_by(attribute_type,attribute) %>% 
  summarize(median=median(recov_kld),
            lowerCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(median=round(median, digits=2),
         lowerCredInt=round(lowerCredInt,digits=2),
         upperCredInt=round(upperCredInt,digits=2))

Attributes_0y$attribute_type<-factor(Attributes_0y$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

RR_0yp<-Attributes_0y %>% 
  ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
  geom_hline(yintercept=90,color="grey70",lty=2)+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
  annotate("text", label = "After 0 y", x = 2, y = 110, size = 4, colour = "black")+
  
  labs(y="Relative recovery (%)",x="")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  coord_cartesian(ylim = c(0, 110))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

#-------------------------#
# Recov. Rates / 20 y     #
#-------------------------#

Attributes_20y<-Recovery_Attributes %>% 
                 dplyr::filter(.,age=="20") %>%
                 mutate(attribute=factor(attribute),
                        recov_kld=100*recov_kld) %>% 
                group_by(attribute_type,attribute) %>% 
                summarize(median=median(recov_kld),
                       lowerCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_low,
                       upperCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_high) %>% 
               mutate(median=round(median, digits=2),
                      lowerCredInt=round(lowerCredInt,digits=2),
                      upperCredInt=round(upperCredInt,digits=2))

Attributes_20y$attribute_type<-factor(Attributes_20y$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

RR_20yp<-Attributes_20y %>% 
         ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
                geom_hline(yintercept=90,color="grey70",lty=2)+
                geom_point(size=2)+
                geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
                scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
         annotate("text", label = "After 20 y", x = 2, y = 110, size = 4, colour = "black")+

         labs(y="Relative recovery (%)",x="")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  coord_cartesian(ylim = c(0, 110))+
         theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                          axis.text.x = element_text(face="bold",size=8),
                          axis.text.y = element_text(face="bold",size=8),
                          legend.position = c(0.90,0.87),
                          legend.title = element_text(size=8,face="bold"),
                          legend.text = element_text(size=7),
                          panel.background = element_rect(fill="transparent", color=NA))

#-----------------------------------#
# Changes in attributes: 20 y - 0y  #
#-----------------------------------#

Attributes_Diff20_0<-Recovery_Attributes %>% 
  select(., attribute, attribute_type,chronosequence, age,recov_kld ) %>% 
  mutate(attribute=factor(attribute),
         recov_kld=100*recov_kld,
         age= paste("y",age,sep="")) %>%
  reshape2::dcast(., attribute_type+ attribute + chronosequence ~ age, value.var="recov_kld", mean ) %>%  
  mutate(recov_diff = y20 - y0) %>% 
  group_by(attribute_type,attribute) %>% 
  summarize(median=median(recov_diff),
            lowerCredInt=ci(recov_diff, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(recov_diff, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(median=round(median, digits=2),
         lowerCredInt=round(lowerCredInt,digits=2),
         upperCredInt=round(upperCredInt,digits=2))

Attributes_Diff20_0$attribute_type<-factor(Attributes_Diff20_0$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

RR_changep<-Attributes_Diff20_0 %>% 
  ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
  #geom_hline(yintercept=0,color="grey70",lty=2)+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
  #annotate("text", label = "After 40 y", x = 2, y = 110, size = 4, colour = "black")+
  
  labs(y="Change in relative recovery (20 y - 0 y, %)",x="")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  coord_cartesian(ylim = c(-5, 110))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

#--------------------------#
# time to recovery 90%     #
#--------------------------#

recov_curve<-recov_curve %>% 
             unite("seq_attr", c("chronosequence","attribute"),sep=":",remove=F)

recov_curvee<-recov_curve %>% 
  group_by(attribute_type, attribute, age) %>% 
  summarize(recov_kldd=median(recov_kld),
            lowerCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_high) 

curves<-ggplot(recov_curvee, aes(x=age, y=recov_kldd,colour=attribute))+
  geom_vline(xintercept=20,color="grey70",lty=2)+
  geom_vline(xintercept=40,color="grey70",lty=2)+
  geom_hline(yintercept=90,color="grey70",lty=2)+
  geom_line(aes(x=age, y=recov_kldd,group=attribute))
  #geom_line(aes(x=age, y=lowerCredInt,group=attribute))
  #geom_line(aes(x=age, y=upperCredInt,group=attribute))
  

# median recovery time
medRecovTime<-ggplot_build(curves)$data[[4]] 

medRecovTime<-medRecovTime %>%
  mutate(y=round(y,digits=0)) %>% 
  dplyr::filter(y==90) %>% 
  group_by(group) %>% 
  summarize(med_age=median(x)) %>% 
  mutate(attribute=levels(recov_curvee$attribute)) %>% 
  arrange(.,-med_age) %>% 
medRecovTime$attribute<-levels(recov_curvee$attribute)

# confidence intervals

RecovTimeCI<-recov_curve %>% 
  mutate(recov_kld=round(recov_kld,digits=0)) %>% 
  group_by(attribute) %>% 
  dplyr::filter(recov_kld==90) %>% 
  group_by(attribute_type, attribute) %>% 
  summarize(recov_kldd=median(age),
            lowerRecov=ci(age, method="HDI", ci=0.95)$CI_low,
            upperRecov=ci(age, method="HDI", ci=0.95)$CI_high)%>% 
   mutate(recov_kldd=round(recov_kldd, digits=2),
          lowerRecov=round(lowerRecov,digits=2),
          upperRecov=round(upperRecov,digits=2))

RecovTimeTog<-medRecovTime %>% 
              left_join(., RecovTimeCI, by="attribute") %>%  
              select(., attribute, attribute_type, med_age, recov_kldd_age=recov_kldd, lowerRecov,upperRecov) %>% 
              mutate(upperRecov=ifelse(med_age>upperRecov,med_age,upperRecov))

RecovTimeTog$attribute_type<-factor(RecovTimeTog$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

rec_time<-RecovTimeTog %>% 
  ggplot(aes(x=reorder(attribute,med_age),y=med_age,color=attribute_type, group=attribute_type))+
  geom_point(size=4)+
  geom_pointrange(aes(ymin=lowerRecov,ymax=upperRecov,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+
  
  labs(y="Time to 90% recovery (y)",x="")+
  scale_y_continuous(breaks=c(0,25,50,75,100,120))+
  coord_cartesian(ylim = c(0, 125))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.10,0.87),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

# together

recov_rates<-cowplot::plot_grid(RR_0yp+theme( axis.text.x = element_text(size=12),
                                                axis.text.y = element_text(size=12),
                                                axis.title.y = element_text(face="bold",size=12),
                                                legend.position = "none"),
                                 RR_20yp+theme(axis.text.y = element_text(size=12),
                                                axis.text.x =element_text(size=12),
                                                axis.title.y = element_text(face="bold",size=12),
                                                legend.position = "none"),
                                 RR_changep+theme(legend.position = "none",
                                                   axis.text.x = element_text(size=12),
                                                   axis.text.y = element_text(face="bold",size=12),
                                                   axis.title.y = element_text(size=12)),
                                 rec_time+theme(legend.position = "none",
                                                        axis.title.y = element_text(face="bold",size=12),
                                                        axis.text.y = element_text(face="bold",size=12),
                                                        axis.text.x = element_text(size=12)),
                                  labels=c("A","B","C","D"),label_size=12,
                                  ncol=2)

legend<-cowplot::get_legend(RR_20yp+theme( legend.direction = "horizontal",
                                      legend.text=element_text(face="bold",size=12),
                                      legend.position=c(0.5,0.5),
                                      legend.title=element_blank()))

recov_rates<-cowplot::plot_grid(legend,recov_rates,  rel_heights = c(0.1,1),ncol=1)

ggsave(recov_rates,filename = file.path("Results", "MultiD_RelRecoveryRates_20y.png"), 
       width    = 37, 
       height   = 27, 
       units    = "cm")

#---------------------#
# 40 y recovery       #
#---------------------#

#-------------------------#
# Recov. rates / 20 y      #
#-------------------------#

Attributes_20yy<-Recovery_Attributes %>% 
  dplyr::filter(.,age=="20") %>%
  mutate(attribute=factor(attribute),
         recov_kld=100*recov_kld) %>% 
  group_by(attribute_type,attribute) %>% 
  summarize(median=median(recov_kld),
            lowerCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(median=round(median, digits=2),
         lowerCredInt=round(lowerCredInt,digits=2),
         upperCredInt=round(upperCredInt,digits=2))

Attributes_20yy$attribute_type<-factor(Attributes_20yy$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

RR_20ypp<-Attributes_20yy %>% 
  ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
  geom_hline(yintercept=90,color="grey70",lty=2)+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
  annotate("text", label = "After 20 y", x = 2, y = 110, size = 4, colour = "black")+
  
  labs(y="Relative recovery (%)",x="")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  coord_cartesian(ylim = c(0, 110))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

#-------------------------#
# Recov. Rates / 40 y     #
#-------------------------#

Attributes_40yy<-Recovery_Attributes %>% 
  dplyr::filter(.,age=="40") %>%
  mutate(attribute=factor(attribute),
         recov_kld=100*recov_kld) %>% 
  group_by(attribute_type,attribute) %>% 
  summarize(median=median(recov_kld),
            lowerCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(recov_kld, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(median=round(median, digits=2),
         lowerCredInt=round(lowerCredInt,digits=2),
         upperCredInt=round(upperCredInt,digits=2))

Attributes_40yy$attribute_type<-factor(Attributes_40yy$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

RR_40yp<-Attributes_40yy %>% 
  ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
  geom_hline(yintercept=90,color="grey70",lty=2)+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
  annotate("text", label = "After 40 y", x = 2, y = 110, size = 4, colour = "black")+
  
  labs(y="Relative recovery (%)",x="")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  coord_cartesian(ylim = c(0, 110))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

#-----------------------------------#
# Changes in attributes: 20 y - 0y  #
#-----------------------------------#

Attributes_Diff40_20<-Recovery_Attributes %>% 
  select(., attribute, attribute_type,chronosequence, age,recov_kld ) %>% 
  mutate(attribute=factor(attribute),
         recov_kld=100*recov_kld,
         age= paste("y",age,sep="")) %>%
  reshape2::dcast(., attribute_type+ attribute + chronosequence ~ age, value.var="recov_kld", mean ) %>%  
  mutate(recov_diff = y40 - y20) %>% 
  group_by(attribute_type,attribute) %>% 
  summarize(median=median(recov_diff),
            lowerCredInt=ci(recov_diff, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(recov_diff, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(median=round(median, digits=2),
         lowerCredInt=round(lowerCredInt,digits=2),
         upperCredInt=round(upperCredInt,digits=2))

Attributes_Diff40_20$attribute_type<-factor(Attributes_Diff40_20$attribute_type, levels=c("Diversity","Function","Soil","Structure"))

RR_changeSIp<-Attributes_Diff40_20 %>% 
  ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
  geom_hline(yintercept=0,color="grey70",lty=2)+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+ 
  #annotate("text", label = "After 40 y", x = 2, y = 110, size = 4, colour = "black")+
  
  labs(y="Change in relative recovery (40 y - 20 y, %)",x="")+
  scale_y_continuous(breaks=c(-5,0,10,20,30,40))+
  coord_cartesian(ylim = c(-6.5, 35))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

# together

recov_ratesSI<-cowplot::plot_grid(RR_20ypp+theme( axis.text.x = element_text(size=12),
                                              axis.text.y = element_text(size=12),
                                              axis.title.y = element_text(face="bold",size=12),
                                              legend.position = "none"),
                                RR_40yp+theme(axis.text.y = element_text(size=12),
                                              axis.text.x =element_text(size=12),
                                              axis.title.y = element_text(face="bold",size=12),
                                              legend.position = "none"),
                                RR_changeSIp+theme(legend.position = "none",
                                                 axis.text.x = element_text(size=12),
                                                 axis.text.y = element_text(face="bold",size=12),
                                                 axis.title.y = element_text(size=12)),
                                labels=c("A","B","C"),label_size=12,
                                ncol=2)

legend<-cowplot::get_legend(RR_20ypp+theme( legend.direction = "horizontal",
                                           legend.text=element_text(face="bold",size=12),
                                           legend.position=c(0.5,0.5),
                                           legend.title=element_blank()))

recov_ratesSI<-cowplot::plot_grid(legend,recov_ratesSI,  rel_heights = c(0.1,1),ncol=1)

ggsave(recov_ratesSI,filename = file.path("Results", "MultiD_RelRecoveryRates_40y.png"), 
       width    = 37, 
       height   = 27, 
       units    = "cm")

#---------------------#
# lambda               #
#---------------------#

Lambda_Attributes<-Lambda_Attributes %>% 
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
         attribute=ifelse(attribute=="AGB","AGB",attribute)) 

Lambda_Attributess<-Lambda_Attributes %>% 
  group_by(attribute_type,attribute) %>% 
  summarize(median=median(value),
            lowerCredInt=ci(value, method="HDI", ci=0.95)$CI_low,
            upperCredInt=ci(value, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(median=round(median, digits=2),
         lowerCredInt=round(lowerCredInt,digits=2),
         upperCredInt=round(upperCredInt,digits=2)) %>% 
  arrange(.,-median)

Lambda_Attributes$attribute_type<-factor(Lambda_Attributes$attribute_type, levels=c("diversity","trait","soil","structure"))

lambda_rec<-Lambda_Attributess %>% 
  ggplot(aes(x=reorder(attribute,-median),y=median,color=attribute_type, group=attribute_type))+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCredInt,ymax=upperCredInt,group=attribute_type,colour=attribute_type),size=1)+
  scale_color_manual(name="Attribute group", 
                     labels=c("diversity"="Diversity", "structure"="Structure","trait"="Function","soil"="Soil"),
                     values=c("soil"="#a6611a","diversity"="#80cdc1","trait"="#7b3294","structure"="#008837"))+
  
  labs(y=expression(bold("Recovery rate"~(lambda))),x="")+
  
  #scale_y_continuous(breaks=c(0,25,50,75,100,125))+
  #coord_cartesian(ylim = c(0, 130))+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),
                   axis.text.x = element_text(face="bold",size=8),
                   axis.text.y = element_text(face="bold",size=8),
                   legend.position = c(0.90,0.80),
                   legend.title = element_text(size=8,face="bold"),
                   legend.text = element_text(size=7),
                   panel.background = element_rect(fill="transparent", color=NA))

ggsave(lambda_rec,filename = file.path("Results", "MultiD_LambdaRecovery_SI.png"), 
       width    = 15, 
       height   = 10, 
       units    = "cm")