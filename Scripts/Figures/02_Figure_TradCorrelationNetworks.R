#-----------------------# 
# Visualise traditional #
# correlation network   #
#-----------------------#
require(tidyverse)
require(corrr)
require(networktools)
require(igraph)
require(ggraph)
require(tidygraph)
require(huge)
library(bayestestR)

load(file="Data/TradCorr_Bootstrap.RData")

#-------------------------#
# make network & estimate #
# expected influence      #
#-------------------------#

# 7 traits

sumCorr_7trts<-Corr_iters %>% 
  mutate(x= ifelse(x=="BIOM","AGB",x),
         y= ifelse(y=="BIOM","AGB",y)) %>% 
  unite("Comb",c(x,y),sep="_") %>% 
  group_by(Comb) %>% 
  summarize(median_r=median(r),
            l95=ci(r, method="HDI", ci=0.95)$CI_low,
            u95=ci(r, method="HDI", ci=0.95)$CI_high) %>%  
  dplyr::filter(., (l95 < 0 & u95 < 0|l95 > 0 & u95 > 0)) %>% 
  group_by(Comb) %>% 
  slice_head(., n=1) %>% 
  separate(Comb, c("x","y"),sep="_") %>% 
  mutate(dupe=duplicated(median_r)) %>% 
  dplyr::filter(., dupe==FALSE)

sumCorr_7trts<-sumCorr_7trts %>% 
  mutate(x_group=ifelse(x=="BIOM","Structure",NA),
         x_group=ifelse(x=="DMAX","Structure",x_group),
         x_group=ifelse(x=="SH","Structure",x_group),
         x_group=ifelse(x=="NF","Function",x_group),
         x_group=ifelse(x=="SD","Diversity",x_group)) %>%
  mutate(y_group=ifelse(y=="SR","Diversity",NA),
         y_group=ifelse(y=="DMAX","Structure",y_group),
         y_group=ifelse(y=="SH","Structure",y_group),
         y_group=ifelse(y=="NF","Function",y_group),
         y_group=ifelse(y=="SD","Diversity",y_group),
         y_group=ifelse(y=="WD","Function",y_group)) %>% 
  select(., x, y, y_group, x_group, median_r, l95, u95) 

CorrNet_c<-graph_from_data_frame(sumCorr_7trts, directed = F, vertices =NULL)
E(CorrNet_c)$width <- (sumCorr_7trts$median_r)
edge.attributes(CorrNet_c)$weight<-sumCorr_7trts$median_r

# communities

TradCorr_7trt<-cluster_spinglass(CorrNet_c, weights=E(CorrNet_c)$width, implementation="neg")

# cluster 1= "AGB"  "DMAX" "NF"   "SH" 
# cluster 2= "SD" "SR" "WD"

CorrNet_cc <- as_tbl_graph(CorrNet_c)
CorrNet_cc<-CorrNet_cc %>% 
  activate(edges) %>%
  activate(nodes) 

#layout

layout7_x<-c(-0.2225209,0.6234898,-0.9009689 ,-0.2225209, 1.0000000,0.6234898,-0.9009689)
layout7_y<-c(0.9749279,0.7818315,0.4338837, -0.9749279, 0.0000000,-0.7818315,-0.4338837)
layout7<-cbind.data.frame(x=layout7_x,y=layout7_y)
layout7$order<-c(1,2,3,4,5,6,7)
layout7$node<-c("AGB","DMAX","NF","SD","SH","SR","WD")


CorrNet7trts_p<-ggraph(CorrNet_cc,layout = layout7) + 
  # geom_edge_link(aes(width = abs(weight),label=round(abs(weight),digits=2)), 
  #                angle_calc="along",label_colour = '#2ca25f',label_dodge=unit(0.35,"cm"),
  #                colour="#2ca25f", alpha = 0.8) + 
  # geom_edge_link(aes(width = abs(mean)), colour ="grey70",alpha = 0.8) + 
  geom_edge_arc(aes(width = abs(median_r), color=abs(median_r)), strength=0.30)+
  #scale_edge_colour_viridis()
  scale_edge_colour_continuous(low="grey90",high="grey30")+
  
  
  #scale_edge_width(range = c(0.25, 2)) +
  geom_node_point(size = 14, colour = c("#008837", "#008837","#7b3294","#80cdc1","#008837","#80cdc1","#7b3294")) +
  geom_node_text(aes(label = name), colour = 'black', fontface="bold", repel=TRUE,
                 
                 vjust = c(-0.7),
                 hjust=c(-2.2)) + 
  theme_graph()+ theme(legend.position="none",
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank())


# Expected influence

ExpInf_iterss<-ExpInf_iters %>% 
  mutate(Attribute= ifelse(Attribute=="BIOM","AGB",Attribute)) %>% 
  group_by(Attribute) %>% 
  summarize(median=median(Expected_infl),
            l95=ci(Expected_infl, method="HDI", ci=0.95)$CI_low,
            u95=ci(Expected_infl, method="HDI", ci=0.95)$CI_high) %>% 
  mutate(group=ifelse(Attribute=="AGB","Structure",NA),
         group=ifelse(Attribute=="SH","Structure",group),
         group=ifelse(Attribute=="DMAX","Structure",group),
         group=ifelse(Attribute=="SR","Diversity",group),
         group=ifelse(Attribute=="SD","Diversity",group),
         group=ifelse(Attribute=="NF","Function",group),
         group=ifelse(Attribute=="WD","Function",group))

ExpInf_p<-ggplot(ExpInf_iterss, aes(x=reorder(Attribute,-median),y=median,group=group,colour=group,fill=group))+
  geom_hline(yintercept = 0, lty=2, colour="grey70")+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=l95, ymax=u95), size=1,fatten=5) +
  #geom_bar(stat="identity")+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function"),
                     values=c("Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+
  scale_fill_manual(name="Attribute group", 
                    labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function"),
                    values=c("Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+
  
  labs(x="",y="Expected influence")+
  theme_bw()+theme(axis.title = element_text(face="bold",size=10),
                   axis.text.x = element_text(face="bold",size=10),
                   axis.text.y = element_text(face="bold",size=10),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=10,face="bold"),
                   legend.text = element_text(size=10),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())

# together

net_tog<-cowplot::plot_grid(CorrNet7trts_p,
                            ExpInf_p+theme(legend.position="none"),
                            labels=c("A","B"),label_size=12,ncol=1)


legend<-cowplot::get_legend(ExpInf_p+theme( legend.direction = "horizontal",
                                            legend.text=element_text(face="bold"),
                                            legend.position=c(0.5,0.5),
                                            legend.title=element_blank()))

net_tog<-cowplot::plot_grid(legend,net_tog,  rel_heights = c(0.025,1),ncol=1)

ggsave(net_tog, file="Results/TradCorr_Network_NodeMetrics_7traits.png",
       height=30, 
       width=15,
       units="cm",
       dpi=500)

save(CorrNet7trts_p,ExpInf_p, file="Data/Network_NodeMetrics_Network_NodeMetrics_7traits.RData")

# All traits together

sumCorr_AllTrts<-CorrAllTrts_iters %>% 
  unite("Comb",c(x,y),sep="_") %>% 
  group_by(Comb) %>% 
  summarize(median_r=median(r),
            l95=ci(r, method="HDI", ci=0.95)$CI_low,
            u95=ci(r, method="HDI", ci=0.95)$CI_high) %>%  
  dplyr::filter(., (l95 < 0 & u95 < 0|l95 > 0 & u95 > 0)) %>% 
  group_by(Comb) %>% 
  slice_head(., n=1) %>% 
  separate(Comb, c("x","y"),sep="_") %>% 
  mutate(dupe=duplicated(median_r)) %>% 
  dplyr::filter(., dupe==FALSE)

CorrNetAll_c<-graph_from_data_frame(sumCorr_AllTrts, directed = F, vertices =NULL)
E(CorrNetAll_c)$width <- (sumCorr_AllTrts$median_r)
edge.attributes(CorrNetAll_c)$weight<-sumCorr_AllTrts$median_r

TradCorr_Alltrt<-cluster_spinglass(CorrNetAll_c, weights=E(CorrNetAll_c)$width, implementation="neg")

# cluster 1= "AGB"  "C"    "DMAX" "N"    "SH"
# cluster 2= "BD" "NF"
# cluster 3= "SC"  "SD"  "SR"  "SLA" "WD"

CorrNetAll_cc <- as_tbl_graph(CorrNetAll_c)
CorrNetAll_cc<-CorrNetAll_cc %>% 
  activate(edges) %>%
  activate(nodes)

layout12_x<-c(6.123234e-17,-1.000000e+00, -8.660254e-01,  5.000000e-01,-5.000000e-01,-5.000000e-01 ,5.000000e-01,8.660254e-01,8.660254e-01  ,1.000000e+00,-8.660254e-01, -1.836970e-16 )
layout12_y<-c(1.000000e+00, 1.224647e-16 ,5.000000e-01, 8.660254e-01,8.660254e-01,-8.660254e-01,-8.660254e-01,-5.000000e-01,5.000000e-01,0.000000e+00,-5.000000e-01,-1.000000e+00)
layout12<-cbind.data.frame(x=layout12_x,y=layout12_y)
layout12$order<-c(1,2,2,4,5,6,7,8,9,10,11,12)
layout12$node<-c("AGB","BD","C","DMAX","N","NF","SC","SD","SH","SR", "SLA","WD")


CorrNetAlltrts_p<-ggraph(CorrNetAll_cc,layout = layout12) + 
  # geom_edge_link(aes(width = abs(weight),label=round(abs(weight),digits=2)), 
  #                angle_calc="along",label_colour = '#2ca25f',label_dodge=unit(0.35,"cm"),
  #                colour="#2ca25f", alpha = 0.8) + 
  #geom_edge_link(aes(width = abs(median_r)), colour ="grey70",alpha = 0.8) + 
  geom_edge_arc(aes(width = abs(median_r), color=abs(median_r)), strength=0.30)+
  #scale_edge_colour_viridis()
  scale_edge_colour_continuous(low="grey90",high="grey30")+
  
  
  #scale_edge_width(range = c(0.25, 2)) +
  geom_node_point(size = 14, colour = c("#008837","#a6611a", "#a6611a","#008837","#a6611a","#7b3294","#80cdc1","#80cdc1","#008837","#80cdc1",
                                        "#7b3294","#7b3294")) +
  geom_node_text(aes(label = name),colour = 'black', fontface="bold", repel=TRUE,
                 
                 vjust = c(0.7),
                 hjust=c(-2.5)) + 
  theme_graph()+ theme(legend.position="none",
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank())



# Expected influence

Expinf_AllTrtss<-ExpInfAllTrts_iters %>% 
  group_by(Attribute) %>% 
  summarize(median=median(Expected_infl),
            l95=ci(Expected_infl, method="HDI", ci=0.95)$CI_low,
            u95=ci(Expected_infl, method="HDI", ci=0.95)$CI_high) %>%  
  mutate(group=ifelse(Attribute=="AGB","Structure",NA),
         group=ifelse(Attribute=="SH","Structure",group),
         group=ifelse(Attribute=="DMAX","Structure",group),
         group=ifelse(Attribute=="SR","Diversity",group),
         group=ifelse(Attribute=="SD","Diversity",group),
         group=ifelse(Attribute=="NF","Function",group),
         group=ifelse(Attribute=="WD","Function",group),
         group=ifelse(Attribute=="C","Soil",group),
         group=ifelse(Attribute=="N","Soil",group),
         group=ifelse(Attribute=="BD","Soil",group),
         group=ifelse(Attribute=="SLA","Function",group),
         group=ifelse(Attribute=="SC","Diversity",group))


ExpInf_AllTraits_p<-ggplot(Expinf_AllTrtss, aes(x=reorder(Attribute,-median),y=median,group=group,colour=group,fill=group))+
  geom_hline(yintercept = 0, lty=2, colour="grey70")+
  geom_point(size=2)+
  geom_pointrange(aes(ymin=l95, ymax=u95), size=1,fatten=5) +
  #geom_bar(stat="identity")+
  scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+
  scale_fill_manual(name="Attribute group", 
                    labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                    values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))+
  labs(x="",y="Expected influence")+
  theme_bw()+theme(axis.title = element_text(face="bold",size=10),
                   axis.text.x = element_text(face="bold",size=10),
                   axis.text.y = element_text(face="bold",size=10),
                   legend.position = c(0.90,0.87),
                   legend.title = element_text(size=10,face="bold"),
                   legend.text = element_text(size=10),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank())

# together

net_tog2<-cowplot::plot_grid(CorrNetAlltrts_p,
                             ExpInf_AllTraits_p+theme(legend.position="none"),
                             labels=c("A","B"),label_size=12,ncol=1)


legendd<-cowplot::get_legend(ExpInf_AllTraits_p+theme( legend.direction = "horizontal",
                                                       legend.text=element_text(face="bold"),
                                                       legend.position=c(0.5,0.5),
                                                       legend.title=element_blank()))

net_tog2<-cowplot::plot_grid(legendd,net_tog2,  rel_heights = c(0.025,1),ncol=1)

ggsave(net_tog2, file="Results/TradCorr_Network_NodeMetrics_AllTraits.png",
       height=30, 
       width=15,
       units="cm",
       dpi=500)

save(CorrNetAlltrts_p,ExpInf_AllTraits_p, file="Data/TradCorr_Network_NodeMetrics_AllTraits.RData")

# print out correlations (full network)

Corr_AllTrts<-CorrAllTrts_iters %>% 
  unite("Comb",c(x,y),sep="_") %>% 
  group_by(Comb) %>% 
  summarize(median=median(r),
            l95=ci(r, method="HDI", ci=0.95)$CI_low,
            u95=ci(r, method="HDI", ci=0.95)$CI_high) %>%  
  #  dplyr::filter(., (l95 < 0 & u95 < 0|l95 > 0 & u95 > 0)) %>% 
  group_by(Comb) %>% 
  slice_head(., n=1) %>% 
  separate(Comb, c("x","y"),sep="_") %>% 
  mutate(dupe=duplicated(median)) %>% 
  dplyr::filter(., dupe==FALSE) %>% 
  select(., "Attribute (x)"=x, "Attribute (y)"=y, r=median, l95,u95) %>% 
  mutate(conf_interv= paste("(", round(l95, digits=2)," - ",round(u95,digits=2),")",sep=""),
         r=round(r, digits=2)) %>% 
  select(.,-l95,-u95) %>% 
  arrange(`Attribute (x)`)

# sample size

pcorr_sampleN<-Hmisc::rcorr(as.matrix(Attributes_20y_Alltrts[,3:15]))

pcorr_sampleN<-setNames(reshape2::melt(pcorr_sampleN$n), c('x', 'y', 'Sample_N'))

pcorr_sampleN<-pcorr_sampleN %>% 
  mutate(x=as.character(x),
         y=as.character(y)) %>% 
  mutate(x=ifelse(x=="nfix","NF",x),
         y=ifelse(y=="nfix","NF",y),
         x=ifelse(x=="q0","SR",x),
         y=ifelse(y=="q0","SR",y),
         x=ifelse(x=="q2","SD",x),
         y=ifelse(y=="q2","SD",y),
         x=ifelse(x=="chao","SC",x),
         y=ifelse(y=="chao","SC",y),
         x=ifelse(x=="AGB","AGB",x),
         y=ifelse(y=="AGB","AGB",y),
         x=ifelse(x=="gini","SH",x),
         y=ifelse(y=="gini","SH",y),
         x=ifelse(x=="maxDBH","DMAX",x),
         y=ifelse(y=="maxDBH","DMAX",y),
         x=ifelse(x=="wd","WD",x),
         y=ifelse(y=="wd","WD",y), 
         x=ifelse(x=="sla","SLA",x),
         y=ifelse(y=="sla","SLA",y),
         x=ifelse(x=="soilbd","BD",x),
         y=ifelse(y=="soilbd","BD",y),
         x=ifelse(x=="C_vol","C",x),
         y=ifelse(y=="C_vol","C",y),
         x=ifelse(x=="N_vol","N",x),
         y=ifelse(y=="N_vol","N",y)) %>% 
  dplyr::filter(., !x=="P_vol") %>% 
  dplyr::filter(., !y=="P_vol") %>% 
  select(.,"Attribute (x)"=x, "Attribute (y)"=y, N=Sample_N )

# join

Corr_AllTrts<-left_join(Corr_AllTrts,pcorr_sampleN, by.y=c("Attribute (x)","Attribute (y)") )

#write.csv(Corr_AllTrts, "Data/AllAttributes_TradPairwiseCorr.csv",row.names=F)

