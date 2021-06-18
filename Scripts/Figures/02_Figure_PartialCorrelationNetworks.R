#----------------------------#
# Recovery rates 20y         #
# network and metric figures #  
#----------------------------#

require(igraph)
require(ggraph)
require(tidyverse)
require(tidygraph)
library(bayestestR)

# Data

load(file= "Data/Networks_Metrics_RecovRates_20y.RData")

# 20 y / 7 traits / 77 sites

Overr<-summary(boot_20y7trt) %>% 
  dplyr::filter(., type=="edge") %>% 
  select(., node1, node2, mean, CIlower, CIupper) %>% 
  mutate(node1=ifelse(node1=="nfix","NF",node1),
         node2=ifelse(node2=="nfix","NF",node2),
         node1=ifelse(node1=="q0","SR",node1),
         node2=ifelse(node2=="q0","SR",node2),
         node1=ifelse(node1=="q2","SD",node1),
         node2=ifelse(node2=="q2","SD",node2),
         node1=ifelse(node1=="AGB","AGB",node1),
         node2=ifelse(node2=="AGB","AGB",node2),
         node1=ifelse(node1=="gini","SH",node1),
         node2=ifelse(node2=="gini","SH",node2),
         node1=ifelse(node1=="maxDBH","DMAX",node1),
         node2=ifelse(node2=="maxDBH","DMAX",node2),
         node1=ifelse(node1=="wd","WD",node1),
         node2=ifelse(node2=="wd","WD",node2)) %>% 
  arrange(.,mean)

Overr$type<-NULL

Overr<-dplyr::select(Overr, node1, node2, mean)

Overr_c<-graph_from_data_frame(Overr, directed = F, vertices =NULL)
E(Overr_c)$width <- Overr$mean
edge.attributes(Overr_c)$weight<-Overr$mean

PartCorr_7trts<-cluster_spinglass(Overr_c, weights=E(Overr_c)$width, implementation="neg")

# cluster 1= SH, DMAX, BIOM
# cluster 2= SD, SR
# cluster 3= NF, 
# cluster 3= WD 

Overr_cc <- as_tbl_graph(Overr_c)
Overr_cc<-Overr_cc %>% 
  activate(edges) %>%
  activate(nodes)

# make layout

layout7_x<-c(-0.2225209,0.6234898,0.6234898,-0.2225209, 1.0000000,-0.9009689,-0.9009689)
layout7_y<-c(0.9749279,-0.7818315,0.7818315,-0.9749279, 0.0000000,-0.4338837,0.4338837)
layout7<-cbind.data.frame(x=layout7_x,y=layout7_y)
layout7$order<-c(1,2,3,4,5,6,7)
layout7$node<-c("AGB","SR","DMAX","SD","SH","NF","WD")

#1.0000000  0.0000000  BIOM   
#0.6234898  0.7818315    SR
#-0.2225209  0.9749279    DMAX 
#-0.9009689  0.4338837  SD 
#-0.9009689 -0.4338837    SH 
#-0.2225209 -0.9749279    NF 
#0.6234898 -0.7818315    WD  
  
Over_c<-ggraph(Overr_cc,layout = layout7) + 
  # geom_edge_link(aes(width = abs(weight),label=round(abs(weight),digits=2)), 
  #                angle_calc="along",label_colour = '#2ca25f',label_dodge=unit(0.35,"cm"),
  #                colour="#2ca25f", alpha = 0.8) + 
  # geom_edge_link(aes(width = abs(mean)), colour ="grey70",alpha = 0.8) + 
  geom_edge_arc(aes(width = abs(mean), color=abs(mean)), strength=0.30)+
  #scale_edge_colour_viridis()
  scale_edge_colour_continuous(low="transparent",high="grey30")+
  
  
  #scale_edge_width(range = c(0.25, 2)) +
  geom_node_point(size = 14, colour = c("#008837", "#80cdc1","#008837","#80cdc1","#008837","#7b3294","#7b3294")) +
  geom_node_text(aes(label = name), colour = 'black', fontface="bold", repel=TRUE,
                 
                 vjust = c(0.7),
                 hjust=c(-2.5)) + 
  theme_graph()+ theme(legend.position="none",
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank())

# network metrics

metrics<-data.frame((boot_20y7trt$bootTable)) %>% 
  dplyr::filter(., type=="expectedInfluence") %>% 
  select(., node1, type,value) %>% 
  group_by(node1,type) %>% 
  summarize(median=median(value),
            lowerCI=ci(value, method="HDI", ci=0.95)$CI_low,
            upperCI=ci(value, method="HDI", ci=0.95)$CI_high)%>% 
  mutate(group=ifelse(node1=="AGB","Structure",NA),
         group=ifelse(node1=="gini","Structure",group),
         group=ifelse(node1=="maxDBH","Structure",group),
         group=ifelse(node1=="q0","Diversity",group),
         group=ifelse(node1=="q2","Diversity",group),
         group=ifelse(node1=="nfix","Function",group),
         group=ifelse(node1=="WD","Function",group)) %>% 
  mutate(node1=ifelse(node1=="nfix","NF",node1),
         node1=ifelse(node1=="q0","SR",node1),
         node1=ifelse(node1=="q2","SD",node1),
         node1=ifelse(node1=="AGB","AGB",node1),
         node1=ifelse(node1=="gini","SH",node1),
         node1=ifelse(node1=="maxDBH","DMAX",node1),
         node1=ifelse(node1=="WD","WD",node1))

node_g<-ggplot(metrics, aes(x=reorder(node1,-median),y=median,group=group,colour=group,fill=group))+
  geom_hline(yintercept = 0, lty=2, colour="grey70")+
  #geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCI, ymax=upperCI), size=1,fatten=5) +
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

net_tog<-cowplot::plot_grid(Over_c,
                            node_g+theme(legend.position="none"),
                            labels=c("A","B"),label_size=12,ncol=1)


legend<-cowplot::get_legend(node_g+theme( legend.direction = "horizontal",
                                          legend.text=element_text(face="bold"),
                                          legend.position=c(0.5,0.5),
                                          legend.title=element_blank()))

net_tog<-cowplot::plot_grid(legend,net_tog,  rel_heights = c(0.025,1),ncol=1)

ggsave(net_tog, file="Results/Network_NodeMetrics_RecovRates20_67Sites.png",
       height=30, 
       width=15,
       units="cm",
       dpi=500)

save(Over_c,node_g, file="Data/Network_NodeMetrics_RecovRates20_67Sites.RData")

#-------------------------------#
# 20y / 8 attributes / 43 sites #
#-------------------------------#

Overr<-summary(boot_20y8trt) %>% 
  dplyr::filter(., type=="edge") %>% 
  select(., node1, node2, mean, CIlower, CIupper) %>% 
  mutate(node1=ifelse(node1=="nfix","NF",node1),
         node2=ifelse(node2=="nfix","NF",node2),
         node1=ifelse(node1=="q0","SR",node1),
         node2=ifelse(node2=="q0","SR",node2),
         node1=ifelse(node1=="q2","SD",node1),
         node2=ifelse(node2=="q2","SD",node2),
         node1=ifelse(node1=="chao","SC",node1),
         node2=ifelse(node2=="chao","SC",node2),
         node1=ifelse(node1=="AGB","AGB",node1),
         node2=ifelse(node2=="AGB","AGB",node2),
         node1=ifelse(node1=="gini","SH",node1),
         node2=ifelse(node2=="gini","SH",node2),
         node1=ifelse(node1=="maxDBH","DMAX",node1),
         node2=ifelse(node2=="maxDBH","DMAX",node2),
         node1=ifelse(node1=="WD","WD",node1),
         node2=ifelse(node2=="WD","WD",node2)) 

Overr$type<-NULL

Overr<-dplyr::select(Overr, node1, node2, mean) 

Overr_c<-graph_from_data_frame(Overr, directed = F, vertices =NULL)
E(Overr_c)$width <- Overr$mean
edge.attributes(Overr_c)$weight<-Overr$mean

PartCorr_8trts<-cluster_spinglass(Overr_c, weights=E(Overr_c)$width, implementation="neg")

# cluster 1= SH, DMAX, BIOM
# cluster 2= SD, SR
# cluster 3= NF, SC,WD

Overr<-dplyr::select(Overr, node1, node2, mean) %>% 
         mutate(mean=ifelse(node1=="SC"& node2=="WD",0,mean),
         mean=ifelse(node1=="SR"& node2=="SD",0,mean))

Overr_c<-graph_from_data_frame(Overr, directed = F, vertices =NULL)
E(Overr_c)$width <- Overr$mean
edge.attributes(Overr_c)$weight<-Overr$mean

Overr_cc <- as_tbl_graph(Overr_c)
Overr_cc<-Overr_cc %>% 
  activate(edges) %>%
  activate(nodes)

layout8_x<-c(6.123234e-17,-7.071068e-01,1.000000e+00, 7.071068e-01 ,-1.000000e+00 ,7.071068e-01,1.836970e-16, -7.071068e-01 )
layout8_y<-c(1.000000e+00,-7.071068e-01,0.000000e+00, 7.071068e-01, 1.224647e-16 ,-7.071068e-01,-1.000000e+00,  7.071068e-01)
layout8<-cbind.data.frame(x=layout8_x,y=layout8_y)
layout8$order<-c(1,2,2,4,5,6,7,8)
layout8$node<-c("AGB","SC","SH","DMAX","NF","SR","SD","WD")

Over_cc<-ggraph(Overr_cc,layout =layout8) + 
  # geom_edge_link(aes(width = abs(weight),label=round(abs(weight),digits=2)), 
  #                angle_calc="along",label_colour = '#2ca25f',label_dodge=unit(0.35,"cm"),
  #                colour="#2ca25f", alpha = 0.8) + 
  # geom_edge_link(aes(width = abs(mean)), colour ="grey70",alpha = 0.8) + 
  geom_edge_arc(aes(width = abs(mean), color=abs(mean)), strength=0.30)+
  #scale_edge_colour_viridis()
  scale_edge_colour_continuous(low="transparent",high="grey30")+
  
  
  #scale_edge_width(range = c(0.25, 2)) +
  geom_node_point(size = 14, colour = c("#008837", "#80cdc1","#008837","#008837","#7b3294","#80cdc1","#80cdc1", "#7b3294")) +
  geom_node_text(aes(label = name), colour = 'black', fontface="bold", repel=TRUE,
                 vjust = c(0.7),
                 hjust=c(-2.5)) + 
  theme_graph()+ theme(legend.position="none",
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank())

# network metrics

metrics<-data.frame((boot_20y8trt$bootTable)) %>% 
  dplyr::filter(., type=="expectedInfluence") %>% 
  select(., node1, type,value) %>% 
  group_by(node1,type) %>% 
  summarize(median=median(value),
            lowerCI=ci(value, method="HDI", ci=0.95)$CI_low,
            upperCI=ci(value, method="HDI", ci=0.95)$CI_high)%>% 
  mutate(group=ifelse(node1=="AGB","Structure",NA),
         group=ifelse(node1=="gini","Structure",group),
         group=ifelse(node1=="maxDBH","Structure",group),
         group=ifelse(node1=="q0","Diversity",group),
         group=ifelse(node1=="q2","Diversity",group),
         group=ifelse(node1=="nfix","Function",group),
         group=ifelse(node1=="WD","Function",group),
         group=ifelse(node1=="chao","Diversity",group)) %>% 
  mutate(node1=ifelse(node1=="nfix","NF",node1),
         node1=ifelse(node1=="q0","SR",node1),
         node1=ifelse(node1=="q2","SD",node1),
         node1=ifelse(node1=="chao","SC",node1),
         node1=ifelse(node1=="AGB","BIOM",node1),
         node1=ifelse(node1=="gini","SH",node1),
         node1=ifelse(node1=="maxDBH","DMAX",node1),
         node1=ifelse(node1=="WD","WD",node1))

node_gg<-ggplot(metrics, aes(x=reorder(node1,-median),y=median,group=group,colour=group,fill=group))+
  geom_hline(yintercept = 0, lty=2, colour="grey70")+
  #geom_point(size=2)+
  geom_pointrange(aes(ymin=lowerCI, ymax=upperCI), size=1,fatten=5) +
  
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

net_togg<-cowplot::plot_grid(Over_cc,
                            node_gg+theme(legend.position="none"),
                            labels=c("A","B"),label_size=12,ncol=1)


legendd<-cowplot::get_legend(node_gg+theme( legend.direction = "horizontal",
                                          legend.text=element_text(face="bold"),
                                          legend.position=c(0.5,0.5),
                                          legend.title=element_blank()))

net_togg<-cowplot::plot_grid(legendd,net_togg,  rel_heights = c(0.025,1),ncol=1)

ggsave(net_togg, file="Results/Network_NodeMetrics_RecovRates20_43Sites.png",
       height=30, 
       width=15,
       units="cm",
       dpi=500)

save(Over_cc,node_gg, file="Data/Network_NodeMetrics_RecovRates20_50Sites8attr.RData")
