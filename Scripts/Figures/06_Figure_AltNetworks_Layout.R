#-----------------------------# 
# Network analysis            #
#  i) trad correlations       #
#      (all traits)           #
#  ii) partial corr. network  #
#-----------------------------#

require(tidyverse)
require(corrr)
require(networktools)
require(igraph)
require(ggraph)
require(tidygraph)
require(huge)
library(bayestestR)

# load data

load(file="Data/RecovRates.RData")
load(file="Data/TradCorr_Bootstrap.RData")

#-------------------------#
# Traditional network     #
#-------------------------#

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


CorrNetAll_cc <- as_tbl_graph(CorrNetAll_c)
CorrNetAll_cc<-CorrNetAll_cc %>% 
  activate(edges) %>%
  activate(nodes)


CorrNetAlltrts_p<-ggraph(CorrNetAll_cc,layout = "dh", cool.fact = 0.5) + 
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

#-------------------------#
# Partial corr.  network  #
#-------------------------#
# 20 y / 7 traits / 77 sites


load(file= "Data/Networks_Metrics_RecovRates_20y.RData")


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


Overr_cc <- as_tbl_graph(Overr_c)
Overr_cc<-Overr_cc %>% 
  activate(edges) %>%
  activate(nodes)


PartialCorr_7Alltrts_p<-ggraph(Overr_cc,layout = "dh", cool.fact=0.5) + 
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

# to get legend

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
                             PartialCorr_7Alltrts_p+theme(legend.position="none"),
                             labels=c("A","B"),label_size=12,ncol=1)


legendd<-cowplot::get_legend(ExpInf_AllTraits_p+theme( legend.direction = "horizontal",
                                                       legend.text=element_text(face="bold"),
                                                       legend.position=c(0.5,0.5),
                                                       legend.title=element_blank()))

net_tog2<-cowplot::plot_grid(legendd,net_tog2,  rel_heights = c(0.05,1),ncol=1)

ggsave(net_tog2, file="Results/Alt_Network_LayOut.png",
       height=30, 
       width=20,
       units="cm",
       dpi=500)

