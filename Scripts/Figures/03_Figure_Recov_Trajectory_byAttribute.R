#------------------------#
# Organize recovery data #
# lambda
# - from B. Herault      #
#------------------------#

require(tidyverse)
require(directlabels)
require(data.table)
#-------------------------#
#         Data            #
#-------------------------#

# Files too big for Github
# load("Data/DBH_distance.Rdata")
# load("Data/AGB_distance.Rdata")
# load("Data/gini_distance.Rdata")
# load("Data/C_vol_distance.Rdata")
# load("Data/chao_Bruno_distance.Rdata")
# load("Data/N_vol_distance.Rdata")
# load("Data/nfix_distance.Rdata")
# load("Data/P_vol_distance.Rdata")
# load("Data/q0_distance.Rdata")
# load("Data/q2_distance.Rdata")
# load("Data/sla_distance.Rdata")
# load("Data/wd_distance.Rdata")
# load("Data/Bulk_density_all_distance.Rdata")

#----------#
# Recovery #
#----------#

agb_recov<-AGB_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=AGB_KLD_pred, recov_lourens= AGB_Lourens_pred, recov_dylan= AGB_Dylan_pred) %>% 
  mutate(attribute="AGB") 

dbh_recov<-DBH_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=DBH_KLD_pred, recov_lourens= DBH_Lourens_pred, recov_dylan= DBH_Dylan_pred) %>% 
  mutate(attribute="maxDBH")
gini_recov<-gini_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=gini_KLD_pred, recov_lourens= gini_Lourens_pred, recov_dylan= gini_Dylan_pred) %>% 
  mutate(attribute="gini") 
C_vol_recov<-C_vol_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=C_vol_KLD_pred, recov_lourens= C_vol_Lourens_pred, recov_dylan= C_vol_Dylan_pred) %>% 
  mutate(attribute="C_vol")

sBD_recov<-Bulk_density_all_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=Bulk_density_all_KLD_pred, recov_lourens= Bulk_density_all_Lourens_pred, recov_dylan= Bulk_density_all_Dylan_pred) %>% 
  mutate(attribute="BD")

N_vol_recov<-N_vol_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=N_vol_KLD_pred, recov_lourens= N_vol_Lourens_pred, recov_dylan= N_vol_Dylan_pred) %>% 
  mutate(attribute="N_vol") 

chao_recov<-mh_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=mh_KLD_pred, recov_lourens= mh_Lourens_pred, recov_dylan= mh_Dylan_pred) %>% 
  mutate(attribute="chao") 

q0_recov<-q0_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=q0_KLD_pred, recov_lourens= q0_Lourens_pred, recov_dylan= q0_Dylan_pred) %>% 
  mutate(attribute="q0") 

q2_recov<-q2_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=q2_KLD_pred, recov_lourens= q2_Lourens_pred, recov_dylan= q2_Dylan_pred) %>% 
  mutate(attribute="q2") 

nfix_recov<-nfix_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=nfix_KLD_pred, recov_lourens= nfix_Lourens_pred, recov_dylan= nfix_Dylan_pred) %>% 
  mutate(attribute="nfix") 

sla_recov<-sla_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=sla_KLD_pred, recov_lourens= sla_Lourens_pred, recov_dylan= sla_Dylan_pred) %>% 
  mutate(attribute="SLA")

wd_recov<-wd_distance_coeffs %>% 
  select(., chronosequence=curve, age, recov_kld=wd_KLD_pred, recov_lourens= wd_Lourens_pred, recov_dylan= wd_Dylan_pred) %>% 
  mutate(attribute="WD")

recov_b<-bind_rows(agb_recov,dbh_recov,gini_recov, chao_recov,q0_recov,q2_recov, sla_recov,wd_recov, nfix_recov, C_vol_recov,sBD_recov, N_vol_recov) %>% 
  mutate(attribute_type=ifelse(attribute=="AGB","Structure",NA),
         attribute_type=ifelse(attribute=="maxDBH","Structure",attribute_type),
         attribute_type=ifelse(attribute=="gini","Structure",attribute_type),
         attribute_type=ifelse(attribute=="chao","Diversity",attribute_type),
         attribute_type=ifelse(attribute=="q0","Diversity",attribute_type),
         attribute_type=ifelse(attribute=="q2","Diversity",attribute_type),
         attribute_type=ifelse(attribute=="SLA","Function",attribute_type),
         attribute_type=ifelse(attribute=="WD","Function",attribute_type),
         attribute_type=ifelse(attribute=="nfix","Function",attribute_type),
         attribute_type=ifelse(attribute=="C_vol","Soil",attribute_type),
         attribute_type=ifelse(attribute=="BD","Soil",attribute_type),
         attribute_type=ifelse(attribute=="N_vol","Soil",attribute_type)) %>% 

  mutate(attribute=ifelse(attribute=="BD","BD",attribute),
         attribute=ifelse(attribute=="N_vol","N",attribute),
         attribute=ifelse(attribute=="C_vol","C",attribute),
         attribute=ifelse(attribute=="sla","SLA",attribute),
         attribute=ifelse(attribute=="WD","WD",attribute),
         attribute=ifelse(attribute=="nfix","NF",attribute),
         attribute=ifelse(attribute=="gini","SH",attribute),
         attribute=ifelse(attribute=="maxDBH","DMAX",attribute),
         attribute=ifelse(attribute=="q0","SR",attribute),
         attribute=ifelse(attribute=="q2","SD",attribute),
         attribute=ifelse(attribute=="chao","SC",attribute),
         attribute=ifelse(attribute=="AGB","AGB",attribute)) %>% 
  mutate(recov_kld=recov_kld*100)

recov_b$attribute<-factor(recov_b$attribute, levels=c("AGB", "DMAX", "SH", "SR", "SD", "SC", "NF","SLA", "WD","BD", "C","N"))

save(recov_b, file="Data/RecovCurvEst.RData")

recov_bb<-recov_b %>% 
          group_by(attribute_type, attribute, age) %>% 
          summarize(recov_kldd=median(recov_kld))

# visualize

curves<-ggplot(recov_bb, aes(x=age, y=recov_kldd,colour=attribute_type))+
        geom_vline(xintercept=20,color="grey70",lty=2)+
        geom_vline(xintercept=40,color="grey70",lty=2)+
        geom_hline(yintercept=90,color="grey70",lty=2)+
        geom_line(aes(group=attribute))+
       #geom_smooth(aes(group=attribute),method = "gam", formula=y~s(x), se=FALSE)+
       scale_color_manual(name="Attribute group", 
                     labels=c("Diversity"="Diversity", "Structure"="Structure","Function"="Function","Soil"="Soil"),
                     values=c("Soil"="#a6611a","Diversity"="#80cdc1","Function"="#7b3294","Structure"="#008837"))

smooth_dat <- setDT(ggplot_build(curves)$data[[4]])
smooth_lab <- smooth_dat[smooth_dat[, .I[x == min(x)], by=colour]$V1]

smooth_lab$label <- levels(recov_bb$attribute)

curvess<-curves + annotate("text", x = smooth_lab$x, y=smooth_lab$y, 
             label=smooth_lab$label,colour=smooth_lab$colour,
             vjust=c(1.5,0.5,0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0.8,0.8,0.8),hjust=c(1,1,1,1,1,1,1,1,1,1,1,1))+

        labs(y="Relative recovery (%)",x="Time (y)")+
        scale_y_continuous(breaks=c(0,20,40,60,80,90,100))+
        scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+
  
        coord_cartesian(ylim = c(0, 105),xlim=c(-3,120))+
        theme_bw()+theme(axis.title = element_text(face="bold",size=14),
                   axis.text.x = element_text(face="bold",size=14),
                   axis.text.y = element_text(face="bold",size=14),
                   legend.position = c(0.85,0.15),
                   legend.title = element_text(size=14,face="bold"),
                   legend.text = element_text(size=14),
                   panel.background = element_rect(fill="transparent", color=NA))

ggsave(curvess, filename = file.path("Results", "MultiD_RelRecovery_Curves.png"), 
       width    = 25, 
       height   = 20, 
       units    = "cm")
