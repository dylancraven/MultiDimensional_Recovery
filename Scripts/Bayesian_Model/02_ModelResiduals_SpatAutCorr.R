#--------------------------#
# Spatial auto-correlation #
#--------------------------#

require(tidyverse)
require(pgirmess)
require(sp)

#-------------------------#
#     Data                #
#-------------------------#

# Data too big for Github
# coord<-read.csv("Data/MDR_site_coordinates.csv")
# coord_soil<-coord %>% 
#             dplyr::filter(., variable=="soil")
# 
# coord_plant<-coord %>% 
#   dplyr::filter(., variable=="vegetation")
# 
# residuals<-read.delim("Data/data_nosoil_nosimilarity_Moran.csv",sep=";",dec = ",")
# 
# res_sim<-read.delim("Data/data_similarity_Moran.csv",sep=";",dec = ",")
# 
# res_soil<-read.delim("Data/data_soil_Moran.csv",sep=";",dec = ",")

# Aboveground attributes
resids<-residuals %>%
        mutate(chronosequence=as.factor(chronosequence),
               plot_2ndfor=as.factor(plot_2ndfor)) %>% 
        mutate_if(is.character,as.numeric) %>% 
        dplyr::select(.,chronosequence, plot_2ndfor, lat_all, long_all, 
                      AGB_Mgha,AGB_pred,
                      maxDBH.95q,dbh_pred,
                      plot_Gini,gini_pred,
                      hillq0_bc,q0_pred, hillq2_bc,q2_pred,
                      perc_Nfix,nfix_pred,
                      cwmSLA_BA,sla_pred,
                      cwmWD_BA,wd_pred) %>% 
                 mutate(AGB_resids= AGB_Mgha-AGB_pred,
                        DMAX_resids=maxDBH.95q-dbh_pred,
                        Gini_resids=plot_Gini-gini_pred,
                        SR_resids=hillq0_bc-q0_pred,
                        SE_resids=hillq2_bc-q2_pred,
                        nFix_resids=perc_Nfix-nfix_pred,
                        SLA_resids=cwmSLA_BA-sla_pred,
                        WD_resids=cwmWD_BA-wd_pred) %>% 
          left_join(., coord_plant,by="chronosequence") %>% 
          dplyr::select(chronosequence,lat_site=lat,long_site=long,
                        AGB=AGB_Mgha,AGB_resids,
                        DMAX=maxDBH.95q,DMAX_resids, 
                        SH=plot_Gini,SH_resids=Gini_resids,
                        SR=hillq0_bc,SR_resids,     
                        SD=hillq2_bc,SD_resids=SE_resids,
                        NF=perc_Nfix,NF_resids=nFix_resids,
                        SLA=cwmSLA_BA,SLA_resids,
                        WD=cwmWD_BA,WD_resids) %>% 
          pivot_longer(.,AGB:WD_resids, names_to="Attributes_Or_Residual",values_to="Value",
                       values_drop_na=TRUE) 

resid_veg<-res_sim %>%
  mutate(chronosequence=as.factor(chronosequence))%>% 
  mutate_if(is.character,as.numeric) %>% 
  dplyr::select(.,chronosequence,  
                medianDist.chao,medianDist.chao_pred) %>% 
  mutate(SC_resids= medianDist.chao-medianDist.chao_pred) %>% 
  left_join(., coord_plant,by="chronosequence") %>% 
  dplyr::select(chronosequence, lat_site=lat,long_site=long,
                SC=medianDist.chao,SC_resids) %>% 
  pivot_longer(.,SC:SC_resids, names_to="Attributes_Or_Residual",values_to="Value",
               values_drop_na=TRUE) %>% 
  bind_rows(., resids)


# soils
resid_soil<-res_soil %>%
  mutate(chronosequence=as.factor(chronosequence),
         plot=as.factor(plot),
         DEPTH=as.factor(DEPTH)) %>% 
  dplyr::filter(.,DEPTH=="0-15") %>% 
  mutate_if(is.character,as.numeric) %>% 
  dplyr::select(.,chronosequence, plot,  
                N_vol,N_vol_pred,
                C_vol,C_vol_pred,
                Bulk_density_all,
                Bulk_density_all_pred) %>% 
  mutate(N_resids= N_vol-N_vol_pred,
         C_resids=C_vol-C_vol_pred,
         BD_resids=Bulk_density_all-Bulk_density_all_pred) %>% 
  left_join(., coord_soil,by="chronosequence") %>% 
  dplyr::select(chronosequence,plot, lat_site=lat,long_site=long,
                N=N_vol,N_resids,
                C=C_vol,C_resids, 
                BD=Bulk_density_all,BD_resids) %>% 
  pivot_longer(.,N:BD_resids, names_to="Attributes_Or_Residual",values_to="Value",
               values_drop_na=TRUE) 

save(resid_veg,resid_soil, file="Data/Residual_fun.RData")

## Spatial autocorrelation (Moran's I)
#load(file="Data/Residual_fun.RData")
# Aboveground attributes

spat_autocor<-list()

for (i in 1:18){
  
 dat<-dplyr::filter(resid_veg, Attributes_Or_Residual==unique(Attributes_Or_Residual)[i])
 
 att<-as.character(unique(dat$Attributes_Or_Residual))
 
 coordz<-cbind(dat$long_site, dat$lat_site)

 ezt<-dat$Value

 corr<-correlog(coords=coordz, z=ezt, method="Moran",nbclass=15)
 corr<-data.frame(corr)

 out<- corr %>% 
       summarize(min_coef=min(coef),max_coef=max(coef), mean_coef=mean(coef)) %>% 
       mutate(Attributes_Or_Residual=att)   
 spat_autocor[[i]]<-rbind.data.frame(out)
}

veg_autocor<-bind_rows(spat_autocor)

# Soil attributes

spat_autocor_s<-list()

for (i in 1:6){
  
  dat<-dplyr::filter(resid_soil, Attributes_Or_Residual==unique(Attributes_Or_Residual)[i])
  att<-as.character(unique(dat$Attributes_Or_Residual))
  
  coordz<-cbind(dat$long_site, dat$lat_site)
  
  ezt<-dat$Value
  
  corr<-correlog(coords=coordz, z=ezt, method="Moran",nbclass=20)
  corr<-data.frame(corr)
  
  out<- corr %>% 
    summarize(min_coef=min(coef),max_coef=max(coef), mean_coef=mean(coef)) %>% 
    mutate(Attributes_Or_Residual=att)   
  spat_autocor_s[[i]]<-rbind.data.frame(out)
}

soil_autocor<-bind_rows(spat_autocor_s)

###

save(veg_autocor,soil_autocor, file="Data/Residual_fun.RData")

load(file="Data/Residual_fun.RData")

tog<-bind_rows(soil_autocor, veg_autocor) %>% 
     mutate(Model_Resids=str_detect(Attributes_Or_Residual, pattern="_resids", negate = FALSE),
            Attributes_Or_Residual=str_replace_all(Attributes_Or_Residual, "_resids",""),
            group=ifelse(Attributes_Or_Residual=="AGB","Structure",NA),
            group=ifelse(Attributes_Or_Residual=="SH","Structure",group),
            group=ifelse(Attributes_Or_Residual=="DMAX","Structure",group),
            group=ifelse(Attributes_Or_Residual=="SR","Diversity",group),
            group=ifelse(Attributes_Or_Residual=="SD","Diversity",group),
            group=ifelse(Attributes_Or_Residual=="NF","Function",group),
            group=ifelse(Attributes_Or_Residual=="WD","Function",group),
            group=ifelse(Attributes_Or_Residual=="C","Soil",group),
            group=ifelse(Attributes_Or_Residual=="N","Soil",group),
            group=ifelse(Attributes_Or_Residual=="BD","Soil",group),
            group=ifelse(Attributes_Or_Residual=="SLA","Function",group),
            group=ifelse(Attributes_Or_Residual=="SC","Diversity",group))

mod_residuals<-tog %>% 
               dplyr::filter(Model_Resids==TRUE) %>% 
               dplyr::select(., Attribute_Group=group,Attributes=Attributes_Or_Residual, meanMoransI=mean_coef, min_MoransI=min_coef,
                             max_MoransI=max_coef)

mod_values<-tog %>% 
  dplyr::filter(Model_Resids==FALSE) %>% 
  dplyr::select(., Attribute_Group=group,Attributes=Attributes_Or_Residual, meanMoransI=mean_coef, min_MoransI=min_coef,
                max_MoransI=max_coef)

write.csv(mod_residuals, "Data/ModelResiduals_MoransI.csv",row.names=FALSE)
write.csv(mod_values, "Data/ModelValues_MoransI.csv",row.names=FALSE)


