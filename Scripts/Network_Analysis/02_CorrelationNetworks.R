#-----------------------------# 
# Correlation network analysis#
#                             #
# - using normal pearson      #
#   correlations              #
#-----------------------------#

require(tidyverse)
require(corrr)
require(networktools)
require(igraph)
require(ggraph)
require(tidygraph)
require(huge)
library(bayestestR)

load(file="Data/MDR_RecoveryRates_All_Attributes.RData")

#-------------------------#
# Transform   data        #
#-------------------------#

# 7 traits

corr_out<-list(); ExpInf_out<-list()

for (i in 1:10000){
  
  dat<-Attributes_20y_7trts %>% 
       slice_sample(., prop=0.80, replace=FALSE)

  dat_trn<-huge.npn(dat[,3:9], npn.func="truncation",verbose=F)
  
  dat_pcorr<-dat_trn %>% correlate(quiet=TRUE, use="pairwise.complete.obs")
  
  dat_pcorr<-stretch(dat_pcorr, na.rm = TRUE, remove.dups = FALSE) # now ready for network
  
  dat_pcorr<-dat_pcorr %>% 
    mutate(x=ifelse(x=="nfix","NF",x),
           y=ifelse(y=="nfix","NF",y),
           x=ifelse(x=="q0","SR",x),
           y=ifelse(y=="q0","SR",y),
           x=ifelse(x=="q2","SD",x),
           y=ifelse(y=="q2","SD",y),
           x=ifelse(x=="AGB","BIOM",x),
           y=ifelse(y=="AGB","BIOM",y),
           x=ifelse(x=="gini","SH",x),
           y=ifelse(y=="gini","SH",y),
           x=ifelse(x=="maxDBH","DMAX",x),
           y=ifelse(y=="maxDBH","DMAX",y),
           x=ifelse(x=="wd","WD",x),
           y=ifelse(y=="wd","WD",y)) 
  
  datNet_c<-graph_from_data_frame(dat_pcorr, directed = F, vertices =NULL)
  E(datNet_c)$width <- dat_pcorr$r
  edge.attributes(datNet_c)$weight<-dat_pcorr$r
  
  Expinf_Dat<-networktools::expectedInf(datNet_c)
  Expinf_Dats<-as.data.frame(matrix(Expinf_Dat$step1))
  Expinf_Dats<-cbind.data.frame(Expinf_Dats, Attribute=names(Expinf_Dat$step1))
  colnames(Expinf_Dats)[1]<-"Expected_infl"

  Expinf_Dats$iter<-i
  dat_pcorr$iter<-i
  
  corr_out[[i]]<-rbind.data.frame(dat_pcorr)
  ExpInf_out[[i]]<-rbind.data.frame(Expinf_Dats)
}
  
Corr_iters<-bind_rows(corr_out)
ExpInf_iters<-bind_rows(ExpInf_out)

#  All traits

corrAlltrts_out<-list(); ExpInf_Alltrts_out<-list()

for (i in 1:10000){
  
  dat<-Attributes_20y_Alltrts %>% 
    slice_sample(., prop=0.80, replace=FALSE)

  pcorr_Alltrts<-dat[,3:15] %>% correlate(quiet=TRUE, use="pairwise.complete.obs")
  pcorr_Alltrts<-stretch(pcorr_Alltrts, na.rm = TRUE, remove.dups = FALSE) # now ready for network
  

  pcorr_Alltrts<-pcorr_Alltrts %>% 
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
    dplyr::filter(., !y=="P_vol")
  
  CorrNetAll_c<-graph_from_data_frame(pcorr_Alltrts, directed = F, vertices =NULL)
  E(CorrNetAll_c)$width <- pcorr_Alltrts$r
  edge.attributes(CorrNetAll_c)$weight<-pcorr_Alltrts$r
  
  Expinf_AllTrts<-networktools::expectedInf(CorrNetAll_c)
  Expinf_AllTrtss<-as.data.frame(matrix(Expinf_AllTrts$step1))
  Expinf_AllTrtss<-cbind.data.frame(Expinf_AllTrtss, Attribute=names(Expinf_AllTrts$step1))
  colnames(Expinf_AllTrtss)[1]<-"Expected_infl"
  
  Expinf_AllTrtss$iter<-i
  pcorr_Alltrts$iter<-i
  
  corrAlltrts_out[[i]]<-rbind.data.frame(pcorr_Alltrts)
  ExpInf_Alltrts_out[[i]]<-rbind.data.frame(Expinf_AllTrtss)
}

CorrAllTrts_iters<-bind_rows(corrAlltrts_out)
ExpInfAllTrts_iters<-bind_rows(ExpInf_Alltrts_out)

save(Corr_iters,ExpInf_iters,CorrAllTrts_iters,ExpInfAllTrts_iters, file="Data/TradCorr_Bootstrap.RData")