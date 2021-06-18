#---------------------------------------#
#  Partial correlation network analysis #
#---------------------------------------#
# loosely based on:
# Epskamp, S. & Fried, E.I. (2018) A Tutorial on Regularized Partial Correlation Networks. Psychological Methods, 23, 617–634.

require(tidyverse)
require(MVN)
require(bootnet)
require(networktools)
require(igraph)
require(ggraph)
require(tidygraph)
require(huge)

load(file="Data/median_RecoveryAttributes.RData")

#-------------------------#
# Multi-variate normal    #
# transform               #
#-------------------------#

# 7 trts (74 sites)
result_trn <- mvn(data = Attributes_20y_7trts[,3:9], mvnTest = "hz", scale=TRUE)

Str_trn<-huge.npn(Attributes_20y_7trts[,3:9], npn.func="truncation",verbose=T)

result_trnn <- mvn(data = Str_trn[,1:7], mvnTest = "hz", scale=TRUE)

# 8 trts (43 sites)
result_trn2 <- mvn(data = Attributes_20y_8trts[,3:10], mvnTest = "hz", scale=TRUE)

Str_trn_2<-huge.npn(Attributes_20y_8trts[,3:10], npn.func="truncation",verbose=T)

result_trnn2 <- mvn(data = Str_trn_2[,1:8], mvnTest = "hz", scale=TRUE)

#---------------------------------------#
# Estimate partial correlation networks #
#---------------------------------------#

net20y_7trts<-estimateNetwork(Str_trn, default="EBICglasso",
                            corMethod="cor_auto", weighted = T, signed=T,directed=F,
                            tuning=0.5, missing="pairwise", 
                            threshold=T, criterion="ebic", verbose=F)

net20y_8trts<-estimateNetwork(Str_trn_2, default="EBICglasso",
                            corMethod="cor_auto", weighted = T, signed=T,directed=F,
                            tuning=0.5, missing="pairwise",
                            threshold=T, criterion="ebic", verbose=F)

#-----------------#  
# bootstrap       #
#-----------------#

comm<-c("Structure","Structure","Structure","Diversity","Diversity","Traits","Traits")

boot_20y7trt_stab<-bootnet(net20y_7trts, statistics=c("strength","edge","expectedInfluence","bridgeStrength","bridgeCloseness","bridgeBetweenness","bridgeExpectedInfluence"),
                           weighted=T,signed=T, directed=F,computeCentrality=T, nCores=2,nBoots=10000,communities=comm, type="case")

corStability(boot_20y7trt_stab) # edge, expectedInfluence, and strength ok (correlation stabiilty coefficients >0.4)

# edge = 0.432
# expected influence = 0.514 

boot_20y7trt<-bootnet(net20y_7trts, statistics=c("strength","edge","expectedInfluence","bridgeStrength","bridgeCloseness","bridgeBetweenness","bridgeExpectedInfluence"),
                      weighted=T,signed=T, directed=F,computeCentrality=T, 
                      caseMin = 0.365, caseMax = 0.595,
                      nCores=2,nBoots=10000,communities=comm)

comm<-c("Structure","Structure","Structure","Diversity","Diversity","Diversity","Traits","Traits")

boot_20y8trt_stab<-bootnet(net20y_8trts, statistics=c("strength","edge","expectedInfluence","bridgeStrength","bridgeCloseness","bridgeBetweenness","bridgeExpectedInfluence"),
                           weighted=T,signed=T, directed=F,computeCentrality=T, nCores=2,nBoots=10000,communities=comm, type="case")

corStability(boot_20y8trt_stab) # no metric had correlation stabiilty coefficients >0.4)
# edge=  0.19
# expected influence = 0.13 

boot_20y8trt<-bootnet(net20y_8trts, statistics=c("strength","edge","expectedInfluence","bridgeStrength","bridgeCloseness","bridgeBetweenness","bridgeExpectedInfluence"),
                      weighted=T,signed=T, directed=F,computseCentrality=T, nCores=2,
                      caseMin = 0.13, caseMax = 0.283,
                      nBoots=10000,communities=comm)

# 
save(boot_20y7trt,boot_20y7trt_stab, boot_20y8trt,boot_20y8trt_stab,
     file="Data/Networks_Metrics_RecovRates_20y.RData")