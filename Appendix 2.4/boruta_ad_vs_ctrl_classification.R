library(plyr)
library(tidyverse)
library(caret)
library(Boruta)

# To run code using multiple threads
#######################################################################################################
#library(parallel)
#library(doParallel)

#cluster <- makeCluster(20)
#registerDoParallel(cluster)
#######################################################################################################
set.seed(123)

path <- "/home/user/data"

setwd(path)


# input data set up - the labels and data files are produced by the "prepare_for_m_learn" script

labels_ctrl_vs_ad <- readRDS("labels_ctrl_vs_ad_ml_class.rds")
data_ctrl_vs_ad <- readRDS("data_ctrl_vs_ad_ml_class.rds")



#######################################################################################################


#CONTROL vs AD
data_labels_ctrl_vs_ad <- cbind(data_ctrl_vs_ad, labels_ctrl_vs_ad)
boruta_ctrl_vs_ad_model=Boruta(labels_ctrl_vs_ad~. , data=data_labels_ctrl_vs_ad, doTrace = 2 , maxRuns = 10000)
saveRDS(boruta_ctrl_vs_ad_model,file="boruta_ctrl_vs_ad_model.rds")
final_boruta_ctrl_vs_ad_model=TentativeRoughFix(boruta_ctrl_vs_ad_model)
ctrl_vs_ad_model_sel_feat=getSelectedAttributes(final_boruta_ctrl_vs_ad_model)
saveRDS(ctrl_vs_ad_model_sel_feat,file="ctrl_vs_ad_model_boruta_sel_feat.rds")

k <-lapply(1:ncol(boruta_ctrl_vs_ad_model$ImpHistory),function(i)
  boruta_ctrl_vs_ad_model$ImpHistory[is.finite(boruta_ctrl_vs_ad_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_ctrl_vs_ad_model$ImpHistory)
boruta_ctrl_vs_ad_model_feature_rank <- sort(sapply(k,median))
write.csv(boruta_ctrl_vs_ad_model_feature_rank,file="boruta_ctrl_vs_ad_model_feature_rank.csv")

