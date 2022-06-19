library(plyr)
library(tidyverse)
library(caret)
library(Boruta)


set.seed(123)

path <- "/home/user/data"

setwd(path)


# input data set up - the labels and data files are produced by the "prepare_for_m_learn" script

labels_dem_vs_nci <- readRDS("labels_dem_vs_nci_ml_class.rds")
data_dem_vs_nci <- readRDS("data_dem_vs_nci_ml_class.rds")



#######################################################################################################



# Example code for running boruta for Dementia vs NCI classification
# the same code is used for classication and regression (just change the colum)


data_labels_dem_vs_nci <- cbind(data_dem_vs_nci, labels_dem_vs_nci)
# Boruta can be slow to run. Running with x = data and y = labels instead of the formula interface can
# lead to shorter run times
boruta_dem_vs_nci_model=Boruta(labels_dem_vs_nci~. , data=data_labels_dem_vs_nci, doTrace = 2 , maxRuns = 10000)
saveRDS(boruta_dem_vs_nci_model,file="boruta_dem_vs_nci_model.rds")
final_boruta_dem_vs_nci_model=TentativeRoughFix(boruta_dem_vs_nci_model)
dem_vs_nci_model_sel_feat=getSelectedAttributes(final_boruta_dem_vs_nci_model)
saveRDS(dem_vs_nci_model_sel_feat,file="dem_vs_nci_model_boruta_sel_feat.rds")

k <-lapply(1:ncol(boruta_dem_vs_nci_model$ImpHistory),function(i)
  boruta_dem_vs_nci_model$ImpHistory[is.finite(boruta_dem_vs_nci_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_dem_vs_nci_model$ImpHistory)
boruta_dem_vs_nci_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
boruta_dem_vs_nci_model_rank_out <- data.frame(features=names(boruta_dem_vs_nci_model_feature_rank), z_score=boruta_dem_vs_nci_model_feature_rank, row.names=NULL)

# outputs all genes ranked by z-score. Genes with z-scores above the ShadowMax value are considered "important"
write.csv(boruta_dem_vs_nci_model_rank_out,file="boruta_dem_vs_nci_model_feature_rank_predict.csv",row.names = FALSE)



#######################################################################################################

# Example of a regression analysis from the gene association analysis
# Here we were finding the others genes that were predictive of PRTN3  levels
# ##### prtn3 ENSG00000196415.9
# 
# 

# read in batch-corrected gene expression data 
data_prtn3 <- gene_data

#covar_prtn3 <- readRDS("prtn3_covar.rds")

# find the ENSG number for your gene of interest
labels_prtn3 <- data_prtn3$ENSG00000196415.9

rm_list <- ("ENSG00000196415.9")
data_prtn3 <- data_prtn3[,!(names(data_prtn3) %in% rm_list)]

#prtn3
data_labels_prtn3 <- cbind(data_prtn3, labels_prtn3)
boruta_prtn3_model=Boruta(labels_prtn3~. , data=data_labels_prtn3, doTrace = 2 , maxRuns = 10000)
saveRDS(boruta_prtn3_model,file="boruta_prtn3_model_genes_predict_2500.rds")
final_boruta_prtn3_model=TentativeRoughFix(boruta_prtn3_model)
prtn3_model_sel_feat=getSelectedAttributes(final_boruta_prtn3_model)
saveRDS(prtn3_model_sel_feat,file="prtn3_model_boruta_sel_feat_genes_predict_2500.rds")

k <-lapply(1:ncol(boruta_prtn3_model$ImpHistory),function(i)
  boruta_prtn3_model$ImpHistory[is.finite(boruta_prtn3_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_prtn3_model$ImpHistory)
boruta_prtn3_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
boruta_prtn3_model_rank_out <- data.frame(features=names(boruta_prtn3_model_feature_rank), z_score=boruta_prtn3_model_feature_rank, row.names=NULL)

write.csv(boruta_prtn3_model_rank_out,file="boruta_prtn3_model_feature_rank_predict_2500.csv",row.names = FALSE)

