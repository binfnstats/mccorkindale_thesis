# msbb data prep

# necessary for machine learning set up
library(edgeR)
library(sva)

setwd("/Users/andrew/Dropbox\ (Sydney\ Uni)/Harari_MSBB/MSBB/")

# set up
covar_clin <- read.csv("2020_10_29_MSBB_clinical.csv", header = TRUE)
covar_tech <- read.csv("2020_10_29_MSBB_technical.csv", header = TRUE)

# # read in the updated file that includes the plaque mean column
covar_plaq <- read.csv("MSBB_individual_metadata_trimmed.csv", header = TRUE)

# #confirm ID columns are the same for the covar files and then combine
identical(covar_clin$individualID,covar_tech$individualID)

# 
covar <- cbind(covar_clin,covar_tech)
covar$BrodmannArea <- as.factor(covar$Tissue)
covar$batch <- as.factor(covar$batch)
covar$sex <- as.factor(covar$sex)
covar$ethnicity  <- as.factor(covar$ethnicity )
covar$Action <- as.factor(covar$Action)
covar$apoeGenotype <- as.factor(covar$apoeGenotype)
covar$pmi_hrs <- (covar$pmi)/60
# 
# # make diagnosis columns - note that coding is different to ROSMAP
covar$path_AD[covar$CERAD==3] <- "AD"
covar$path_AD[covar$CERAD==2] <- "AD"
covar$path_AD[covar$CERAD==1] <- "CTRL"
covar$path_AD[covar$CERAD==4] <- "CTRL"
covar$path_AD <- as.factor(covar$path_AD)
# 
# 

covar$cogdx[covar$CDR==0] <- "NCI" 
covar$cogdx[covar$CDR==0.5] <- "MCI" 
covar$cogdx[covar$CDR==1] <- "DEM"
covar$cogdx[covar$CDR==2] <- "DEM"
covar$cogdx[covar$CDR==3] <- "DEM"
covar$cogdx[covar$CDR==4] <- "DEM"
covar$cogdx[covar$CDR==5] <- "DEM"
covar$cogdx <- as.factor(covar$cogdx)

# binary cognitive diagnosis
covar$demNoDem[covar$CDR==0] <- "noDEM"
covar$demNoDem[covar$CDR==0.5] <- "noDEM"
covar$demNoDem[covar$CDR==1] <- "DEM"
covar$demNoDem[covar$CDR==2] <- "DEM"
covar$demNoDem[covar$CDR==3] <- "DEM"
covar$demNoDem[covar$CDR==4] <- "DEM"
covar$demNoDem[covar$CDR==5] <- "DEM"
covar$demNoDem <- as.factor(covar$demNoDem)
# 
# 
## mark samples to be removed based on quality as in previous studies
# covar$remove[covar$RIN<4] <- 1
# covar$remove[covar$Action=="Exclude"] <- 1
# covar$remove[covar$Action=="Remap"] <- 1


# make a covariate file for each brain region
# covar_bm10 <- covar[covar$BrodmannArea=="BM10",]
# covar_bm22 <- covar[covar$BrodmannArea=="BM22",]
# covar_bm36 <- covar[covar$BrodmannArea=="BM36",]
# covar_bm44 <- covar[covar$BrodmannArea=="BM44",]
# 
# write.csv(covar_bm10,"covar_bm10.csv")
# write.csv(covar_bm22,"covar_bm22.csv")
# write.csv(covar_bm36,"covar_bm36.csv")
# write.csv(covar_bm44,"covar_bm44.csv")
# 

#######################################################################################################
## BM10
#######################################################################################################
setwd("/inDir/")
reads_bm10 <- read.delim("BM10/gene_count_matrix.tsv",sep ="\t" , header = TRUE,row.names = 1)

covar_bm10 <- read.csv("covar_bm10.csv")
# remove the cases that meet the exclusion criteria
covar_bm10 <- covar_bm10[!covar_bm10$remove==1,]


# remove samples based on PCA and hierarchical clustering
rmlist_bm10 <- c("hB_RNA_12744","hB_RNA_13477","hB_RNA_13330","hB_RNA_13373_L43C014","hB_RNA_12680")
reads_bm10 <- reads_bm10[,!(names(reads_bm10) %in% rmlist_bm10 )]
covar_bm10=covar_bm10[!(covar_bm10$sampleIdentifier %in% rmlist_bm10 ), ]

reads_bm10 <- reads_bm10[,order(names(reads_bm10))]
covar_bm10 <- covar_bm10[order(covar_bm10$sampleIdentifier), ]

# PCA identified two distinct batches
covar_bm10$batch_adj[covar_bm10$batch=="S107B355"] <- "B"
covar_bm10$batch_adj[covar_bm10$batch=="S108B355"] <- "B"
covar_bm10$batch_adj[covar_bm10$batch=="S109B355"] <- "B"
covar_bm10$batch_adj[covar_bm10$batch=="S111B394"] <- "B"
covar_bm10$batch_adj[covar_bm10$batch=="S112B394"] <- "B"

covar_bm10$batch_adj[covar_bm10$batch=="E007C014"] <- "C"
covar_bm10$batch_adj[covar_bm10$batch=="K75rC014"] <- "C"
covar_bm10$batch_adj[covar_bm10$batch=="L43C014"] <- "C"
covar_bm10$batch_adj[covar_bm10$batch=="P60C014"] <- "C"
covar_bm10$batch_adj[covar_bm10$batch=="P62C014"] <- "C"
covar_bm10$batch_adj <- as.factor(covar_bm10$batch_adj)

# use as input to diffex etc
reads_wo_ctrl=reads_bm10[,which(!colSums(reads_bm10)<200000)] 
covar <- covar_bm10


#######################################################################################################
## BM22
#######################################################################################################
setwd("/Users/andrew/Dropbox\ (Sydney\ Uni)/Harari_MSBB/MSBB/")
reads_bm22 <- read.delim("BM22/gene_count_matrix.tsv",sep ="\t" , header = TRUE,row.names = 1)

covar_bm22 <- read.csv("covar_bm22.csv")

covar_bm22 <- covar_bm22[!covar_bm22$remove==1,]

reads_bm22 <- reads_bm22[,names(reads_bm22) %in% covar_bm22$sampleIdentifier]

# remove samples based on PCA and hierarchical clustering
rmlist_bm22 <- c("hB_RNA_7925","hB_RNA_7955","hB_RNA_7995","hB_RNA_8015","hB_RNA_8115","hB_RNA_9190","hB_RNA_9193","hB_RNA_9199")

reads_bm22 <- reads_bm22[,!(names(reads_bm22) %in% rmlist_bm22 )]
covar_bm22 <- covar_bm22[!(covar_bm22$sampleIdentifier %in% rmlist_bm22 ), ]

reads_bm22 <- reads_bm22[,order(names(reads_bm22))]
covar_bm22 <- covar_bm22[order(covar_bm22$sampleIdentifier), ]

# recode batch variable because the batches clearly clustered into two groups
# based on if the batch ended with B394/5 or C014

covar_bm22$batch_adj[covar_bm22$batch=="H154B394"] <- "B"
covar_bm22$batch_adj[covar_bm22$batch=="S112B394"] <- "B"
covar_bm22$batch_adj[covar_bm22$batch=="S113B355"] <- "B"
covar_bm22$batch_adj[covar_bm22$batch=="S113B394"] <- "B"
covar_bm22$batch_adj[covar_bm22$batch=="S114B394"] <- "B"
covar_bm22$batch_adj[covar_bm22$batch=="S115B394"] <- "B"

covar_bm22$batch_adj[covar_bm22$batch=="E007C014"] <- "C"
covar_bm22$batch_adj[covar_bm22$batch=="E2C014"] <- "C"
covar_bm22$batch_adj[covar_bm22$batch=="E3C014"] <- "C"
covar_bm22$batch_adj[covar_bm22$batch=="K75rC014"] <- "C"
covar_bm22$batch_adj[covar_bm22$batch=="K79C014"] <- "C"
covar_bm22$batch_adj[covar_bm22$batch=="L43C014"] <- "C"
covar_bm22$batch_adj <- as.factor(covar_bm22$batch_adj)

# use as input to diffex etc
reads_wo_ctrl=reads_bm22[,which(!colSums(reads_bm22)<200000)] 
covar <- covar_bm22
setwd("/Users/andrew/Dropbox\ (Sydney\ Uni)/Harari_MSBB/MSBB/BM22/STAR/results/")

#######################################################################################################
## BM36
#######################################################################################################
setwd("/Users/andrew/Dropbox\ (Sydney\ Uni)/Harari_MSBB/MSBB/")
reads_bm36 <- read.delim("BM36/STAR/genecounts/gene_count_matrix.tsv",sep ="\t" , header = TRUE,row.names = 1)

covar_bm36 <- read.csv("covar_bm36.csv")
covar_bm36 <- covar_bm36[!covar_bm36$remove==1,]

reads_bm36  <- reads_bm36[,names(reads_bm36) %in% covar_bm36$sampleIdentifier]

# remove samples based on PCA and hierarchical clustering or duplicates
rmlist_bm36 <- c("hB_RNA_10372","hB_RNA_10822_resequenced","hB_RNA_10912","hB_RNA_10972",
                 "hB_RNA_11023","hB_RNA_11042","hB_RNA_11072","hB_RNA_12222")
reads_bm36 <- reads_bm36[,!(names(reads_bm36) %in% rmlist_bm36 )]
covar_bm36=covar_bm36[!(covar_bm36$sampleIdentifier %in% rmlist_bm36 ), ]

reads_bm36 <- reads_bm36[,order(names(reads_bm36))]
covar_bm36 <- covar_bm36[order(covar_bm36$sampleIdentifier), ]

covar_bm36$batch_adj[covar_bm36$batch=="P19B648"] <- "B"
covar_bm36$batch_adj[covar_bm36$batch=="S109B355"] <- "B"
covar_bm36$batch_adj[covar_bm36$batch=="S110B355"] <- "B"
covar_bm36$batch_adj[covar_bm36$batch=="S111B355"] <- "B"
covar_bm36$batch_adj[covar_bm36$batch=="S115B355"] <- "B"
covar_bm36$batch_adj[covar_bm36$batch=="S151B648"] <- "B"

covar_bm36$batch_adj[covar_bm36$batch=="E007C014"] <- "C"
covar_bm36$batch_adj[covar_bm36$batch=="K76C014"] <- "C"
covar_bm36$batch_adj[covar_bm36$batch=="K77C014"] <- "C"
covar_bm36$batch_adj[covar_bm36$batch=="K79C014"] <- "C"
covar_bm36$batch_adj[covar_bm36$batch=="L43C014"] <- "C"
covar_bm36$batch_adj <- as.factor(covar_bm36$batch_adj)

# use as input to diffex etc
reads_wo_ctrl=reads_bm36[,which(!colSums(reads_bm36)<200000)] 
covar <- covar_bm36


#######################################################################################################
## BM44
#######################################################################################################
setwd("/Users/andrew/Dropbox\ (Sydney\ Uni)/Harari_MSBB/MSBB/")
reads_bm44 <- read.delim("BM44/STAR/genecounts/gene_count_matrix.tsv",sep ="\t" , header = TRUE,row.names = 1)

covar_bm44 <- read.csv("covar_bm44.csv")
covar_bm44 <- covar_bm44[!covar_bm44$remove==1,]

reads_bm44 <- reads_bm44[,names(reads_bm44) %in% covar_bm44$sampleIdentifier]

# remove samples based on PCA and hierarchical clustering
rmlist_bm44 <- c("hB_RNA_16555","hB_RNA_16735","hB_RNA_16775","hB_RNA_16895","hB_RNA_16905","hB_RNA_4398",
                 "hB_RNA_4720_L43C014","hB_RNA_4782","hB_RNA_4862_L43C014","hB_RNA_4881","hB_RNA_4891")
reads_bm44 <- reads_bm44[,!(names(reads_bm44) %in% rmlist_bm44 )]
covar_bm44=covar_bm44[!(covar_bm44$sampleIdentifier %in% rmlist_bm44 ), ]

reads_bm44 <- reads_bm44[,order(names(reads_bm44))]
covar_bm44 <- covar_bm44[order(covar_bm44$sampleIdentifier), ]

covar_bm44$batch_adj[covar_bm44$batch=="S107B355"] <- "B"
covar_bm44$batch_adj[covar_bm44$batch=="S108B355"] <- "B"
covar_bm44$batch_adj[covar_bm44$batch=="S109B355"] <- "B"
covar_bm44$batch_adj[covar_bm44$batch=="S111B394"] <- "B"
covar_bm44$batch_adj[covar_bm44$batch=="S112B394"] <- "B"

covar_bm44$batch_adj[covar_bm44$batch=="E007C014"] <- "C"
covar_bm44$batch_adj[covar_bm44$batch=="K75rC014"] <- "C"
covar_bm44$batch_adj[covar_bm44$batch=="L43C014"] <- "C"
covar_bm44$batch_adj[covar_bm44$batch=="P60C014"] <- "C"
covar_bm44$batch_adj[covar_bm44$batch=="P62C014"] <- "C"
covar_bm44$batch_adj <- as.factor(covar_bm44$batch_adj)

#######################################################################################################
## filtering and normalisation done separately in each region - example is for BM10
#######################################################################################################

setwd("/BM10")


reads_wo_ctrl=reads_bm10[,which(!colSums(reads_bm10)<200000)] 
covar <- covar_bm10

# filter based on total count per row = roughly an average of 5 across all samples
reads_fil <- reads_wo_ctrl[which(!rowSums(reads_wo_ctrl)<=1500),]
#reads_fil <- reads_wo_ctrl[which(!rowSums(reads_wo_ctrl)<=2000),]
#reads_fil <- reads_wo_ctrl[which(!rowSums(reads_wo_ctrl)<=2500),]
dim(reads_fil)

# create differential expression object in edgeR
y=DGEList(counts=reads_fil)

# normalise in edgeR using TMM method
y=calcNormFactors(y,method = "TMM")

# output a table with the normalised reads for other functions if desired
reads_norm<-cpm(y, normalized.lib.sizes=T)
datExpr0=t(reads_norm)
write.csv(reads_norm, "bm10_normalised_reads_1500.csv")

# for individual regional analysis - batch correction was done
# for the combined analysis - counts were log-transformed and scaled before being
# combined and then batch correction performed

data_BM10 <- data.frame(scale(log2(reads_norm + 1)))


#######################################################################################################
# Testing out combining MSBB regions with log2 transformation
#######################################################################################################

# all final covariate files were combined and saved
ldat_BM_ALL <- readRDS("covar_combined.RDS")


# only keep  common genes in MSBB
keep_list <- Reduce(intersect,list(names(data_BM10),names(data_BM22),names(data_BM36),names(data_BM44)))
data_BM10 <- data_BM10[,names(data_BM10) %in% keep_list]
data_BM22 <- data_BM22[,names(data_BM22) %in% keep_list]
data_BM36 <- data_BM36[,names(data_BM36) %in% keep_list]
data_BM44 <- data_BM44[,names(data_BM44) %in% keep_list]

data_MSBB_ALL <-rbind(data_BM10,data_BM22,data_BM36,data_BM44)

data_MSBB_ALL <- t(data_MSBB_ALL)
combat <- ComBat(dat=data_MSBB_ALL,batch=ldat_BM_ALL$batch_adj,mod=NULL,par.prior=TRUE,prior.plots = FALSE)
data_MSBB_comb <- data.frame(t(combat))

# output is batch corrected counts
# PCA plots showed that ComBat removed the moderate batch effect and lack of regional effect


# to run an analysis within a brain region then use the following to extract
# 
ldat_bm10 <- ldat_BM_ALL[ldat_BM_ALL$BrodmannArea=="BM10",]
ldat_bm22 <- ldat_BM_ALL[ldat_BM_ALL$BrodmannArea=="BM22",]
ldat_bm36 <- ldat_BM_ALL[ldat_BM_ALL$BrodmannArea=="BM36",]
ldat_bm44 <- ldat_BM_ALL[ldat_BM_ALL$BrodmannArea=="BM44",]

# #######################################################################################################
# # BM10
# #######################################################################################################
# 
data_bm10 <- data_MSBB_comb[rownames(data_MSBB_comb) %in% ldat_bm10$sampleIdentifier,]

labels_bm10 <- ldat_bm10$path_AD

# run boruta
boruta_bm10_model=Boruta(x=data_bm10 , y=labels_bm10, doTrace = 2 , maxRuns = 10000,num.threads=35)
saveRDS(boruta_bm10_model,file="boruta_bm10_model_pathAD_combat.rds")
final_boruta_bm10_model=TentativeRoughFix(boruta_bm10_model)
bm10_model_sel_feat=getSelectedAttributes(final_boruta_bm10_model)
saveRDS(bm10_model_sel_feat,file="bm10_model_boruta_sel_feat_pathAD_combat.rds")

k <-lapply(1:ncol(boruta_bm10_model$ImpHistory),function(i)
  boruta_bm10_model$ImpHistory[is.finite(boruta_bm10_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_bm10_model$ImpHistory)
boruta_bm10_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
boruta_bm10_model_rank_out <- data.frame(features=names(boruta_bm10_model_feature_rank), z_score=boruta_bm10_model_feature_rank, row.names=NULL)

write.csv(boruta_bm10_model_rank_out,file="boruta_bm10_model_feature_pathAD.csv",row.names = FALSE)

# #######################################################################################################
# Regions combined
# #######################################################################################################
# 

data_MSBB_comb <- data_MSBB_comb[rownames(data_MSBB_comb) %in% ldat_BM_ALL$sampleIdentifier,]
labels_BM_ALL <- ldat_BM_ALL$path_AD

boruta_ROS_MSBB_model=Boruta(x=data_MSBB_comb , y=labels_BM_ALL, doTrace = 2 , maxRuns = 10000,num.threads=35)
saveRDS(boruta_ROS_MSBB_model,file="boruta_MSBBALL_model_Braak_pathAD.rds")
final_boruta_ROS_MSBB_model=TentativeRoughFix(boruta_ROS_MSBB_model)
ROS_MSBB_model_sel_feat=getSelectedAttributes(final_boruta_ROS_MSBB_model)
saveRDS(ROS_MSBB_model_sel_feat,file="MSBBALL_model_boruta_sel_feat_pathAD.rds")

k <-lapply(1:ncol(boruta_ROS_MSBB_model$ImpHistory),function(i)
  boruta_ROS_MSBB_model$ImpHistory[is.finite(boruta_ROS_MSBB_model$ImpHistory[,i]),i])
names(k) <- colnames(boruta_ROS_MSBB_model$ImpHistory)
boruta_ROS_MSBB_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
boruta_ROS_MSBB_model_rank_out <- data.frame(features=names(boruta_ROS_MSBB_model_feature_rank), z_score=boruta_ROS_MSBB_model_feature_rank, row.names=NULL)

write.csv(boruta_ROS_MSBB_model_rank_out,file="boruta_MSBBALL_model_feature_pathAD.csv",row.names = FALSE)

# #######################################################################################################
# Regions combined predict GWAS genes - same code as for ROSMAP
# #######################################################################################################
# 

library(plyr)
library(tidyverse)
library(Boruta)

set.seed(123)

path <- "/home/genes_pred_gwas"
setwd(path)

dat <- readRDS("../data_MSBB_comb.rds")

covar <- readRDS("../covar_combined.RDS")

# read in list of genes that you want to do
gene_lists <- read.csv("gwas_hits_int.csv")
gwas_list <- gene_lists$gwas_hits

# some genes from the reference are duplicates with PAR_Y on the end
# this removes them
dat <- dat[,!(grepl("PAR_Y", names(dat)))]
names(dat) <- gsub("\\..*","",names(dat))



#######################################################################################################
# input is normalised reads - genes in columns, subjects in rows
# make gene of interest the label and then remove it from the input data

for (ctrast in gwas_list) {
  
  label_int <- dat[,names(dat) %in% ctrast]
  dat_int <- dat[,!names(dat) %in% ctrast]
  #set max threads to appropriate amount
  boruta_int_model=Boruta(x=dat_int,y=label_int, doTrace = 2 , maxRuns = 1000, num.threads=25)
  
  saveRDS(boruta_int_model,file=paste("boruta_",ctrast,"_model_genes_predict_ALL.rds",sep = ""))
  final_boruta_int_model=TentativeRoughFix(boruta_int_model)
  grn_model_sel_feat=getSelectedAttributes(final_boruta_int_model)
  saveRDS(grn_model_sel_feat,file=paste(ctrast,"_model_boruta_sel_feat_genes_predict_ALL.rds",sep = ""))
  
  k <-lapply(1:ncol(boruta_int_model$ImpHistory),function(i)
    boruta_int_model$ImpHistory[is.finite(boruta_int_model$ImpHistory[,i]),i])
  names(k) <- colnames(boruta_int_model$ImpHistory)
  boruta_int_model_feature_rank <- sort(sapply(k,median),decreasing = TRUE)
  boruta_int_model_rank_out <- data.frame(features=names(boruta_int_model_feature_rank), z_score=boruta_int_model_feature_rank, row.names=NULL)
  
  write.csv(boruta_int_model_rank_out,file=paste(ctrast,"_model_feature_rank_predict_ALL.csv",sep = ""),row.names = FALSE)
  
}


