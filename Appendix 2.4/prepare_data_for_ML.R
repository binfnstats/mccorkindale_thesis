library(edgeR)
library(sva)

setwd("/home/user/data")

# read in covar file and reads file
covar <- read.delim("full_mRNA_cohort.csv", sep = ",", header = TRUE,row.names = 1)
reads <- read.delim("rsem_gene_counts_638.tsv",row.names = 1,sep="\t", header = TRUE, stringsAsFactors = FALSE)


# order covar and read files so that they match
covar=covar[order(covar$sampleID), ]
reads1 <- reads[,order(names(reads))]

# read in list of cases to exclude then remove them from covar and reads
remove_list_ml <- readRDS("ml_exclusions.rds")
reads1 <- reads1[,!(names(reads1) %in% remove_list_ml)]
covar <- covar[!covar$sampleID %in% remove_list_ml, ]
# rename read files so that diagnosis is included in sample names (AD_sampleID,CTRL_sampleID)
names(reads1) <- covar$Dx_sampleID

# check for columns with abnormally low counts
dim(reads1)
reads1=reads1[,which(!colSums(reads1)<200000)] 
dim(reads1)
dim(covar)

# filter the reads to only keep genes that have a total sum of 2500 across all subjects (roughly an average of 5)
reads_fil <- reads_wo_ctrl[which(!rowSums(reads_wo_ctrl)<=2500),] 

dim(reads_fil)
# normalise the reads using the TMM method using edgeR
y=DGEList(counts=reads_fil)
y=calcNormFactors(y,method = "TMM")
reads_norm<-cpm(y, normalized.lib.sizes=T)


# log transform the reads to make the counts closer to a normal distribution
dat1 <- log2(reads_norm+1)

# apply ComBat batch correction to output batch corrected counts (input for machine learning)
combat <- ComBat(dat=dat1,batch=covar$batch,mod=NULL,par.prior=TRUE,prior.plots = FALSE)
dat <- t(data.frame(combat))



# save gene set and covar from full ml cohort

saveRDS(dat, "gene_input_full_ml_cohort.rds")
saveRDS(covar, "covar_full_ml_cohort.rds")



# prepare label and data objects for ml classification and regression
# Full cohort used for ml regression. Subset with equal numbers of AD and control used for classification

#select covariate columns for machine learning
vars <- c("Dx_sampleID","cogdx","episodic_memory","global_cognition","fast_decline_episodic","fast_decline_global")
ml_covar <- covar[,names(covar) %in% vars]

# confirm covar and reads are in the same order
covar_names <- covar$Dx_sampleID
reads_names <- rownames(dat)
identical(reads_names, covar_names)

# combine genes and covariates into one object
covar_genes <- cbind(dat,ml_covar)

# save full cohort with covariate columns for machine learning regression
saveRDS(covar_genes, "ml_regression_input_full_cohort.rds")

# balanced sub cohort for machine learning classification

# reag1_reag0
set.seed(123)
mldata <- covar_genes

#mldata$keep[mldata$cogdx==1 | mldata$cogdx==2] <- 1
#mldata <- mldata[complete.cases(mldata),]
mldata$labels[mldata$cogdx==1] <- "NCI"
mldata$labels[mldata$cogdx==4] <- "DEM"
table(mldata$cogdx)

cogdx1 <- mldata[mldata$cogdx==1,]
cogdx4 <- mldata[mldata$cogdx==4,]
cogdx4_172 <- cogdx1[sample(nrow(cogdx1), 172), ]

mldata <- rbind(cogdx1,cogdx4_172) 

mldata$labels[mldata$cogdx==1] <- "NCI"
mldata$labels[mldata$cogdx==4] <- "DEM"

mldata <- mldata[complete.cases(mldata),]
mldata <- droplevels(mldata)
contrast_labels <- mldata$labels

rm_cols <- vars <- c("Dx_sampleID","cogdx","episodic_memory","global_cognition","fast_decline_episodic","fast_decline_global")

mldata <- mldata[,!(names(mldata) %in% rm_cols)]


saveRDS(mldata,"data_dem_vs_nci_ml_class.rds")
saveRDS(contrast_labels,"labels_dem_vs_nci_ml_class.rds")