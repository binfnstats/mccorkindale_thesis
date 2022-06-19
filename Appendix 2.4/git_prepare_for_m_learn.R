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

# filter the reads to only keep genes that have a count of at least 10 in 5 or more samples
keep_ad=apply(reads1[,which(grepl("AD",colnames(reads1)))],1,function(x) length(x[x>=10])>=5)
keep_ctrl=apply(reads1[,which(grepl("CTRL",colnames(reads1)))],1,function(x) length(x[x>=10])>=5) 

# subset reads based on filtered values - check dimensions of the data set
keep=keep_ctrl|keep_ad
reads_fil=reads1[keep,]

dim(reads1)
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
vars <- c("Dx_sampleID","ad_reagan","amyloid","amyloid_sf","tangles","tangles_sf")
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
mldata$labels[mldata$ad_reagan==0] <- "CTRL"
mldata$labels[mldata$ad_reagan==1] <- "AD"
table(mldata$ad_reagan)

ad_reagan0 <- mldata[mldata$ad_reagan==0,]
ad_reagan1 <- mldata[mldata$ad_reagan==1,]
adr1_189 <- ad_reagan1[sample(nrow(ad_reagan1), 189), ]

mldata <- rbind(ad_reagan0,adr1_189) 

mldata$labels[mldata$ad_reagan==0] <- "CTRL"
mldata$labels[mldata$ad_reagan==1] <- "AD"

mldata <- mldata[complete.cases(mldata),]
mldata <- droplevels(mldata)
contrast_labels <- mldata$labels

rm_cols <- vars <- c("Dx_sampleID","amyloid","amyloid_sf","tangles","tangles_sf")

mldata <- mldata[,!(names(mldata) %in% rm_cols)]


saveRDS(mldata,"data_ctrl_vs_ad_ml_class.rds")
saveRDS(contrast_labels,"labels_ctrl_vs_ad_ml_class.rds")