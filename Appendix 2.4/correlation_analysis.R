

#######################################################################################################
## APOE4 path correlation
#######################################################################################################

# read in data and covar file - code assumes that the sample names
# are in the same order across the two files
covar <- covariate_file
dat <- filtered_normalised_reads


# remove genes with PAR_Y in the name as they are duplicates and then 
# remove version number of ENSGs if desired
dat <- dat[,!(grepl("PAR_Y", names(dat)))]
names(dat) <- gsub("\\..*","",names(dat))

# make sure the inp
dat$apoe4d <- covar$apoe4d
dat_noE4 <- dat[dat$apoe4d==0,]
dat_E4 <- dat[dat$apoe4d==1,]
dat_E4$apoe4d <- NULL
dat_noE4$apoe4d <- NULL

covar_noE4 <- covar[covar$apoe4d==0,]
covar_E4 <- covar[covar$apoe4d==1,]

covar <- covar_noE4
covar <- covar_E4
covar <- droplevels(covar)
# create a dataframe of the variables of interest
#ROSMAP
dat_cor <- data.frame(covar$amyloid_mf,covar$tangles_mf,covar$cogn_ep_random_slope,
                      covar$np,as.numeric(covar$braaksc),as.numeric(covar$ceradsc))
names(dat_cor) <- c("amyloid_mf", "tangles_mf","cogn_ep_random_slope", "np","braaksc","ceradsc")
# MSBB
dat_cor <- data.frame(covar$plaquemean,covar$cerad_num,covar$Braak,covar$CDR)
names(dat_cor) <- c("plaquemean","cerad_num","Braak","CDR")

# prefer the spearman method as it is non-linear
dat <- dat_E4
#dat <- dat_noE4
gene_cor <- cor(x=dat,y =dat_cor,use='pairwise.complete.obs',method = "spearman")
gene_cor <- data.frame(gene_cor)

# make a column with the absolute correlation value for each variable
gene_cor$abs_amyloid_mf <- abs(gene_cor$covar.amyloid_mf)
gene_cor$abs_tangles_mf <- abs(gene_cor$covar.tangles_mf)
gene_cor$abs_cogn_ep_random_slope <- abs(gene_cor$covar.cogn_ep_random_slope)
gene_cor$abs_np <- abs(gene_cor$covar.np)
gene_cor$abs_braaksc <- abs(gene_cor$as.numeric.covar.braaksc.)
gene_cor$abs_ceradsc <- abs(gene_cor$as.numeric.covar.ceradsc.)

# add in gene information for the each transcript in "dat" using bioMart
# change host to desired version (best to specify as default version changes over time)
mart=useMart(biomart = "ENSEMBL_MART_ENSEMBL",host="https://Feb2021.archive.ensembl.org", dataset="hsapiens_gene_ensembl")

gl <- names(dat)
gt <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "entrezgene_id","gene_biotype",
                           "description"), filters = "ensembl_gene_id", values = gl,mart = mart)

idx=match(gl,gt$ensembl_gene_id)
ensembl_gene_id=gt$ensembl_gene_id[idx]
hgnc_symbol=gt$hgnc_symbol[idx]
entrezgene=gt$entrezgene[idx]
description=gt$description[idx]
gene_biotype=gt$gene_biotype[idx]

out <- cbind(ensembl_gene_id,hgnc_symbol,entrezgene,description,gene_biotype)
out <- data.frame(out)

dat <- cbind(out,gene_cor)


write.csv(gene_cor_info,"genes_cor_spearman_apoE4.csv")

# for chapter five. correlations were z-score transformed to compare between
# the ROSMAP and MSBB studies


#ROSMAP
dat$z_amyloid_mf <- (dat$amyloid_mf-mean(dat$amyloid_mf))/sd(dat$amyloid_mf)
dat$z_tangles_mf<- (dat$tangles_mf-mean(dat$tangles_mf))/sd(dat$tangles_mf)
dat$z_braaksc <- (dat$braaksc-mean(dat$braaksc))/sd(dat$braaksc)
dat$z_ceradsc <- (dat$ceradsc-mean(dat$ceradsc))/sd(dat$ceradsc)
dat$z_np <- (dat$np-mean(dat$np))/sd(dat$np)
dat$z_cogn_ep_random_slope <- (dat$cogn_ep_random_slope-mean(dat$cogn_ep_random_slope))/sd(dat$cogn_ep_random_slope)


#MSBB
dat$z_plaquemean <- (dat$plaquemean-mean(dat$plaquemean))/sd(dat$plaquemean)
dat$z_cerad_num <- (dat$cerad_num-mean(dat$cerad_num))/sd(dat$cerad_num)
dat$z_Braak <- (dat$Braak-mean(dat$Braak))/sd(dat$Braak)
dat$z_CDR <- (dat$CDR-mean(dat$CDR))/sd(dat$CDR)

