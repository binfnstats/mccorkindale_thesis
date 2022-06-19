library("VennDiagram")
library("RUVSeq")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("topGO")
library("sleuth")
library("biomaRt")
library(WGCNA)
library("ggbiplot")
library(reshape)
library(gplots)
library(ops)
library(calibrate)
library(biomaRt)
library(sva)
library(ggplot2)
library("corrplot")
library(gage)
library("pathview")

setwd("/home/user/data")

# read in covar file and reads file
covar <- read.delim("full_mRNA_cohort.csv", sep = ",", header = TRUE,row.names = 1)
reads <- read.delim("rsem_gene_counts_638.tsv",row.names = 1,sep="\t", header = TRUE, stringsAsFactors = FALSE)


# order covar and read files so that they match
covar=covar[order(covar$sampleID), ]
reads1 <- reads[,order(names(reads))]

# read in list of cases to exclude then remove them from covar and reads
remove_list_de <- readRDS("diffex_exclusions.rds")
reads1 <- reads1[,!(names(reads1) %in% remove_list_de)]
covar <- covar[!covar$sampleID %in% remove_list_de, ]
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

# output a table with the normalised reads for other functions - Used for violin plots in Figure 1D
reads_norm<-cpm(y, normalized.lib.sizes=T)
datExpr0=t(reads_norm)

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

###PCA analaysis - usually overlay it with batch to see if there is a batch effect
rsem.pca <- prcomp(t(reads_norm), scale = TRUE)

pdf(file = "PCA_rsem_batch.pdf")
ggbiplot(rsem.pca, choices = 1:2, obs.scale = 1, var.scale = 1,groups = covar$Batch, ellipse = TRUE, circle = TRUE, var.axes = F,labels=covar$N)
scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()


# ANOVA to find factors contributing to variation in the data
z=reads_fil
m=melt(z)
colnames(m) <- c("sample_ID","counts")



ROS_or_MAP <- rep(as.numeric(covar$studyn), each=nrow(z))
Cognitive_Status <- rep(as.numeric(covar$cogdx), each=nrow(z))
Age <- rep(as.numeric(covar$age_death), each=nrow(z))
Education <- rep(as.numeric(covar$educ), each=nrow(z))
Sex <- rep(as.numeric(covar$msex), each=nrow(z))
PMI <- rep(as.numeric(covar$pmi), each=nrow(z))
Batch <- rep(as.numeric(covar$batch), each=nrow(z))
RIN <- rep(as.numeric(covar$rin), each=nrow(z))
Diagnosis <- rep(as.numeric(covar$ad_reagan), each=nrow(z))
APOE_Status <- rep(as.numeric(covar$apoe_num), each=nrow(z))


matrix <- data.frame(m, PMI,ROS_or_MAP,Age,Batch,RIN,Diagnosis,APOE_Status,Sex,Education,Cognitive_Status)
fit1 <- lm(counts ~ Diagnosis + PMI+ROS_or_MAP+Age+Cognitive_Status+ Batch+RIN+APOE_Status+Sex+Education, data=matrix)
a <- anova(fit1)
nfac <- length(a[,1])-1
maxval = 100
pdf(file="Anovar.pdf")
nfac <- length(a[,1])-1
barplot(100*a$"Sum Sq"[1:nfac]/sum(a$"Sum Sq"[1:nfac]),names.arg=rownames(a[1:nfac,]),ylim=c(0,maxval),cex.names=0.6,las=3,space=0.3,width=1,xlim=c(1,20))
dev.off()

# create design matrix for differential expression analysis
design <- model.matrix(~0 + ad_reagan + rin + batch + study+ msex + age_death, data=covar)
mart=useMart(biomart = "ENSEMBL_MART_ENSEMBL",host="http://Dec2017.archive.ensembl.org", dataset="hsapiens_gene_ensembl")

# set up kegg pathway analysis
kg.hsa=kegg.gsets("hsa")
#kg.hsa=kegg.gsets(species="hsa")
# select signaling and metabolic pathways
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
data(egSymb) # added as it was returning an error
kegg.gs.sym <- lapply(kegg.gs, eg2sym)
colors=brewer.pal(8,"Set3")
x=covar$ad_reagan


#######################################################################################################
## ctrl vs ad 
#######################################################################################################


# produce differential expression object
y=estimateGLMCommonDisp(y,design,verbose=T)
y=estimateGLMTagwiseDisp(y,design)
fit=glmFit(y,design)
pdf(file="BCV_10in5_ctrl_vs_ad.pdf")
plotBCV(y)
dev.off()
saveRDS(y, "10in5_diffex_y_ctrl_vs_ad.rds")
saveRDS(fit, "10in5_diffex_fit_ctrl_vs_ad.rds")
#contrasts with DEGs between the different groups - in this case AD vs control

# reag0 vs reag1
ctrl_vs_ad=glmLRT(fit, contrast=c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0))

saveRDS(adpath_relto_ctrl,"glmLRT_ctrl_vs_ad.rds")

# this line tells how many DEGs you have
degsum <- summary(dt_adpath_relto_ctrl <- decideTestsDGE(adpath_relto_ctrl, p=0.05, adjust="BH"))
degsum


summary(dt_adpath_relto_ctrl <- decideTestsDGE(adpath_relto_ctrl, p=0.05, adjust="BH"))
isDE_adpath_relto_ctrl <- as.logical(dt_adpath_relto_ctrl)
DE_adpath_relto_ctrl <- rownames(y)[isDE_adpath_relto_ctrl]
adpath_relto_ctrl$table$FDR <- p.adjust(adpath_relto_ctrl$table$PValue, method="BH")
DEG_adpath_relto_ctrl <- adpath_relto_ctrl$table[adpath_relto_ctrl$table$FDR < 0.05,]
DEG_adpath_relto_ctrl=DEG_adpath_relto_ctrl[DEG_adpath_relto_ctrl$FDR > 0,]
DEG_adpath_relto_ctrl=DEG_adpath_relto_ctrl[DEG_adpath_relto_ctrl$PValue < 0.05,]
#DEG_adpath_relto_ctrl=DEG_adpath_relto_ctrl[(DEG_adpath_relto_ctrl$logFC<(-0.5))|(DEG_adpath_relto_ctrl$logFC>0.5),]

results_DEG_adpath_relto_ctrl <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene", "gene_biotype","transcript_biotype","description", "chromosome_name", "start_position", "end_position", "strand","go_id"), filters = "ensembl_gene_id", values = sub("\\..+","",rownames(DEG_adpath_relto_ctrl)),mart = mart)

idx_DEG_adpath_relto_ctrl=match(sub("\\..+","",rownames(DEG_adpath_relto_ctrl)),results_DEG_adpath_relto_ctrl$ensembl_gene_id)

DEG_adpath_relto_ctrl$ensembl_gene_id=results_DEG_adpath_relto_ctrl$ensembl_gene_id[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$ensembl_transcript_id=results_DEG_adpath_relto_ctrl$ensembl_transcript_id[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$hgnc_symbol=results_DEG_adpath_relto_ctrl$hgnc_symbol[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$entrezgene=results_DEG_adpath_relto_ctrl$entrezgene[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$gene_biotype=results_DEG_adpath_relto_ctrl$gene_biotype[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$transcript_biotype=results_DEG_adpath_relto_ctrl$transcript_biotype[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$description=results_DEG_adpath_relto_ctrl$description[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$chromosome_name=results_DEG_adpath_relto_ctrl$chromosome_name[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$start_position=results_DEG_adpath_relto_ctrl$start_position[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$end_position=results_DEG_adpath_relto_ctrl$end_position[idx_DEG_adpath_relto_ctrl]
DEG_adpath_relto_ctrl$strand=results_DEG_adpath_relto_ctrl$strand[idx_DEG_adpath_relto_ctrl]
####DEG_adpath_relto_ctrl$go_id=results_DEG_adpath_relto_ctrl$go_id[idx_DEG_adpath_relto_ctrl]

# produces table of differentially expressed genes
write.table(DEG_adpath_relto_ctrl, file="adpath_relto_ctrl_DEG.tsv", sep="\t")
top_adpath_relto_ctrl=rownames(DEG_adpath_relto_ctrl)
adpath_relto_ctrl_cpm <- cpm(y)[top_adpath_relto_ctrl, ]
z_score_adpath_relto_ctrl <- ((adpath_relto_ctrl_cpm - rowMeans(adpath_relto_ctrl_cpm))/apply(adpath_relto_ctrl_cpm,1,sd))

write.table(adpath_relto_ctrl_cpm, file="adpath_relto_ctrl_DEG_cpm.tsv", sep="\t", row.names =T)

unique_DEG_adpath_relto_ctrl=unique(DEG_adpath_relto_ctrl$entrezgene)


pdf("adpath_relto_ctrl_heatmap.pdf")
pheatmap(z_score_adpath_relto_ctrl, scale="row", cluster_cols=T, cluster_rows=T, show_rownames=F,color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100))
dev.off()

pdf("adpath_relto_ctrl_MA_plot.pdf")
plotSmear(adpath_relto_ctrl, de.tags=DE_adpath_relto_ctrl)
abline(h=c(-1,1), col="blue")
dev.off()

pdf("adpath_relto_ctrl_pvalues.pdf")
hist(adpath_relto_ctrl$table$PValue, breaks=seq(0, 1, 0.05))
dev.off()

pdf("adpath_relto_ctrl_lrt_RLE.pdf")
plotRLE(adpath_relto_ctrl$fitted.values, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

anno_bp <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl")
allGenes_bp <- unique(unlist(anno_bp))
uniqueDEG_adpath_relto_ctrl <- unique(results_DEG_adpath_relto_ctrl$ensembl_gene_id)
geneList_adpath_relto_ctrl <- factor(as.integer(allGenes_bp %in% uniqueDEG_adpath_relto_ctrl))
names(geneList_adpath_relto_ctrl) <- allGenes_bp

GOdata_bp_adpath_relto_ctrl <- new("topGOdata", ontology = "BP", allGenes = geneList_adpath_relto_ctrl, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")
resultFisher_bp_adpath_relto_ctrl <- runTest(GOdata_bp_adpath_relto_ctrl, algorithm = "classic", statistic = "fisher")
resultKS_bp_adpath_relto_ctrl <- runTest(GOdata_bp_adpath_relto_ctrl, algorithm = "classic", statistic = "ks")
resultKS.elim_bp_adpath_relto_ctrl <- runTest(GOdata_bp_adpath_relto_ctrl, algorithm = "elim", statistic = "ks")

pdf("adpath_relto_ctrl_gotree_BP.pdf")
showSigOfNodes(GOdata_bp_adpath_relto_ctrl, score(resultFisher_bp_adpath_relto_ctrl), firstSigNodes = 10, useInfo = "all")
dev.off()

allRes_BP_adpath_relto_ctrl <- GenTable(GOdata_bp_adpath_relto_ctrl, classicFisher = resultFisher_bp_adpath_relto_ctrl, classicKS = resultKS_bp_adpath_relto_ctrl, elimKS = resultKS.elim_bp_adpath_relto_ctrl, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
for (i in 1:nrow(allRes_BP_adpath_relto_ctrl)){ allRes_BP_adpath_relto_ctrl$hgnc_symbol[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$hgnc_symbol)[which(grepl(allRes_BP_adpath_relto_ctrl$GO.ID[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
for (i in 1:nrow(allRes_BP_adpath_relto_ctrl)){ allRes_BP_adpath_relto_ctrl$ensg[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$ensembl_gene_id)[which(grepl(allRes_BP_adpath_relto_ctrl$GO.ID[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
write.table(allRes_BP_adpath_relto_ctrl, file="adpath_relto_ctrl_toptree_BP.tsv", sep="\t", row.names =F)

anno_cc <- annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "ensembl")
allGenes_cc <- unique(unlist(anno_cc))
uniqueDEG_adpath_relto_ctrl <- unique(results_DEG_adpath_relto_ctrl$ensembl_gene_id)
geneList_adpath_relto_ctrl <- factor(as.integer(allGenes_cc %in% uniqueDEG_adpath_relto_ctrl))
names(geneList_adpath_relto_ctrl) <- allGenes_cc

GOdata_cc_adpath_relto_ctrl <- new("topGOdata", ontology = "CC", allGenes = geneList_adpath_relto_ctrl, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")
resultFisher_cc_adpath_relto_ctrl <- runTest(GOdata_cc_adpath_relto_ctrl, algorithm = "classic", statistic = "fisher")
resultKS_cc_adpath_relto_ctrl <- runTest(GOdata_cc_adpath_relto_ctrl, algorithm = "classic", statistic = "ks")
resultKS.elim_cc_adpath_relto_ctrl <- runTest(GOdata_cc_adpath_relto_ctrl, algorithm = "elim", statistic = "ks")

pdf("adpath_relto_ctrl_gotree_CC.pdf")
showSigOfNodes(GOdata_cc_adpath_relto_ctrl, score(resultFisher_cc_adpath_relto_ctrl), firstSigNodes = 10, useInfo = "all")
dev.off()

allRes_cc_adpath_relto_ctrl <- GenTable(GOdata_cc_adpath_relto_ctrl, classicFisher = resultFisher_cc_adpath_relto_ctrl, classicKS = resultKS_cc_adpath_relto_ctrl, elimKS = resultKS.elim_cc_adpath_relto_ctrl, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
for (i in 1:nrow(allRes_cc_adpath_relto_ctrl)){ allRes_cc_adpath_relto_ctrl$hgnc_symbol[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$hgnc_symbol)[which(grepl(allRes_cc_adpath_relto_ctrl$GO.ID[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
for (i in 1:nrow(allRes_cc_adpath_relto_ctrl)){ allRes_cc_adpath_relto_ctrl$ensg[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$ensembl_gene_id)[which(grepl(allRes_cc_adpath_relto_ctrl$GO.ID[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
write.table(allRes_cc_adpath_relto_ctrl, file="adpath_relto_ctrl_toptree_CC.tsv", sep="\t", row.names =F)

anno_mf <- annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "ensembl")
allGenes_mf <- unique(unlist(anno_mf))
uniqueDEG_adpath_relto_ctrl <- unique(results_DEG_adpath_relto_ctrl$ensembl_gene_id)
geneList_adpath_relto_ctrl <- factor(as.integer(allGenes_mf %in% uniqueDEG_adpath_relto_ctrl))
names(geneList_adpath_relto_ctrl) <- allGenes_mf

GOdata_mf_adpath_relto_ctrl <- new("topGOdata", ontology = "MF", allGenes = geneList_adpath_relto_ctrl, nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ensembl")
resultFisher_mf_adpath_relto_ctrl <- runTest(GOdata_mf_adpath_relto_ctrl, algorithm = "classic", statistic = "fisher")
resultKS_mf_adpath_relto_ctrl <- runTest(GOdata_mf_adpath_relto_ctrl, algorithm = "classic", statistic = "ks")
resultKS.elim_mf_adpath_relto_ctrl <- runTest(GOdata_mf_adpath_relto_ctrl, algorithm = "elim", statistic = "ks")

pdf("adpath_relto_ctrl_gotree_MF.pdf")
showSigOfNodes(GOdata_mf_adpath_relto_ctrl, score(resultFisher_mf_adpath_relto_ctrl), firstSigNodes = 10, useInfo = "all")
dev.off()

allRes_mf_adpath_relto_ctrl <- GenTable(GOdata_mf_adpath_relto_ctrl, classicFisher = resultFisher_mf_adpath_relto_ctrl, classicKS = resultKS_mf_adpath_relto_ctrl, elimKS = resultKS.elim_mf_adpath_relto_ctrl, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
for (i in 1:nrow(allRes_mf_adpath_relto_ctrl)){ allRes_mf_adpath_relto_ctrl$hgnc_symbol[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$hgnc_symbol)[which(grepl(allRes_mf_adpath_relto_ctrl$GO.ID[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
for (i in 1:nrow(allRes_mf_adpath_relto_ctrl)){ allRes_mf_adpath_relto_ctrl$ensg[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$ensembl_gene_id)[which(grepl(allRes_mf_adpath_relto_ctrl$GO.ID[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
write.table(allRes_mf_adpath_relto_ctrl, file="adpath_relto_ctrl_toptree_MF.tsv", sep="\t", row.names =F)

unique_DEG_adpath_relto_ctrl=unique(DEG_adpath_relto_ctrl$entrezgene)
go_DEG_adpath_relto_ctrl=goana(unique_DEG_adpath_relto_ctrl)

top_BP_DEG_adpath_relto_ctrl=topGO(go_DEG_adpath_relto_ctrl,ont="BP",n=30)
for (i in 1:nrow(top_BP_DEG_adpath_relto_ctrl)){ top_BP_DEG_adpath_relto_ctrl$hgnc_symbol[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$hgnc_symbol)[which(grepl(row.names(top_BP_DEG_adpath_relto_ctrl)[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
for (i in 1:nrow(top_BP_DEG_adpath_relto_ctrl)){ top_BP_DEG_adpath_relto_ctrl$ensg[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$ensembl_gene_id)[which(grepl(row.names(top_BP_DEG_adpath_relto_ctrl)[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
write.table(top_BP_DEG_adpath_relto_ctrl, file="adpath_relto_ctrl_top1_BP.tsv", sep="\t")

top_MF_DEG_adpath_relto_ctrl=topGO(go_DEG_adpath_relto_ctrl,ont="MF",n=30)
for (i in 1:nrow(top_MF_DEG_adpath_relto_ctrl)){ top_MF_DEG_adpath_relto_ctrl$hgnc_symbol[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$hgnc_symbol)[which(grepl(row.names(top_MF_DEG_adpath_relto_ctrl)[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
for (i in 1:nrow(top_MF_DEG_adpath_relto_ctrl)){ top_MF_DEG_adpath_relto_ctrl$ensg[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$ensembl_gene_id)[which(grepl(row.names(top_MF_DEG_adpath_relto_ctrl)[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
write.table(top_MF_DEG_adpath_relto_ctrl, file="adpath_relto_ctrl_top1_MF.tsv", sep="\t")

top_CC_DEG_adpath_relto_ctrl=topGO(go_DEG_adpath_relto_ctrl,ont="CC",n=30)
for (i in 1:nrow(top_CC_DEG_adpath_relto_ctrl)){ top_CC_DEG_adpath_relto_ctrl$hgnc_symbol[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$hgnc_symbol)[which(grepl(row.names(top_CC_DEG_adpath_relto_ctrl)[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
for (i in 1:nrow(top_CC_DEG_adpath_relto_ctrl)){ top_CC_DEG_adpath_relto_ctrl$ensg[i]=paste(as.vector(results_DEG_adpath_relto_ctrl$ensembl_gene_id)[which(grepl(row.names(top_CC_DEG_adpath_relto_ctrl)[i],results_DEG_adpath_relto_ctrl$go_id))],collapse=",")}
write.table(top_CC_DEG_adpath_relto_ctrl, file="adpath_relto_ctrl_top1_CC.tsv", sep="\t")



all_ensg_mart <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", values = sub("\\..+","",rownames(adpath_relto_ctrl)),mart = mart)
all_ensg_mart= all_ensg_mart[!is.na(all_ensg_mart$entrezgene),]
all_ensg_mart = all_ensg_mart[!duplicated(all_ensg_mart$entrezgene), ]

idx_adpath_relto_ctrl=match(all_ensg_mart$ensembl_gene_id, sub("\\..+","",rownames(adpath_relto_ctrl$table)))
exp.fc=adpath_relto_ctrl$table$logFC[idx_adpath_relto_ctrl]
names(exp.fc)=all_ensg_mart$entrezgene[idx_adpath_relto_ctrl]
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
write.table(rbind(fc.kegg.p$greater, fc.kegg.p$less), file = "adpath_relto_ctrl_fc.kegg.p.tsv", sep = "\t")
fc.kegg.sig=sigGeneSet(fc.kegg.p)
write.table(rbind(fc.kegg.sig, fc.kegg.sig), file = "adpath_relto_ctrl.kegg.sig.tsv", sep ="\t")

fc.kegg.p.2p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL,same.dir = F)
write.table(rbind(fc.kegg.p.2p$greater, fc.kegg.p.2p$less), file = "adpath_relto_ctrl_fc.kegg.p.2p.tsv", sep = "\t")
fc.kegg.2p.sig=sigGeneSet(fc.kegg.p.2p)
write.table(rbind(fc.kegg.2p.sig, fc.kegg.2p.sig), file = "adpath_relto_ctrl.kegg.2d.sig.tsv", sep ="\t")


sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
pv.out.list <- sapply(path.ids2[1:20], function(pid) pathview( gene.data = exp.fc, pathway.id = pid, species = "hsa", out.suffix="adpath_relto_ctrl"))


sel <- fc.kegg.p.2p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p.2p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p.2p$greater)[sel]
sel.l <- fc.kegg.p.2p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p.2p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p.2p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
#pv.out.list <- sapply(path.ids2[1:20], function(pid) pathview( gene.data = exp.fc, pathway.id = pid, species = "hsa", out.suffix="2p.adpath_relto_ctrl"))

#adpath_relto_ctrl.kegg.esg.up <- esset.grp(adpath_relto_ctrl.kegg.p, adpath_relto_ctrl, gsets = kegg.gs, ref = hn, samp = dcis, test4up = T, output = T, outname = "adpath_relto_ctrl.kegg.up", make.plot = F)
#adpath_relto_ctrl.kegg.esg.dn <- esset.grp(adpath_relto_ctrl.kegg.p, adpath_relto_ctrl, gsets = kegg.gs, ref = hn, samp = dcis,test4up = F, output = T, outname = "adpath_relto_ctrl.kegg.dn", make.plot = F)
