


########################################################################################################
# code to read in all boruta lists and sort them then find overlaps
#######################################################################################################

setwd("boruta_results_dir")
#
gene_lists <-read.csv("genes_of_int.csv")
gwas_list <- gene_lists$ensembl_gene_id
#gwas_list <- gwas_lists$ensembl_gene_id[c(1:10)]
for (ctrast in gwas_list) {
  # read in boruta results files for each gene
  cont_int <- read.csv(paste(ctrast,"model_feature_rank_predict.csv",sep = "_"),header = T)
  
  # select all with a z-score >= to 2 
  cont_int <- cont_int[cont_int$z_score>=2,]
  # alternatively you could select an arbitrary number e.g., the top 20
  #gl <- as.character(cont_int$features)[c(1:20)]
  # select just the gene lists
  gl <- as.character(cont_int$features)
  # extract only genes (remove shadoxMax/shadowMean/shadowMin)
  gl <- gl[grep("ENSG",gl)]
  #save the gene list in memory named by the ENSG
  assign(ctrast,gl)
  print(ctrast)
}


# generate the code required for the list below as all lists are named according to ENSG
# for (ctr in gwas_list) {
#   cat(paste(ctr,'=',ctr,',',sep=""))
# }

# 
gwas_top200 <- list(ENSG00000165029=ENSG00000165029,ENSG00000064687=ENSG00000064687,ENSG00000108798=ENSG00000108798,
                    ENSG00000159640=ENSG00000159640,ENSG00000137845=ENSG00000137845,ENSG00000154734=ENSG00000154734,
                    ENSG00000087116=ENSG00000087116,ENSG00000158859=ENSG00000158859,ENSG00000188157=ENSG00000188157,
                    ENSG00000138613=ENSG00000138613,ENSG00000130203=ENSG00000130203,ENSG00000142192=ENSG00000142192,
                    ENSG00000136717=ENSG00000136717,ENSG00000087589=ENSG00000087589,ENSG00000078487=ENSG00000078487,
                    ENSG00000108091=ENSG00000108091,ENSG00000198087=ENSG00000198087,ENSG00000105383=ENSG00000105383,
                    ENSG00000108556=ENSG00000108556,ENSG00000120885=ENSG00000120885,ENSG00000174469=ENSG00000174469,
                    ENSG00000203710=ENSG00000203710,ENSG00000134463=ENSG00000134463,ENSG00000229153=ENSG00000229153,
                    ENSG00000073712=ENSG00000073712,ENSG00000115641=ENSG00000115641,ENSG00000030582=ENSG00000030582,
                    ENSG00000135077=ENSG00000135077,ENSG00000163666=ENSG00000163666,ENSG00000196126=ENSG00000196126,
                    ENSG00000002587=ENSG00000002587,ENSG00000157368=ENSG00000157368,ENSG00000168918=ENSG00000168918,
                    ENSG00000174628=ENSG00000174628,ENSG00000103510=ENSG00000103510,ENSG00000131042=ENSG00000131042,
                    ENSG00000012223=ENSG00000012223,ENSG00000110079=ENSG00000110079,ENSG00000110077=ENSG00000110077,
                    ENSG00000071051=ENSG00000071051,ENSG00000142233=ENSG00000142233,ENSG00000166924=ENSG00000166924,
                    ENSG00000073921=ENSG00000073921,ENSG00000197943=ENSG00000197943,ENSG00000125534=ENSG00000125534,
                    ENSG00000196415=ENSG00000196415,ENSG00000120899=ENSG00000120899,ENSG00000029725=ENSG00000029725,
                    ENSG00000100599=ENSG00000100599,ENSG00000161929=ENSG00000161929,ENSG00000179526=ENSG00000179526,
                    ENSG00000140090=ENSG00000140090,ENSG00000137642=ENSG00000137642,ENSG00000066336=ENSG00000066336,
                    ENSG00000106460=ENSG00000106460,ENSG00000145901=ENSG00000145901,ENSG00000095970=ENSG00000095970,
                    ENSG00000265148=ENSG00000265148,ENSG00000148429=ENSG00000148429,ENSG00000186153=ENSG00000186153,
                    ENSG00000110514=ENSG00000110514)


gwas_top200 <- gwas_top200[order(names(gwas_top200))]

# make a list of the gene ontology results

# for (ctr in w_list) {
#   cat(paste(ctr,'=',ctr,',',sep=""))
# }

library(UpSetR)
library(tidyverse)

data_AD_C_cond = map(gwas_top200, ~ data.frame(target = .x, present = 1L,stringsAsFactors = F)) %>%
  bind_rows(.id = "source") %>%
  spread(source, present, 0) 

gl <- data_AD_C_cond$target
#gl <- names(dat)
gt <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "entrezgene_id","gene_biotype",
                           "description"), filters = "ensembl_gene_id", values = gl,mart = mart)

idx=match(gl,gt$ensembl_gene_id)
ensembl_gene_id=gt$ensembl_gene_id[idx]
hgnc_symbol=gt$hgnc_symbol[idx]
entrezgene=gt$entrezgene[idx]
description=gt$description[idx]
gene_biotype=gt$gene_biotype[idx]

#out <- cbind(ensembl_gene_id,hgnc_symbol,entrezgene,description,gene_biotype)
#out <- data.frame(out)

out_upset <- cbind(ensembl_gene_id,hgnc_symbol,entrezgene,description,gene_biotype)
out_upset <- data.frame(out_upset)

# use biomart to get gene information for the "target column and then join objects together
# change the column names to HGNC and then sort alphabetically

# ALL
names(data_AD_C_cond) <- c("target",gwas_lists$hgnc_symbol)
data_AD_C_cond <- data_AD_C_cond[,order(names(data_AD_C_cond))]
# numbers depend on the size of the list of input genes
# ROSMAP
# move "target" column to be the first column and then add a total column
data_AD_C_cond <- data_AD_C_cond[,c(55,1:54,56:62)]
data_AD_C_cond$Total <- rowSums(data_AD_C_cond[,c(2:62)])
# MSBB
data_AD_C_cond <- data_AD_C_cond[,c(53,1:52,54:60)]
data_AD_C_cond$Total <- rowSums(data_AD_C_cond[,c(2:60)])

gwas_upset_out <- cbind(out_upset,data_AD_C_cond)
gwas_upset_out <- gwas_upset_out[order(gwas_upset_out$Total,decreasing = TRUE),]


#write.csv(gwas_upset_out,"gwas_top400_Ros_All_ABS_upset.csv")
write.csv(gwas_upset_out,"gwas_boruta_MSBB_combat_upset.csv",row.names = FALSE)
#write.csv(gwas_upset_out,"top100_diseaseROS_MSBB_combat_upset.csv",row.names = FALSE)




#######################################################################################################
# Add gene information to Boruta output file and optionally run enrichment analysis
#######################################################################################################

library(biomaRt)
library(WebGestaltR)
# set the host to desired release. better to set the host for reproducibility
mart=useMart(biomart = "ENSEMBL_MART_ENSEMBL",host="https://Feb2021.archive.ensembl.org", dataset="hsapiens_gene_ensembl")

gwas_list <- list_of_vars_of_interest
for (ctrast in gwas_list) {
  # enter gene list
  setwd(path)
  # make sure the file name matches that of the Boruta output csv files
  cont_int <- read.csv(paste(ctrast,"model_feature_rank_predict.csv",sep = "_"),header = T)
  # extract ensembl IDs
  # can change the z_score threshold to something higher/lower
  cont_int <- cont_int[cont_int$z_score>=2,]
  cont_int <- cont_int[!cont_int$features=="shadowMax",]
  gl <- as.character(cont_int$features)
  gl <- gl[grep("ENSG",gl)]
  gl <- as.vector(gl)
  
  # use biomart to get gene information
  gt <- getBM(attributes = c("ensembl_gene_id", "refseq_mrna","hgnc_symbol", "entrezgene_id",
                             "description", "chromosome_name", "start_position", "end_position",
                             "strand","go_id"), filters = "ensembl_gene_id", values = gl,mart = mart)
  idx=match(gl,gt$ensembl_gene_id)
  
  #make a dataframe using the gene and z_score columns
  d <- cont_int
  d$ensembl_gene_id=gt$ensembl_gene_id[idx]
  d$ensembl_transcript_id=gt$ensembl_transcript_id[idx]
  d$hgnc_symbol=gt$hgnc_symbol[idx]
  d$entrezgene=gt$entrezgene[idx]
  d$gene_biotype=gt$gene_biotype[idx]
  d$transcript_biotype=gt$transcript_biotype[idx]
  d$description=gt$description[idx]
  d$chromosome_name=gt$chromosome_name[idx]
  d$start_position=gt$start_position[idx]
  d$end_position=gt$end_position[idx]
  d$strand=gt$strand[idx]
  
  # output the table containing the extra information from the mart columns
  write.csv(d,paste("boruta_gene",ctrast,"table.csv",sep="_"), row.names = FALSE)
  
  # run enrichment analysis using webGestalt on the gene list of interest
  # WebGestaltR(enrichMethod = "ORA",organism = "hsapiens",minNum = 10,maxNum=2000,enrichDatabase = "pathway_KEGG",
  #             interestGene = gl,interestGeneType = "ensembl_gene_id",referenceSet = "genome",fdrThr = 0.1,
  #             outputDirectory = getwd(),projectName = paste("kegg",ctrast,sep = "_"))
  # 
  # # GO BP
  # WebGestaltR(enrichMethod = "ORA",organism = "hsapiens",minNum = 10,maxNum=2000,enrichDatabase = "geneontology_Biological_Process",
  #             interestGene = gl,interestGeneType = "ensembl_gene_id",referenceSet = "genome",fdrThr = 0.1,
  #             outputDirectory = getwd(),projectName = paste("goBP",ctrast,sep = "_"))
  # 
  # # GO CC
  # WebGestaltR(enrichMethod = "ORA",organism = "hsapiens",minNum = 10,maxNum=2000,enrichDatabase = "geneontology_Cellular_Component",
  #             interestGene = gl,interestGeneType = "ensembl_gene_id",referenceSet = "genome",fdrThr = 0.1,
  #             outputDirectory = getwd(),projectName = paste("goCC",ctrast,sep = "_"))
  # 
  # 
  # # GO BPnoRedundant
  # WebGestaltR(enrichMethod = "ORA",organism = "hsapiens",minNum = 10,maxNum=2000,enrichDatabase = "geneontology_Biological_Process_noRedundant",
  #             interestGene = gl,interestGeneType = "ensembl_gene_id",referenceSet = "genome",fdrThr = 0.1,
  #             outputDirectory = getwd(),projectName = paste("goBP_no_redund",ctrast,sep = "_"))
  # 
  # # GO CCnoRedundant
  # WebGestaltR(enrichMethod = "ORA",organism = "hsapiens",minNum = 10,maxNum=2000,enrichDatabase = "geneontology_Cellular_Component_noRedundant",
  #             interestGene = gl,interestGeneType = "ensembl_gene_id",referenceSet = "genome",fdrThr = 0.1,
  #             outputDirectory = getwd(),projectName = paste("goCC_no_redund",ctrast,sep = "_"))
  # 
}


