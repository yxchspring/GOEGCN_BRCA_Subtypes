### Obatin the weighted DEGs for classificaiton and pathway enrichment analysis
# set your path
setwd("your folder path")
# setwd("E:/RstudioProjects/Subtypes/Subtypes_Rcodes")
rm(list = ls())
# load the required packages
### before R version 3.5
# source("http://bioconductor.org/biocLite.R")
# biocLite("CancerSubtypes")

###  R version 3.5 or greater
# # if exists "Rtools is required to build R packages ...", 
# # please change the mirrors in RStudio: Tools>global options>packages.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("CancerSubtypes")

library(CancerSubtypes)
############## Step 0: folders preparation
if (file.exists("./classifiers")==FALSE){
  dir.create("./classifiers")
}
if (file.exists("./classifiers_ROC")==FALSE){
  dir.create("./classifiers_ROC")
}
if (file.exists("./pathways_PEGCN")==FALSE){
  dir.create("./pathways_PEGCN")
}
if (file.exists("./ROC_graphs")==FALSE){
  dir.create("./ROC_graphs")
}

############## Step 1: Load data
### load the brca tumor FPKM data
load("./Data/ds1_brca_tumor_fpkm.Rdata")
### load the processed clinical data
load("./Data/ds1_brca_clinic.Rdata")

### load the brca tumor Counts data
# This data has the same samples and genes as the fpkm's
load("./Data/ds1_brca_tumor_counts.Rdata")

cat("The work for conducting the differential expression analysis starts!\n")

############## Step 2: Conduct the differential expression analysis
################## Step 2.1: obatin the new clinic, fpkm, counts, and classes label data
subtypes_names <- c("Basal","Her2","LumA","LumB","Normal")
### screen the brca_clinic$subtype_BRCA_Subtype_PAM50 and remove the items whose values are NA
sample_na <- which(is.na(brca_clinic$subtype_BRCA_Subtype_PAM50))
brca_clinic_new <- brca_clinic[-sample_na,]

### get the new the tumor data for fpkm
brca_tumor_new <- brca_tumor[,-sample_na]

### get the new the tumor data for counts
brca_tumor_counts_new <- brca_tumor_counts[,-sample_na]

### set the labels for the new tumor data
classes_tumor <- brca_clinic_new$subtype_BRCA_Subtype_PAM50

################## Step 2.2: conduct the differential expression analysis using edgeR
### save the DEGs after filtering the logFC and pvalue cutoff
subtype_DEGs_filter_list <- list()
subtype_count <- 0 #counter
for (ii in subtypes_names) {
  subtype_count <- subtype_count + 1
  ### obtain the index for experiment group
  index_subtype <- which(classes_tumor==ii) 
  
  ### get the results of DEGs using vomm+limma for RNA-Seq Counts data
  subtype_results=DiffExp.limma(Tumor_Data=brca_tumor_counts_new[,-index_subtype], # control group
                                Normal_Data=brca_tumor_counts_new[,index_subtype], # experiment group
                                group=NULL,topk=NULL,RNAseq=TRUE)
  ### retain the genes whose logFC >=0.5 and p-vlalue <-0.05
  index_logFC <- which(abs(subtype_results[[1]]$logFC)>=0.5)
  # index_Pvalue <- which(abs(subtype_results[[1]]$P.Value)<=0.05)
  index_Pvalue <- which(abs(subtype_results[[1]]$adj.P.Val)<=0.01)
  
  ### conduct the intersection operation
  index_logFC_Pvalue <- intersect(index_logFC,index_Pvalue) 
  ### get the initial DEGs
  DEGs_filter <- subtype_results[[1]]$ID[index_logFC_Pvalue] 
  ### save the initial to the list variable
  subtype_DEGs_filter_list[[subtype_count]] <- DEGs_filter 
}
names(subtype_DEGs_filter_list) <- subtypes_names

cat("The work for retrieving weights of genes from gene regulatroy network starts!\n")

############## Step 3: Acuqire the weights of genes from gene regulatroy network
############## Step 3.1: get the weights(ranking) for all the genes of tumor data
### load the gene regulatory network
data(Ranking)
### find the corresponding index of all genes of tumor data from Ranking
index_genes <- match(rownames(brca_tumor_new),Ranking$mRNA_TF_miRNA.v21_SYMBOL)
### remove the NA value for the rownames(brca_tumor_new)
index_genes_nonna <- index_genes[which(!(is.na(index_genes)))] 
### get the dataframe including the gene names and ranking(weight) data
DEGs_ranking <- data.frame(GeneName = rownames(brca_tumor_new)[which(!(is.na(index_genes)))],
                           Ranking[index_genes_nonna,],stringsAsFactors=FALSE)
### reorder the DEGs_ranking by descending order of ranking
DEGs_ranking_order <- DEGs_ranking[order(DEGs_ranking$ranking_default,decreasing = TRUE),]

### assgin the minimum to the iterm whose ranking_default is NA
index_rank=which(is.na(DEGs_ranking_order$ranking_default))
DEGs_ranking_order$ranking_default[index_rank]=min(DEGs_ranking_order$ranking_default,na.rm =TRUE)

############## Step 3.2:conduct the intersction between inital DEGs and top 2000 genes with high weightes, 
##############          and get the weighted DEGs for classification 
### the list for saving the weighted DEGs for classification
subtype_weighted_DEGs_classification <- list()
for (kk in c(1:length(subtypes_names))) {
  subtype_weighted_DEGs_classification[[kk]] <- intersect(subtype_DEGs_filter_list[[kk]],DEGs_ranking_order$GeneName[c(1:1000)])
}
names(subtype_weighted_DEGs_classification) <- subtypes_names

############## Step 3.3:conduct the intersction between inital DEGs and top 4000 genes with high weightes, 
##############          and get the weighted DEGs for classification 
### the list for saving the weighted DEGs for pathway enrichment analysis
subtype_weighted_DEGs_pathway <- list()
for (kk in c(1:length(subtypes_names))) {
  subtype_weighted_DEGs_pathway[[kk]] <- intersect(subtype_DEGs_filter_list[[kk]],DEGs_ranking_order$GeneName[c(1:3000)])
}
names(subtype_weighted_DEGs_pathway) <- subtypes_names

cat("The work for save results data starts!\n")
############## Step 4: save the final results
### save the data whose brca_clinic$subtype_BRCA_Subtype_PAM50 are not NA
save(brca_tumor_new,classes_tumor,file = "./Data/ds2_brca_tumor_new.Rdata")
save(brca_clinic_new,file = "./Data/ds2_brca_clinic_new.Rdata")

### save the data of weighted DEGs for classification
save(subtype_weighted_DEGs_classification,file = "./Data/ds2_weighted_DEGs_classification.Rdata")

### save the data of weighted DEGs for pathway enrichment analysis
save(subtype_weighted_DEGs_pathway,file = "./Data/ds2_weighted_DEGs_pathway.Rdata")
cat("The work is completed!\n")












