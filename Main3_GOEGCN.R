### Conduct pathway enrichment analysis for the BRCA subtypes using GOEGCN with weighted DEGs
# set your path
setwd("your folder path")
rm(list = ls())
# load the required packages
library(openxlsx)

############## Step 1: Load data
### load the brca_tumor_new data
load("./ds2_brca_tumor_new.Rdata")

### load the weighted DEGs for pathway enrichment analysis
load("./ds2_weighted_DEGs_pathway.Rdata")

### The fucntion for recalculate the p-values for GOEGCN and collect and reorder the final enriched pathways
source('./main3_GeneCouples_pathways.R')
# Basal   Her2   LumA   LumB Normal 
# 192     82    564    207     40 

subtypes <- c("Basal","Her2", "LumA", "LumB", "Normal")

############## Step 2: Get the results using GOEGCN with weighted DEGs
GOEGCN_results_list <- list() # save the GOEGCN results
savePath_GOEGCN <- "./pathways_GOEGCN/"
### Obtain the enriched pathways using GOEGCN with weighted DEGs for each subtype
for (ii in c(1:5)) {
  GOEGCN_results_list_ii <- GeneCouples_pathways(brca_tumor_new, subtype_weighted_DEGs_pathway[[ii]],
                                                classes_tumor,subtypes[ii], pvalueCutoff=0.05,corr_cutoff=0.3,
                                                ifsetdiff=TRUE)
  cat("The work for the the",subtypes[ii] ," ends! \n")
  ### save the results of subtype ii into xlsx file
  write.xlsx(GOEGCN_results_list_ii, file=paste0(savePath_GOEGCN,"ds4_",ii,"_GOEGCN_results.xlsx"))
  ### save the results into the list varible
  GOEGCN_results_list[[ii]] <- GOEGCN_results_list_ii
}
save(GOEGCN_results_list,file=paste0(savePath_GOEGCN,"ds4_",ii,"_GOEGCN_results.xlsx"))
cat("The work for the the ifsetdiff==true ends! \n")
