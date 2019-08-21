### Conduct pathway enrichment analysis for the BRCA subtypes using PEGCN with weighted DEGs
# set your path
# setwd("your folder path")
rm(list = ls())
# load the required packages
library(openxlsx)

############## Step 1: Load data
### load the brca_tumor_new data
load("./Data/ds2_brca_tumor_new.Rdata")

### load the weighted DEGs for pathway enrichment analysis
load("./Data/ds2_weighted_DEGs_pathway.Rdata")

### The fucntion for recalculate the p-values for PEGCN and collect and reorder the final enriched pathways
source('./main3_GeneCouples_pathways.R')
# Basal   Her2   LumA   LumB Normal 
# 192     82    564    207     40 

subtypes <- c("Basal","Her2", "LumA", "LumB", "Normal")

############## Step 2: Get the results using PEGCN with weighted DEGs
PEGCN_results_list <- list() # save the PEGCN results
savePath_PEGCN <- "./pathways_PEGCN/"
### Obtain the enriched pathways using PEGCN with weighted DEGs for each subtype
for (ii in c(1:5)) {
  PEGCN_results_list_ii <- GeneCouples_pathways(brca_tumor_new, subtype_weighted_DEGs_pathway[[ii]],
                                                classes_tumor,subtypes[ii], pvalueCutoff=0.05,corr_cutoff=0.3,
                                                ifsetdiff=TRUE)
  cat("The work for the the",subtypes[ii] ," ends! \n")
  ### save the results of subtype ii into xlsx file
  write.xlsx(PEGCN_results_list_ii, file=paste0(savePath_PEGCN,"ds4_",ii,"_PEGCN_results.xlsx"))
  ### save the results into the list varible
  PEGCN_results_list[[ii]] <- PEGCN_results_list_ii
}
save(PEGCN_results_list,file=paste0(savePath_PEGCN,"ds4_",ii,"_PEGCN_results.xlsx"))
cat("The work for the the ifsetdiff==true ends! \n")


















