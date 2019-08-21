GeneCouples_pathways <- function(Data, WDEGs,classes_tumor,subtype_ii, pvalueCutoff=0.05,corr_cutoff=0.3,ifsetdiff=TRUE) {
  
  # BiocManager::install(c("GO.db","org.Hs.eg.db","clusterProfiler"))
  
  library(CancerSubtypes)
  library(GO.db)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  # library(psych)

  
  ###################### Step 1: load the data
  # load("s6_pathway_DEGs.Rdata")
  # Basal   Her2   LumA   LumB Normal 
  # 192     82    564    207     40 
  
  ###################### Step 2: conduct the pathway analysis for subtype_ii
  ###################### Step 2.1: obtain the WDGS for the RNA-Seq based data
  Data_WDEGs_Tumor <- Data[WDEGs,]
  ###################### Step 2.2: reassign the labels for the RNA-Seq based data
  # data with the WDGS
  # Data_WDEGs_Tumor <- Data[WDEGs,]
  Data_WDEGs_Tumor_norm <- data.normalization(Data_WDEGs_Tumor,type="feature_Median",log2=TRUE)
  
  # data for the subtype_ii
  Data_WDEGs_Tumor_experi <- Data_WDEGs_Tumor_norm[,which(classes_tumor==subtype_ii)]
  # data for the non subtype_ii
  Data_WDEGs_Tumor_ctrl <- Data_WDEGs_Tumor_norm[,which(classes_tumor!=subtype_ii)]
  
  ###################### Step 2.3: conduct the GO analysis for the subtype_ii with WDGS
  # conduct the enrichGO analysis
  enrichGO_pathways <- enrichGO(gene = WDEGs, 
                                ont = "BP", 
                                keyType  = 'SYMBOL',
                                OrgDb ="org.Hs.eg.db",
                                pvalueCutoff = pvalueCutoff)
  enrichGO_pathways@result <- enrichGO_pathways@result[which(enrichGO_pathways@result$p.adjust<= pvalueCutoff),]
  
  # simplify the results of the enrichGO analysis
  enrichGO_pathways_sim_09  <- simplify(enrichGO_pathways, 
                                        cutoff=0.9, 
                                        by="p.adjust", 
                                        select_fun=min)
  
  enrichGO_pathways_sim_09.result <- enrichGO_pathways_sim_09@result[which(enrichGO_pathways_sim_09@result$p.adjust <= pvalueCutoff),]
  
  ###################### Step 2.4: construct the PCC newtwork for the control group
  # The subtype_ii is regarded as the experimental group and the non subtype_ii is regared as the control groups
  # corr_cutoff <- 0.3
  
  ###################### Step 2.4.1: construct the similarity matrix
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  ### obatin PCC (Pearson correalation coefficient ) for non subtype_ii data (Control group)
  # transpose the dataframe into sample * features
  matrix_ctrl <- t(Data_WDEGs_Tumor_ctrl)
  
  # calculate the PCC network
  PCC_ctrl <- cor(matrix_ctrl, method="pearson")
  # set the diag elements to be 0
  diag(PCC_ctrl) <- 0 
  
  # the correlation is set to be positive
  PCC_ctrl_corr <- abs(PCC_ctrl) 
  # the elements who are NA is set to be 0
  PCC_ctrl_corr[is.na(PCC_ctrl_corr)] <- 0 
  # get the upper triangle matrix
  PCC_ctrl_corr_upper <- get_upper_tri(PCC_ctrl_corr)
  
  # get the significant positions (coordinate)
  # Sig_location_ctrl_upper <- which(PCC_ctrl_corr_upper >= corr_cutoff)
  
  # get the matrix with significant positions (coordiante)
  length_WDGS <- length(WDEGs)
  matrix_Sig_location_ctrl_upper <-  matrix(NA,ncol = length_WDGS,nrow = length_WDGS)
  rownames(matrix_Sig_location_ctrl_upper) <- WDEGs
  colnames(matrix_Sig_location_ctrl_upper) <- WDEGs
  
  # matrix_Sig_location_ctrl_upper[Sig_location_ctrl_upper] <- 100
  
  ###################### Step 2.5: construct the PCC newtwork for the experimental group
  # The subtype_ii is regarded as the experimental group and the non subtype_ii is regared as the control groups
  ###################### Step 2.5.1: construct the similarity matrix
  ### obatin PCC (Pearson correalation coefficient ) for subtype_ii data (experimental group)
  matrix_experi <- t(Data_WDEGs_Tumor_experi)
  
  PCC_experi <- cor(matrix_experi, method="pearson")
  diag(PCC_experi) <- 0
  
  PCC_experi_corr <- abs(PCC_experi) # the correlation is set to be positive
  PCC_experi_corr[is.na(PCC_experi_corr)] <- 0 # the elements who are NA is set to be 0
  
  # get the upper triangle matrix
  PCC_experi_corr_upper <- get_upper_tri(PCC_experi_corr)
  
  # get the significant positions (coordinate)
  # Sig_location_experi_upper <- which(PCC_experi_corr_upper >= corr_cutoff)
  
  # get the matrix with significant positions (coordiante)
  # length_WDGS <- length(WDEGs)
  matrix_Sig_location_experi_upper <-  matrix(NA,ncol = length_WDGS,nrow = length_WDGS)
  rownames(matrix_Sig_location_experi_upper) <- WDEGs
  colnames(matrix_Sig_location_experi_upper) <- WDEGs
  
  # matrix_Sig_location_experi_upper[Sig_location_experi_upper] <- 100
  
  ###################### Step 3: compute the novel p-value for pathway analysis of non subtype_ii
  if(ifsetdiff==TRUE){
    # get the significant positions (coordinate) for control group
    Sig_location_ctrl_upper <- which(PCC_ctrl_corr_upper >= corr_cutoff)
    # get the significant positions (coordinate) for the experimental group
    Sig_location_experi_upper <- which(PCC_experi_corr_upper >= corr_cutoff)
    # get the shared locations between Sig_location_ctrl_upper and Sig_location_experi_upper
    shared_location <- intersect(Sig_location_ctrl_upper,Sig_location_experi_upper)
    
    ### Then, we get the new location variables
    Sig_location_ctrl_upper <- setdiff(Sig_location_ctrl_upper,shared_location)
    Sig_location_experi_upper <- setdiff(Sig_location_experi_upper,shared_location)
    ### Then, we get the new location matrix
    matrix_Sig_location_ctrl_upper[Sig_location_ctrl_upper] <- 100
    matrix_Sig_location_experi_upper[Sig_location_experi_upper] <- 100
    
  }else{
    # get the significant positions (coordinate)
    Sig_location_ctrl_upper <- which(PCC_ctrl_corr_upper >= corr_cutoff)
    matrix_Sig_location_ctrl_upper[Sig_location_ctrl_upper] <- 100
    
    # get the significant positions (coordinate) for experimental group
    Sig_location_experi_upper <- which(PCC_experi_corr_upper >= corr_cutoff)
    matrix_Sig_location_experi_upper[Sig_location_experi_upper] <- 100
  }
  
  
  ###################### Step 3.1: get the all parameters ready for non subtype_ii: N,K,n,k
  
  # obtain the N for ctrl
  N_ctrl <- (length_WDGS-1)*length_WDGS/2
  # obtain the K for ctrl
  K_ctrl <- length(Sig_location_ctrl_upper)
  NK_ctrl <- paste0(K_ctrl,"/",N_ctrl) # background
  
  # obtain the n and k for each item of pathway result
  nk_ctrl <- list()
  GOcount_ctrl <- list()
  GO_pvalue_ctrl <- list()
  counter <- 0
  length_GO_reuslt <- dim(enrichGO_pathways_sim_09.result)[1]
  for (ii in c(1:length_GO_reuslt)) {
    counter <- counter + 1
    # get the GO terms one by one
    geneID_GO <- unlist(strsplit(enrichGO_pathways_sim_09.result$geneID[ii] ,split="/"))
    # retrieve the each row and col from the matrix_Sig_location_ctrl_upper matrix
    matrix_Sig_upper_GO <- matrix_Sig_location_ctrl_upper[geneID_GO,geneID_GO]
    ### set the n, k for each term of GO term set
    len_geneID_GO <- length(geneID_GO)
    n_ctrl <- (len_geneID_GO - 1) * len_geneID_GO/2
    k_ctrl <- length(which(matrix_Sig_upper_GO==100))
    # M <- K
    # m+n <- N
    # x <- k
    # k <- n
    # calcualte the p-value for the novel method
    ctrl_pvalue <- 1 - phyper(k_ctrl-1, K_ctrl, N_ctrl-K_ctrl, n_ctrl)
    # gene couples ratio
    nk_ctrl[[counter]] <- paste(k_ctrl,"/",n_ctrl,sep = "")
    # gene couples count
    GOcount_ctrl[[counter]] <- k_ctrl
    # p-value for each term
    GO_pvalue_ctrl[[counter]] <- ctrl_pvalue
    
  }
  GO_pvalue_ctrl <- unlist(GO_pvalue_ctrl) 
  GO_padjust_ctrl <- NULL
  GO_padjust_ctrl <- p.adjust(GO_pvalue_ctrl, method = "BH", n = length(GO_pvalue_ctrl)) ############ p.adjust
  
  ###################### Step 3.2: contruct a dataframe to save the novel pathway results
  # until now, we have obtained the following parameters for computeing p-value
  # 1) N_ctrl, K_ctrl, NK_ctrl, n_ctrl, k_ctrl, nk_ctrl
  # 2) GOcount_ctrl, GO_pvalue_ctrl
  # Now we should create the new dataframe
  
  GeneRatio_ctrl <- unlist(nk_ctrl)
  BgRation_ctrl <- rep(NK_ctrl,length_GO_reuslt)
  # GO_pvalue_ctrl
  # GO_padjust_ctrl
  GOcount_ctrl <- unlist(GOcount_ctrl)
  
  Genecouples_pathways_ctrl <- data.frame("ID"=enrichGO_pathways_sim_09.result$ID,
                                               "Description"=enrichGO_pathways_sim_09.result$Description,
                                               "GeneRatio"=GeneRatio_ctrl,
                                               "BgRatio"=BgRation_ctrl,
                                               "pvalue"=GO_pvalue_ctrl,
                                               "p.adjust"=GO_padjust_ctrl,
                                               "count"=GOcount_ctrl
  )
  
  rownames(Genecouples_pathways_ctrl) <- Genecouples_pathways_ctrl$ID
  
  ###################### Step 3.3: screen the results whose pvaue is less than pvalueCutoff
  # Genecouples_pathways_ctrl <- Genecouples_pathways_ctrl[Genecouples_pathways_ctrl$pvalue <= pvalueCutoff,]
  Genecouples_pathways_ctrl <- Genecouples_pathways_ctrl[Genecouples_pathways_ctrl$p.adjust <= pvalueCutoff,]
  Genecouples_pathways_ctrl_order <- Genecouples_pathways_ctrl[order(Genecouples_pathways_ctrl$p.adjust,
                                                                                   decreasing = FALSE),]

  cat("The pathway analysis reust for ctrl is saved!\n")
  
  ###################### Step 4: compute the novel p-value for pathway analysis of subtype_ii
  ###################### Step 4.1: get the all parameters ready for subtype_ii: N,K,n,k
  # obtain the N for subtype_ii
  N_experi <- (length_WDGS-1)*length_WDGS/2
  # obtain the K for subtype_ii
  K_experi <- length(Sig_location_experi_upper)
  NK_experi <- paste0(K_experi,"/",N_experi) # background
  
  # obtain the n and k for each item of pathway result
  nk_experi <- list()
  GOcount_experi <- list()
  GO_pvalue_experi <- list()
  counter <- 0
  # length_GO_reuslt <- dim(enrichGO_pathways_sim_09.result)[1]
  for (ii in c(1:length_GO_reuslt)) {
    counter <- counter + 1
    # get the GO terms one by one
    geneID_GO <- unlist(strsplit(enrichGO_pathways_sim_09.result$geneID[ii] ,split="/"))
    # retrieve the each row and col from the matrix_Sig_location_subtype_ii_upper matrix
    matrix_Sig_upper_GO <- matrix_Sig_location_experi_upper[geneID_GO,geneID_GO]
    ### set the n, k for each term of GO term set
    len_geneID_GO <- length(geneID_GO)
    n_experi <- (len_geneID_GO - 1) * len_geneID_GO/2
    k_experi <- length(which(matrix_Sig_upper_GO==100))
    # M <- K
    # m+n <- N
    # x <- k
    # k <- n
    # calcualte the p-value for the novel method
    experi_pvalue <- 1 - phyper(k_experi-1, K_experi, N_experi-K_experi, n_experi)
    # gene couples ratio
    nk_experi[[counter]] <- paste(k_experi,"/",n_experi,sep = "")
    # gene couples count
    GOcount_experi[[counter]] <- k_experi
    # p-value for each term
    GO_pvalue_experi[[counter]] <- experi_pvalue
    
  }
  GO_pvalue_experi <- unlist(GO_pvalue_experi) 
  GO_padjust_experi <- NULL
  GO_padjust_experi <- p.adjust(GO_pvalue_experi, method = "BH", n = length(GO_pvalue_experi)) ############ p.adjust
  
  ###################### Step 4.2: contruct a dataframe to save the novel pathway results
  # until now, we have obtained the following parameters for computeing p-value
  # 1) N_experi, K_experi, NK_experi, n_experi, k_experi, nk_experi
  # 2) GOcount_experi, GO_pvalue_experi
  # Now we should create the new dataframe
  
  GeneRatio_experi <- unlist(nk_experi)
  BgRation_experi <- rep(NK_experi,length_GO_reuslt)
  # GO_pvalue_experi
  # GO_padjust_experi
  GOcount_experi <- unlist(GOcount_experi)
  Genecouples_pathways_experi <- data.frame("ID"=enrichGO_pathways_sim_09.result$ID,
                                                 "Description"=enrichGO_pathways_sim_09.result$Description,
                                                 "GeneRatio"=GeneRatio_experi,
                                                 "BgRatio"=BgRation_experi,
                                                 "pvalue"=GO_pvalue_experi,
                                                 "p.adjust"=GO_padjust_experi,
                                                 "count"=GOcount_experi
  )
  
  rownames(Genecouples_pathways_experi) <- Genecouples_pathways_experi$ID
  
  ###################### Step 4.3: screen the results whose pvaue is less than pvalueCutoff
  Genecouples_pathways_experi <- Genecouples_pathways_experi[Genecouples_pathways_experi$p.adjust <= pvalueCutoff,]
  Genecouples_pathways_experi_order <- Genecouples_pathways_experi[order(Genecouples_pathways_experi$p.adjust,
                                                                                       decreasing = FALSE),]

  cat("The pathway analysis results for experi is saved!\n")
  
  ###################### Step 5: save the Control, Experiment, Common groups of enriched pathways using PEGCN with weighted DEGs
  ### get the shared GO terms ID
  experi_ID <- Genecouples_pathways_experi_order$ID
  ctrl_ID <- Genecouples_pathways_ctrl_order$ID
  shared_ID <- intersect(experi_ID,ctrl_ID)
  
  ### save the final result
  # for the control group
  Genecouples_pathway_ctrl_unique <- Genecouples_pathways_ctrl_order[setdiff(ctrl_ID,shared_ID),]
  
  # for experimental group
  Genecouples_pathway_experi_unique <- Genecouples_pathways_experi_order[setdiff(experi_ID,shared_ID),]
  
  # for the shared ones:
  Genecouples_pathway_shared <- Genecouples_pathways_experi_order[shared_ID,]
  cat(paste0("The pathway analysis results using PEGCN with weighted DEGS for ",subtype_ii, " is completed!\n"))
  return(list(Control=Genecouples_pathway_ctrl_unique,
              Experiment=Genecouples_pathway_experi_unique,
              Common=Genecouples_pathway_shared))
  
}

















