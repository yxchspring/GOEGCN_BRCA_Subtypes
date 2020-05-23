### Obatin the weighted DEGs for classificaiton using ROC
# set your path
# setwd("your folder path")

# reference websites:
# https://cloud.tencent.com/developer/ask/42426
# https://stackoverflow.com/questions/17150183/plot-multiple-lines-in-one-graph
# https://www.cnblogs.com/li-20151130/p/9038524.html


rm(list = ls())
# load the required packages
library(caret)
library(klaR)
library(DMwR)
library(ROCR)

############## Step 1: Load data
### load the required data from "Main1_weighted_DEGs.R"

### load the data: brca_tumor_new, classes_tumor
load("./Data/ds2_brca_tumor_new.Rdata")

### load the weighted DEGs for classification: subtype_weighted_DEGs_classification
load("./Data/ds2_weighted_DEGs_classification.Rdata")
savepath <- "./classifiers_ROC/"
############## Step 2: constuct the binary classification models using manchine learning for each subtype
############## Step 2.1: The shared content is as follows,
split=0.60 # the coefficient for spliting traning and testing set
methods <- c("nb", "rf", "svmRadial") # four machine learning methods
repeats_num <- 1 # the repeat number for 5-fold cross-validation

############## Step 2.2: learn the binary classifier for each subtype
select_pred_subtype <- function(do.this){
  switch (do.this,
          Basal = {pred2_subtype_ROC <- prediction(pred_subtype_ROC$Basal, subtype_test_labels)},
          Her2 = {pred2_subtype_ROC <- prediction(pred_subtype_ROC$Her2, subtype_test_labels)},
          LumA = {pred2_subtype_ROC <- prediction(pred_subtype_ROC$LumA, subtype_test_labels)},
          LumB = {pred2_subtype_ROC <- prediction(pred_subtype_ROC$LumB, subtype_test_labels)},
          Normal = {pred2_subtype_ROC <- prediction(pred_subtype_ROC$Normal, subtype_test_labels)},
          stop("It should be one of {(Basal,Her2,LumA,LumB,Normal}")
  )
  return(pred2_subtype_ROC)
}


subtypes_names <- c("Basal","Her2","LumA","LumB","Normal")
subtypes_names_alias <- c("Basal-like","Her2","LumA","LumB","Normal-like")
subtype_count <- 0 #counter

for (ii in subtypes_names){
  cat("The work for ",ii," starts!\n")
  subtype_count <- subtype_count + 1
  
  # load the model_results_list and cfM_results_list
  load(paste0(savepath,"s3_classifiers_results_ROC_",ii,".Rdata"))
  
  # model_results_list <- list() # this variable is used to save the classification results for subtype ii
  # cfM_results_list <- list() # save the confusionMatrix results for subtype ii
  
  ############## Step 2.2.1: spliting the data into traning and testing set
  ### obtain the tumor data for subtype ii with weighted DEGs for classification
  brca_tumor_new_subtype <- brca_tumor_new[subtype_weighted_DEGs_classification[[subtype_count]],]
  brca_tumor_new_subtype_t <- as.data.frame(t(brca_tumor_new_subtype))
  # subtype_t_class <- ifelse(classes_tumor == "Basal","Basal","non_Basal")
  ### transfer the multiple classes labels into binary classes labels
  subtype_t_class <- ifelse(classes_tumor == ii,ii,paste0("non_",ii))
  brca_tumor_new_subtype_t$classes <- as.factor(subtype_t_class)
  
  set.seed(256)
  ### acquire the training and testing sets
  trainIndex_subtype <- createDataPartition(brca_tumor_new_subtype_t$classes, p=split, list=FALSE)
  subtype_train <- brca_tumor_new_subtype_t[ trainIndex_subtype,]
  subtype_test <- brca_tumor_new_subtype_t[-trainIndex_subtype,]
  
  ############## Step 2.2.2: conduct the binary classification model training
  # Basal   Her2   LumA   LumB Normal 
  # 192     82    564    207     40 
  
  ############## method_counter = 1
  model_subtype <- model_results_list[[1]] # load the model
  cfn_Mtx_subtype <-  cfM_results_list[[1]] # load the confusionMatrix
  
  ############## Step 2.2.3: plot the ROC curves and save these figures
  pred_subtype_ROC <- predict(model_subtype, subtype_test, type="prob")
  subtype_test_labels <- ifelse(subtype_test$classes==ii,1,0)
  
  # pred2_subtype_ROC <- prediction(pred_subtype_ROC$Basal, subtype_test_labels)
  pred2_subtype_ROC <- select_pred_subtype(ii)
  
  perf2_subtype_ROC <- performance(pred2_subtype_ROC,"tpr","fpr")
  auc2_subtype_ROC1 <- performance(pred2_subtype_ROC,measure = "auc",x.measure = "cutoff")
  
  ######## save the roc curve
  ### please create the ROC_graphs in advance
  file_subtype_ROC <- paste("./ROC_graphs/",ii,"_ROC.png",sep = "")
  # png(filename = file_subtype_ROC,width = 5,height = 5,res = 300)
  plot(perf2_subtype_ROC,colorize=FALSE, col="red", type = "l", lty=1, lwd=2, 
       cex=1.5, cex.main=1.2, cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)

  
  
  
  ############## method_counter = 2
  model_subtype <- model_results_list[[2]] # load the model
  cfn_Mtx_subtype <-  cfM_results_list[[2]] # load the confusionMatrix
  
  ############## Step 2.2.3: plot the ROC curves and save these figures
  pred_subtype_ROC <- predict(model_subtype, subtype_test, type="prob")
  subtype_test_labels <- ifelse(subtype_test$classes==ii,1,0)
  
  # pred2_subtype_ROC <- prediction(pred_subtype_ROC$Basal, subtype_test_labels)
  pred2_subtype_ROC <- select_pred_subtype(ii)
  
  perf2_subtype_ROC <- performance(pred2_subtype_ROC,"tpr","fpr")
  auc2_subtype_ROC2 <- performance(pred2_subtype_ROC,measure = "auc",x.measure = "cutoff")
  
  ######## save the roc curve
  
  plot(perf2_subtype_ROC,add = TRUE,colorize=FALSE, col="green", type = "l", lty=2, lwd=2, 
       cex=1.5, cex.main=1.2, cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)
  
  
  ############## method_counter = 3
  model_subtype <- model_results_list[[3]] # load the model
  cfn_Mtx_subtype <-  cfM_results_list[[3]] # load the confusionMatrix
  
  ############## Step 2.2.3: plot the ROC curves and save these figures
  pred_subtype_ROC <- predict(model_subtype, subtype_test, type="prob")
  subtype_test_labels <- ifelse(subtype_test$classes==ii,1,0)
  
  # pred2_subtype_ROC <- prediction(pred_subtype_ROC$Basal, subtype_test_labels)
  pred2_subtype_ROC <- select_pred_subtype(ii)
  
  perf2_subtype_ROC <- performance(pred2_subtype_ROC,"tpr","fpr")
  auc2_subtype_ROC3 <- performance(pred2_subtype_ROC,measure = "auc",x.measure = "cutoff")
  
  ######## save the roc curve
  plot(perf2_subtype_ROC,add = TRUE, colorize=FALSE, col="blue", type = "l", lty=3, lwd=2,
       cex=1.5, cex.main=1.2, cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)
  
  
  auc2_ROC_prefix <- c("AUC of nb: ", "AUC of rf: ", "AUC of svmRadial: ")
  auc2_ROC_values <- c(round(unlist(auc2_subtype_ROC1@y.values),digits = 4),
              round(unlist(auc2_subtype_ROC2@y.values),digits = 4),
              round(unlist(auc2_subtype_ROC3@y.values),digits = 4))
  

  legend("bottomright",legend=paste0(auc2_ROC_prefix, auc2_ROC_values),
         col = c("red","green","blue"),lty=c(1,2,3),bg="white",lwd=1, cex=1
         )
  
  dev.copy(png, file_subtype_ROC, width=6, height=6,res=300,units="in")
  dev.off()
  
  cat("The work for ",ii," ends!\n")

  
}

############## Step 3: save the results of the binary classification models for each subtype

cat("The work for binary classifiers learning for all BRCA subtypes is completed!\n")