### Obatin the weighted DEGs for classificaiton
# set your path
# setwd("your folder path")
rm(list = ls())
# load the required packages
library(caret)
library(klaR)
library(DMwR)

############## Step 1: Load data
### load the required data from "Main1_weighted_DEGs.R"

### load the data: brca_tumor_new, classes_tumor
load("./Data/ds2_brca_tumor_new.Rdata")

### load the weighted DEGs for classification: subtype_weighted_DEGs_classification
load("./Data/ds2_weighted_DEGs_classification.Rdata")
savepath <- "./classifiers/"
############## Step 2: constuct the binary classification models using manchine learning for each subtype
############## Step 2.1: The shared content is as follows,
split=0.60 # the coefficient for spliting traning and testing set
methods <- c("nb", "rf", "svmRadial") # four machine learning methods
repeats_num <- 100 # the repeat number for 5-fold cross-validation

############## Step 2.2: learn the binary classifier for each subtype
subtypes_names <- c("Basal","Her2","LumA","LumB","Normal")
subtype_count <- 0 #counter
for (ii in subtypes_names){
  cat("The work for ",ii," starts!\n")
  subtype_count <- subtype_count + 1
  model_results_list <- list() # this variable is used to save the classification results for subtype ii
  cfM_results_list <- list() # save the confusionMatrix results for subtype ii
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
  method_counter <- 0
  for (kk in methods) {
    cat("The work for ",kk,":",ii," starts!\n")
    method_counter <- method_counter + 1
    ### there is no need of sampling for subtype LumA
    if(method_counter==3){
      # define training control
      train_control_subtype <- trainControl(method = "repeatedcv",
                                          number = 5,
                                          repeats = repeats_num,
                                          verboseIter = FALSE,
                                          sampling = NULL
                                          # sampling = "smote"
                                          )
      
    }else{
      # define training control
      train_control_subtype <- trainControl(method = "repeatedcv",
                                          number = 5,
                                          repeats = repeats_num,
                                          verboseIter = FALSE,
                                          # sampling = NULL
                                          sampling = "smote")
      
    }
    # train the model
    set.seed(825)
    model_subtype <- caret::train(classes ~ .,
                                data = subtype_train,
                                method = kk,
                                # method = "svmLinear3",
                                preProcess = c("scale", "center"),
                                trControl = train_control_subtype)
    ### conduct the prediciton using the model_subtype
    pred_subtype <- predict(model_subtype, subtype_test[,-ncol(subtype_test)])
    ### get the confusion matrix
    cfn_Mtx_subtype <- confusionMatrix(pred_subtype, subtype_test$classes,positive = ii)
    
    model_results_list[[method_counter]] <- model_subtype # save the model
    cfM_results_list[[method_counter]] <- cfn_Mtx_subtype # save the confusionMatrix
    
    cat("The work for ",kk,":",ii," ends!\n")
  }
  save(model_results_list,cfM_results_list,file = paste0(savepath,"s3_classifiers_results_",ii,".Rdata"))
  
}

############## Step 3: save the results of the binary classification models for each subtype

cat("The work for binary classifiers learning for all BRCA subtypes is completed!\n")














### Obatin the weighted DEGs for classificaiton
# set your path
setwd("your folder path")
rm(list = ls())
# load the required packages
library(caret)
library(klaR)
library(DMwR)

############## Step 1: Load data
### load the required data from "Main1_weighted_DEGs.R"

### load the data: brca_tumor_new, classes_tumor
load("ds2_brca_tumor_new.Rdata")

### load the weighted DEGs for classification: subtype_weighted_DEGs_classification
load("ds2_weighted_DEGs_classification.Rdata")
savepath <- "./classifiers/"
############## Step 2: constuct the binary classification models using manchine learning for each subtype
############## Step 2.1: The shared content is as follows,
split=0.60 # the coefficient for spliting traning and testing set
methods <- c("nb", "rf", "svmRadial") # four machine learning methods
repeats_num <- 100 # the repeat number for 5-fold cross-validation

############## Step 2.2: learn the binary classifier for each subtype
subtypes_names <- c("Basal","Her2","LumA","LumB","Normal")
subtype_count <- 0 #counter
for (ii in subtypes_names){
  cat("The work for ",ii," starts!\n")
  subtype_count <- subtype_count + 1
  model_results_list <- list() # this variable is used to save the classification results for subtype ii
  cfM_results_list <- list() # save the confusionMatrix results for subtype ii
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
  method_counter <- 0
  for (kk in methods) {
    cat("The work for ",kk,":",ii," starts!\n")
    method_counter <- method_counter + 1
    ### there is no need of sampling for subtype LumA
    if(method_counter==3){
      # define training control
      train_control_subtype <- trainControl(method = "repeatedcv",
                                          number = 5,
                                          repeats = repeats_num,
                                          verboseIter = FALSE,
                                          sampling = NULL
                                          # sampling = "smote"
                                          )
      
    }else{
      # define training control
      train_control_subtype <- trainControl(method = "repeatedcv",
                                          number = 5,
                                          repeats = repeats_num,
                                          verboseIter = FALSE,
                                          # sampling = NULL
                                          sampling = "smote")
      
    }
    # train the model
    set.seed(825)
    model_subtype <- caret::train(classes ~ .,
                                data = subtype_train,
                                method = kk,
                                # method = "svmLinear3",
                                preProcess = c("scale", "center"),
                                trControl = train_control_subtype)
    ### conduct the prediciton using the model_subtype
    pred_subtype <- predict(model_subtype, subtype_test[,-ncol(subtype_test)])
    ### get the confusion matrix
    cfn_Mtx_subtype <- confusionMatrix(pred_subtype, subtype_test$classes,positive = ii)
    
    model_results_list[[method_counter]] <- model_subtype # save the model
    cfM_results_list[[method_counter]] <- cfn_Mtx_subtype # save the confusionMatrix
    
    cat("The work for ",kk,":",ii," ends!\n")
  }
  save(model_results_list,cfM_results_list,file = paste0(savepath,"s3_classifiers_results_",ii,".Rdata"))
  
}

############## Step 3: save the results of the binary classification models for each subtype

cat("The work for binary classifiers learning for all BRCA subtypes is completed!\n")
