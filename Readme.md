####################################################
Paper:
RNA-Seq-based breast cancer subtype classification using machine learning approaches

####################################################
This proceduce is provided by Xiangchun Yu
Date:
December 24, 2018
Email:
yxchspring0209@foxmail.com

###################### 1. Data ##############################
1) ds1_brca_clinic.Rdata:
The clinic data of BRCA


2) ds1_brca_tumor_counts.Rdata:
The BRCA tumor Counts data with the same genes and samples as the BRCA tumor FPKM data of 3)

3) ds1_brca_tumor_fpkm.Rdata:
The tumor data of BRCA after filtering out the genes whose mean values are less than 0.2 and variations are less than 2 across the tumor samples.


Note: The data 2) and 3) are too big on GitHub, so please refer to the following website for the details,
Link: https://pan.baidu.com/s/1pDIGDjasZ79Q82h8lPcoUQ
Extraction code: dk5d

4) ref_brca_normal_counts.Rdata:
The BRCA normal Counts data with same genes as 2), the data is for just for reference, and itis not used in this paper.

5) ref_brca_normal_fpkm.Rdata:
The BRCA normal FPKM data with same genes as 3), the data is for just for reference, and itis not used in this paper.

###################### 2. Main functions ####################
There are 5 main functions and the corresponding descriptions are as follows,
1) Main1_weighted_DEGs.R
This main function is used to obtain the final weighted DEGs for classification and pathway enrichment analysis.

2) Main2_classifier.R
This function is adopted to conduct the BRCA subtypes classification using weighted DEGs (for classification)
Learn the corresponding classifiers in training set and evaluate the performance in the testing set for each subtype

3) Main2_classifier_ROC.R
This function is similar with the 2), the difference is the "ROC" metric is adopted when training the model.
This function can plot the corresponding ROC curves and the ROC curves for each subtype are saved at directory "ROC_graphs".

4) main3_GeneCouples_pathways.R
This is the subfunction which is called by 5) Main3_PEGCN.R.
This function is used to conduct thepathway enrichment analysis using PEGCN with weighted DEGs (for pathway enrichment analysis)

5) Main3_PEGCN.R
This function is used to conduct the pathway enrichment analysis using PEGCN with weighted DEG, and the save the corresponding enriched pathways.
The results for each subtype is saved into Excel file and the list results are saved into Rdata file.