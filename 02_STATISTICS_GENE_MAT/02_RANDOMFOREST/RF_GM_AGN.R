require(tidyverse)
require(vegan)
require(caret)
require(ggplot2)
require(foreach)
require(doParallel)

# Import Data
Gene_mat = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_FreeLiving.rds")
metadata_GeneMat_KNN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")

# Select only CTD samples
metadata_GeneMat_KNN=metadata_GeneMat_KNN[which(metadata_GeneMat_KNN$ACE_seq_name %in% as.factor(row.names(Gene_mat))),]

# Match order of rows
rownames(metadata_GeneMat_KNN)=metadata_GeneMat_KNN$ACE_seq_name
metadata_GeneMat_KNN = select(metadata_GeneMat_KNN,-c("ACE_seq_name"))
metadata_GeneMat_KNN = metadata_GeneMat_KNN[match(rownames(Gene_mat),rownames(metadata_GeneMat_KNN)),]

cores=55
cl=makeCluster(cores)
registerDoParallel(cl)

RF.genes <- foreach(i=1:ncol(Gene_mat), .combine = rbind) %dopar% {
  
  require(caret)
  require(ranger)
  
  # Subset abundance data to one gene
  gene <- as.data.frame(Gene_mat[,i])
  colnames(gene) <- "Abundance" # Avoid also loading dplyr
  
  mtry.grid <-  expand.grid(mtry = 4:7, min.node.size = 5, splitrule = "variance") # Only tune mtry 
  
  repeatedCV <- trainControl(method = "repeatedcv", # 3 repeats of 4-fold cross validation
                             number = 4,
                             repeats = 3,
                             verboseIter = FALSE)
  
  set.seed(1994)
  mod <- caret::train(
    x = metadata_GeneMat_KNN,
    y = gene$Abundance,
    method = 'ranger',
    num.trees = 501,
    tuneGrid = mtry.grid,
    trControl = repeatedCV,
    importance = "permutation"
  )
  
  # Tuned parameter
  mtry <- mod$finalModel$mtry
  
  # Out of bag estimatin of R2
  rsquared <- mod$finalModel$r.squared # R squared of best model using OOB samples from RF
  mse <- mod$finalModel$prediction.error  # Mean Squared Error using OOB samples from RF
  r2 <- mean(mod$resample$Rsquared) # Same as fit$results$Rsquared[which(fit$results$RMSE==min(fit$results$RMSE))], the rsquared calculated by caret for the final model
  
  # Variable importance
  importance <- mod$finalModel$variable.importance # Mean Decrease in Accuracy (same order as in input table)
  
  new.row <- c(names(Gene_mat)[i], mtry, rsquared, r2, mse, importance)
  new.row
}

RF.genes <- as.data.frame(RF.genes)
colnames(RF.genes) <- c("AGC_ID","mtry","Rsquared_rf","R2_caret", "MSE",paste("Importance",colnames(metadata_GeneMat_KNN),sep="."))
saveRDS(RF.genes, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/random_forests_AGC_NZVuniquecut20_FreeLiving.rds")

