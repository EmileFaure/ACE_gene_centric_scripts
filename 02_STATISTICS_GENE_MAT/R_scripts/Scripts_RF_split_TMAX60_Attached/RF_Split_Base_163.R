library(tidyverse)
library(caret)
library(ranger)

Gene_mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached_163.rds")
metadata_GeneMat_KNN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")

metadata_GeneMat_KNN=metadata_GeneMat_KNN[which(metadata_GeneMat_KNN$ACE_seq_name %in% as.factor(row.names(Gene_mat))),]
rownames(metadata_GeneMat_KNN)=metadata_GeneMat_KNN$ACE_seq_name
metadata_GeneMat_KNN = select(metadata_GeneMat_KNN,-c("ACE_seq_name"))
metadata_GeneMat_KNN = metadata_GeneMat_KNN[match(rownames(Gene_mat),rownames(metadata_GeneMat_KNN)),]

new.row = as.data.frame(matrix(nrow = ncol(Gene_mat), ncol = 55))

for(i in c(1:ncol(Gene_mat))) {
  gene <- as.data.frame(Gene_mat[,i])
  colnames(gene) <- "Abundance"
  mtry.grid <-  expand.grid(mtry = 5:8, min.node.size = 5, splitrule = "variance") # Only tune mtry 

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

  # new.row <- c(names(Gene_mat)[i], mtry, rsquared, r2, mse, importance)
  new.row[i,] <- c(names(Gene_mat)[i], mtry, rsquared, r2, mse, importance)
  cat('Processing AGC', i, 'of', ncol(Gene_mat),'\n')
}

new.row = as.data.frame(new.row)
colnames(new.row) <- c("AGC_ID","mtry","Rsquared_rf","R2_caret", "MSE",paste("Importance",colnames(metadata_GeneMat_KNN),sep="."))
rownames(new.row) = new.row$AGC_ID
new.row = new.row[,-1]

new.row[] = lapply(new.row, function(x) {
  as.numeric(as.character(x))
})
new.row = as.data.frame(new.row)

write.table(new.row, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/random_forests_AGC_TMAX60GQ_NZVuniquecut20_Attached_163.txt", quote=FALSE,sep="\t")



