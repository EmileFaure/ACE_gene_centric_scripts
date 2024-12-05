Gene_mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached.rds")
k=1
j=1
while (k<ncol(Gene_mat)){
  Gene_mat_split=Gene_mat[,c(k:min(k+9999,ncol(Gene_mat)))]
  file_name=paste("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached",j,sep="_")
  file_name=paste(file_name,"rds",sep=".")
  saveRDS(Gene_mat_split,file_name)
  k=k+10000
  j=j+1
  }

