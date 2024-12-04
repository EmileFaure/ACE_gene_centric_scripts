
# Function inspired by the test for enrichment python function developped by Minot and Willis
# in this paper : https://link.springer.com/article/10.1186/s40168-019-0722-6

# From two tables of genes annotations, the function will compute contingency tables for each unique annotation
# Then, it will compute Fisher tests to identify annotations enriched in list 1 versus list 2
# Finally, p-values obtained are corrected using an adjustment method at the choice of the user

# User has to define the annotation column that is used (e.g. KEGG, EggNOG, etc...).

# The function returns a data frame with on eline for each unique annotation, and the corresponding p-values (raw and corrected)

# Author : Emile Faure

test_for_enrichment <- function(Annotation_Table_1,
  Annotation_Table_2,
  Annotation_Type,
  alpha=0.01,
  fdr_method="BH",
  minimum.occ=1000) {

  suppressPackageStartupMessages(require(tidyverse))
  
  # Count occurences of each unique annotation
  list1_vc <- Annotation_Table_1 %>% count({{Annotation_Type}})
  list2_vc <- Annotation_Table_2 %>% count({{Annotation_Type}})

  # Keep track of total number of genes in each list
  list1_len <- nrow(Annotation_Table_1)
  list2_len <- nrow(Annotation_Table_2)

  output = data.frame(matrix(NA,nrow=nrow(list1_vc),ncol=3))
  names(output)=c("Annotation","Odds_Ratio","p_value")

  for (i in c(1:nrow(list1_vc))){
    Annotation_i = as.character(select(list1_vc,{{Annotation_Type}}) %>% slice(i)) # Focal annotation
    Count_i_l1 = as.numeric(select(list1_vc,n) %>% slice(i)) # Number of occurences in list 1
    Count_i_l2 = as.numeric(list2_vc %>% filter({{Annotation_Type}}==Annotation_i) %>% select(n)) # Number of occurences in list 2
    if (is.na(Count_i_l2)) {
      Count_i_l2=0
    }

    Length_i_l1 = list1_len - Count_i_l1 # Number of lines that are NOT corresponding to focal annotation in list 1
    Length_i_l2 = list2_len - Count_i_l2 # Number of lines that are NOT corresponding to focal annotation in list 2

    if (Count_i_l1 + Count_i_l2 >= minimum.occ) {
      Contingency_Tab = matrix(c(Count_i_l1,Length_i_l1,Count_i_l2,Length_i_l2),nrow=2,ncol=2,byrow=F)
      fishtest = fisher.test(Contingency_Tab,alternative="greater")
      
      output[i,]=c(Annotation_i,fishtest$estimate,fishtest$p.value)
      } else {
      output[i,]=c(Annotation_i,NA,NA)
      }
    cat('Processing annotation',i,'from',nrow(list1_vc),'\n')
    }
  
  output$Odds_Ratio=as.numeric(output$Odds_Ratio)
  output$p_value=as.numeric(output$p_value)
  
  output$Corrected_p_value=p.adjust(output$p_value, method=fdr_method)
  output$Enriched=case_when(output$Corrected_p_value > alpha ~ "FALSE", output$Corrected_p_value <= alpha ~ "TRUE")
  output=output[order(output$Corrected_p_value),]  
  return(output)
  }


