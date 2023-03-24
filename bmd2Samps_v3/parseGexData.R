##loadGeneExpressionDAta

library(readxl)
library(dplyr)
library(tidyr)

tab<-readxl::read_xlsx('data/ZF_gex.xlsx',skip=1)

ind_res<-tab%>% tidyr::pivot_longer(ends_with("Indication"),names_to='condition',values_to='indication')%>%
  dplyr::select(GeneID,Gene,condition,indication)
fc_res<-tab%>% 
  tidyr::pivot_longer(ends_with("Log2FoldChange"),names_to='condition',values_to='Log2FoldChange')%>%
  dplyr::select(GeneID,Gene,condition,Log2FoldChange)

pval_res<-tab%>%
  tidyr::pivot_longer(ends_with("Adj_p_value"),names_to='condition',values_to='adj_p_value')%>%
  dplyr::select(GeneID,Gene,condition,adj_p_value)

  
  