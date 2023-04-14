library(dplyr)
library(data.table)

#' Format transcriptomics data 
#' 
#' @param transcripts_path Path to the transcriptomics data. The columns should be 
#'    "GeneID", "Gene", "condition", "indication", "Log2FoldChage", "adj_p_value"
#' @param chemical_ids_path Path to the key file that maps comparisons to their chemical ids. 

tranascripts_path <- "~/Git_Repos/srpAnalytics/bmd2Samps_v3/srpDEGstats.tsv"
chemical_ids_path <- "~/Git_Repos/srpAnalytics/transcriptomics/chemical_table.txt"

transcripts <- fread(transcripts_path)