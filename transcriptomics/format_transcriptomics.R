library(dplyr)
library(data.table)

#' Format transcriptomics data 
#' 
#' @param transcripts_path Path to the transcriptomics data. The columns should be 
#'    "GeneID", "Gene", "condition", "indication", "Log2FoldChage", "adj_p_value"
#' @param chemical_ids_path Path to the key file that maps comparisons to their chemical ids. 
#'    Contains "condition", "CASNO", and "Concentration" 
#' @param out_path Path for resulting CSV to be written to
format_transcriptomics <- function(transcripts_path = "~/Git_Repos/srpAnalytics/bmd2Samps_v3/srpDEGstats.tsv",
                                   chemical_ids_path = "~/Git_Repos/srpAnalytics/transcriptomics/chemical_table.txt",
                                   out_path = "~/Downloads/transcriptomics_out.txt") {
  ###############
  ## READ DATA ##
  ###############
  
  transcripts <- fread(transcripts_path)
  chemical_table <- fread(chemical_ids_path)
  
  ##################
  ## FORMAT TABLE ##
  ##################
  
  # Merge chemical information, remove missing gene information
  transcript_data <- left_join(transcripts, chemical_table, by = "condition") %>%
    filter(!is.na(GeneID))
  
  ####################
  ## OUTPUT RESULTS ##
  ####################
  
  fwrite(transcript_data, out_path, sep = "\t", row.names = F, quote = F)
  
}

format_transcriptomics()
  