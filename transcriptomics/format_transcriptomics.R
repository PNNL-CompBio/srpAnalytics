library(dplyr)
library(data.table)
library(zebrafish.db)

#' Format transcriptomics data 
#' 
#' @param transcripts_path Path to the transcriptomics data. The columns should be 
#'    "GeneID", "Gene", "condition", "indication", "Log2FoldChage", "adj_p_value"
#' @param out_path Path for resulting CSV to be written to
format_transcriptomics <- function(transcripts_path = "~/Git_Repos/srpAnalytics/bmd2Samps_v3/srpDEGstats.tsv",
                                   out_path = "~/Downloads/") {
  ###############
  ## READ DATA ##
  ###############
  
  transcripts <- fread(transcripts_path)
  
  ##################
  ## FORMAT TABLE ##
  ##################
  
  # Merge chemical information, remove missing gene information, add flag 
  transcript_data <- transcripts %>%
    filter(!is.na(GeneID)) %>%
    mutate(
      flag = factor(ifelse(adj_p_value <= 0.05, 1, 0) * ifelse(Log2FoldChange >= 0, 1, -1), levels = c(-1, 0, 1))
    )
  
  ##############
  ## ONTOLOGY ##
  ##############
  
  # Pull list of GO terms
  GO_Terms <- as.list(zebrafishGO)
  
  # Pull entrez IDs using zebra fish db
  Entrez <- as.list(zebrafishENTREZID)
  Entrez_DB <- data.table(
    ManufacturerID = names(Entrez) %>% unlist(),
    EntrezID = Entrez %>% unlist()
  ) %>%
    filter(!is.na(EntrezID)) %>%
    mutate(
      EntrezID = paste0("GeneID:", EntrezID),
      GOTerm = purrr::map(ManufacturerID, function(x) {
        names(GO_Terms[[x]]) %>% paste0(collapse = " ")
      }) %>% unlist()
    ) 
  
  # Pivot longer and unique
  Entrez_DB <- data.table(
    EntrezID = lapply(1:nrow(Entrez_DB), function(x) {
      theLength <- Entrez_DB$GOTerm[x] %>% strsplit(" ") %>% length()
      if (theLength == 0) {theLength == 1}
      rep(Entrez_DB$EntrezID[x], theLength)
    }) %>% unlist(),
    GOTerm = lapply(1:nrow(Entrez_DB), function(x) {
      theGos <- Entrez_DB$GOTerm[x]
      ifelse(theGos != "", strsplit(theGos, " ") %>% unlist(), NA)
    }) %>% unlist() 
  ) %>% 
    unique()
  
  browser()


  ####################
  ## OUTPUT RESULTS ##
  ####################
  
  fwrite(transcript_data %>% dplyr::select(GeneID, Gene, condition, Log2FoldChange, adj_p_value, flag), 
         file.path(out_path, "transcriptomics_out.txt"), sep = "\t", row.names = F, quote = F)
  
}

format_transcriptomics()
  