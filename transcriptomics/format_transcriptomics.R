# Basic formatting libraries 
library(dplyr)
library(data.table)

# Zebrafish database information
library(zebrafish.db)

# Gene Ontology Enrichment Package
library(clusterProfiler)

#' Format transcriptomics data 
#' 
#' @param transcripts_path Path to the transcriptomics data. The columns should be 
#'    "GeneID", "Gene", "condition", "indication", "Log2FoldChage", "adj_p_value"
#' @param chemical_table_path A metadata file with "condition", "Chemical_ID", "Concentration"
#'    written in uM, and "Species" 
#' @param out_path Path for resulting CSVs to be written to
format_transcriptomics <- function(transcripts_path = "~/Git_Repos/srpAnalytics/bmd2Samps_v3/srpDEGstats.tsv",
                                   chemical_table_path = "~/Git_Repos/srpAnalytics/transcriptomics/chemical_table.txt",
                                   out_path = "~/Downloads/") {
  ###############
  ## READ DATA ##
  ###############
  
  # Read data 
  transcripts <- fread(transcripts_path)
  chemical_table <- fread(chemical_table_path)
  
  ################################
  ## MAKE GENE CONVERSION TABLE ##
  ################################
  
  # Expand the provided Entrez ID's to include manufacturer and zfin
  
  # Pull entrez and zfin
  Entrez <- as.list(zebrafishENTREZID)
  Zfin <- as.list(zebrafishZFIN)
  
  # Build the entrez table
  Entrez_Table <- data.table(
    ManufacturerID = names(Entrez) %>% unlist(),
    EntrezID = Entrez %>% unlist()
  ) %>%
    filter(!is.na(EntrezID)) %>%
    mutate(
      GeneID = paste0("GeneID:", EntrezID)
    ) %>% 
    dplyr::select(-EntrezID)
  
  # Build the zfin table
  Zfin_Table <- data.table(
    ManufacturerID = names(Zfin) %>% unlist(),
    ZfinID = Zfin %>% unlist()
  )
  
  # Make the Gene Conversion Table 
  Gene_Conversion_Table <- left_join(Entrez_Table, Zfin_Table, by = "ManufacturerID") 
  
  ##################################
  ## FORMAT TRANSCRIPTOMICS TABLE ##
  ##################################
  
  # Extract unique Entrez Gene - Zfin page combinations
  EntrezZfin <- Gene_Conversion_Table %>%
    dplyr::select(-ManufacturerID) %>%
    unique() %>%
    filter(!is.na(ZfinID)) %>%
    mutate(Gene_URL = paste0("https://zfin.org/", ZfinID)) %>%
    dplyr::select(-ZfinID)
  
  # 1. Merge chemical information
  # 2. Remove missing gene information
  # 3. Add flag for directionality of expression
  # 4. Add Chemical_ID, Concentration, and Species
  # 5. Create new "Comparison" variable for tracking all unique experiments (i.e. this chemical at this concentration with this species)
  # 6. Merge zfin urls 
  # 7. Clean up variable names and order
  transcript_data <- transcripts %>%
    filter(!is.na(GeneID)) %>%
    mutate(
      Flag = factor(ifelse(adj_p_value <= 0.05, 1, 0) * ifelse(Log2FoldChange >= 0, 1, -1), levels = c(-1, 0, 1))
    ) %>%
    dplyr::left_join(chemical_table, "condition") %>%
    dplyr::mutate(Comparison = paste(Chemical_ID, Concentration, Species)) %>%
    dplyr::left_join(EntrezZfin, EntrezZfin, by = "GeneID") %>%
    dplyr::rename(AdjPValue = adj_p_value) %>%
    dplyr::select(Chemical_ID, Comparison, Species, Concentration, GeneID, Gene, 
                  Gene_URL, Log2FoldChange, AdjPValue, Flag)
  
  ##############
  ## ONTOLOGY ##
  ##############
  
  # Pull list of GO terms
  GO_Terms <- as.list(zebrafishGO)
  
  # Make GO_Term Table 
  GO_Table <- data.table(
    ManufacturerID = names(GO_Terms) %>% unlist()
  ) %>%
    mutate(GO_ID = lapply(ManufacturerID, function(x) {
      names(GO_Terms[[x]]) %>% unlist() %>% paste0(collapse = " ")
    }) %>% unlist()) %>%
    filter(GO_ID != "")
  
  # Pivot longer and unique
  GO_Table <- data.table(
    ManufacturerID = lapply(1:nrow(GO_Table), function(x) {
      theLength <- GO_Table$GO_ID[x] %>% strsplit(" ") %>% length()
      if (theLength == 0) {theLength == 1}
      rep(GO_Table$ManufacturerID[x], theLength)
    }) %>% unlist(),
    GOTerm = lapply(1:nrow(GO_Table), function(x) {
      theGos <- GO_Table$GO_ID[x]
      ifelse(theGos != "", strsplit(theGos, " ") %>% unlist(), NA)
    }) %>% unlist() 
  ) %>% 
    unique() %>%
    left_join(Gene_Conversion_Table, by = "ManufacturerID")

  # 1. Remove any p-values greater than 0.05
  
  browser()
  
  transcript_data %>% 
    filter(AdjPValue <= 0.05) %>%
    dplyr::select(Comparison, GeneID) %>%
    dplyr::mutate(GeneID = gsub("GeneID:", "", GeneID)) 
  
  #res <- enrichGO(de, 'zebrafish.db', ont="ALL", pvalueCutoff=0.01, pAdjustMethod = "bonferroni")
  # res@result %>% dplyr::select(ID, Description, Count, ONTOLOGY) %>% arrange(-Count) %>% head(10)

  ####################
  ## OUTPUT RESULTS ##
  ####################
  
  # Write transcript table
  fwrite(transcript_data %>% dplyr::select(GeneID, Gene, condition, Log2FoldChange, adj_p_value, flag), 
         file.path(out_path, "transcriptomics_out.txt"), sep = "\t", row.names = F, quote = F)
  
  # Write ontology table
  
  # Write pathway table 
  
}

format_transcriptomics()
  