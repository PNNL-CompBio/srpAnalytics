# Basic formatting libraries 
library(dplyr)
library(data.table)

# Zebrafish database information
library(zebrafish.db)

# Gene Ontology & Pathway Enrichment Package
library(clusterProfiler)

#' Format transcriptomics data 
#' 
#' @param transcripts_path Path to the transcriptomics data. The columns should be 
#'    "GeneID", "Gene", "condition", "indication", "Log2FoldChage", "adj_p_value"
#' @param chemical_table_path A metadata file with "condition", "Chemical_ID", "Concentration"
#'    written in uM, and "Species" 
#' @param out_path Path for resulting CSVs to be written to
#' @param messages TRUE/FALSE to indicate whetere 
format_transcriptomics <- function(transcripts_path = "~/Git_Repos/srpAnalytics/bmd2Samps_v3/srpDEGstats.tsv",
                                   chemical_table_path = "~/Git_Repos/srpAnalytics/transcriptomics/chemical_table.txt",
                                   out_path = "~/Downloads/",
                                   messages = TRUE) {
  ###############
  ## READ DATA ##
  ###############
  
  if (messages) {message("...Reading data")}
  
  # Read data 
  transcripts <- fread(transcripts_path)
  chemical_table <- fread(chemical_table_path)
  
  ################################
  ## MAKE GENE CONVERSION TABLE ##
  ################################
  
  if (messages) {message("...Making transcript table")}
  
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
  
  ##############################################
  ## ISOLATE SIGNIFICANT GENES PER COMPARISON ##
  ##############################################

  # 1. Remove any p-values greater than 0.05
  # 2. Returns gene per comparison without GeneID tag
  GenesOfInterest <- transcript_data %>% 
    filter(AdjPValue <= 0.05) %>%
    dplyr::select(Comparison, GeneID) %>%
    dplyr::mutate(GeneID = gsub("GeneID:", "", GeneID)) 
    
  ##############
  ## ONTOLOGY ##
  ##############

  # Run an enrichment analysis per comparison 
  if (messages) {message("...Running gene ontology enrichment")}
  
  GO_Enrichment <- do.call(rbind, lapply(unique(GenesOfInterest$Comparison), function(comp) {
    
    if (messages) {message(paste0("......on comparison: ", comp))}
    
    # Select the genes 
    DE <- GenesOfInterest %>% filter(Comparison == comp) %>% select(GeneID) %>% unlist()
    
    # Run the enrichment analysis
    Enriched <- enrichGO(DE, 'zebrafish.db', ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "bonferroni")
    
    if (is.null(Enriched) || nrow(Enriched@result) == 0 || nrow(Enriched@result[Enriched@result$p.adjust <= 0.05,]) == 0) {
      return(
        data.frame(Comparison = comp, Ontology = NA, OntologyCount = NA, OntologyType = NA, OntologyPValue = NA)
      )
    }
    
    # Build the appropriate table of the top 10 hits  
    EnrichTable <- Enriched@result %>% 
      dplyr::select(Description, Count, ONTOLOGY, p.adjust) %>% 
      rename(Ontology = Description, OntologyCount = Count, OntologyType = ONTOLOGY, OntologyPValue = p.adjust) %>%
      arrange(-OntologyPValue) %>% 
      filter(OntologyPValue <= 0.05) %>%
      head(10) %>%
      mutate(OntologyType = ifelse(OntologyType == "BP", "Biological Process",
                            ifelse(OntologyType == "MF", "Molecular Function", "Cellular Component")),
             Comparison = comp
      ) %>%
      dplyr::select(Comparison, Ontology, OntologyCount, OntologyType, OntologyPValue)
    
    # Remove gene ontology IDs
    row.names(EnrichTable) <- 1:nrow(EnrichTable)
    
    return(EnrichTable)
    
  }))
  
  ##############
  ## PATHWAYS ##
  ##############
  
  # Run an enrichment analysis per comparison 
  if (messages) {message("...Running pathway enrichment")}
  
  Pathway_Enrichment <- do.call(rbind, lapply(unique(GenesOfInterest$Comparison), function(comp) {
    
    if (messages) {message(paste0("......on comparison: ", comp))}
    
    # Select the genes 
    DE <- GenesOfInterest %>% filter(Comparison == comp) %>% select(GeneID) %>% unlist()
    
    # Run the enrichment analysis
    Pathways <- enrichKEGG(DE, organism = "dre", pvalueCutoff = 0.05, pAdjustMethod = "bonferroni")
    
    if (is.null(Pathways) || nrow(Pathways@result) == 0 || nrow(Pathways@result[Pathways@result$p.adjust <= 0.05,]) == 0) {
      return(
        data.frame(Comparison = comp, Pathway = NA, PathwayCount = NA, PathwayPValue = NA)
      )
    }
    
    # Build the appropriate table of the top 10 hits  
    PathwayTable <- Pathways@result %>% 
      dplyr::select(Description, Count, p.adjust) %>% 
      rename(Pathway = Description, PathwayCount = Count, PathwayPValue = p.adjust) %>%
      arrange(PathwayPValue) %>% 
      filter(PathwayPValue <= 0.05) %>%
      head(10) %>%
      mutate(Comparison = comp) %>%
      dplyr::select(Comparison, Pathway, PathwayCount, PathwayPValue)
    
    # Remove gene ontology IDs
    row.names(PathwayTable) <- 1:nrow(PathwayTable)
    
    return(PathwayTable)
    
  }))
  
  ####################
  ## OUTPUT RESULTS ##
  ####################
  
  # Write transcript table
  fwrite(transcript_data, file.path(out_path, "gene.txt"), sep = "\t", row.names = F, quote = F, na = "NA")
  
  # Write ontology table
  fwrite(GO_Enrichment, file.path(out_path, "ontology.txt"), sep = "\t", row.names = F, quote = F, na = "NA")
  
  # Write pathway table 
  fwrite(Pathway_Enrichment, file.path(out_path, "pathway.txt"), sep = "\t", row.names = F, quote = F, na = "NA")
  
}

format_transcriptomics()
  