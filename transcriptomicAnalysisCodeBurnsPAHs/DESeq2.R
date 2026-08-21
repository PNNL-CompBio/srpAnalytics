# Load DESeq2 library
library(DESeq2)

# Create function to generate normalized counts and log2FC data files
generateData <- function(){
  setwd("C:/Users/burnkyle/Box/Manuscript files/Chapter 2 Manuscript/RNA sequencing/GitHub files")
  counts <- read.csv("Raw_Counts/All_Raw_Counts_Final.csv", row.names = 1)
  meta <- read.csv("Raw_Counts/condition_info.csv",sep = ',',check.names = F)
  
  # Make sure that count column order matches metadata row order
  counts <- counts[,meta$x]
  rownames(meta) <- meta$x
  
  # Generate the DESeq data object, rounding counts to integers so DESeq can handle the data
  dds <- DESeqDataSetFromMatrix(countData=as.matrix(round(counts)), 
                                colData=meta, 
                                design=~condition)
  # Remove low read count genes
  dds <- dds[rowSums(counts(dds)) > 10, ]
  
  # Run the analysis
  dds <- DESeq(dds)
  
  # Save the analysis
  saveRDS(dds,file = 'dds_object_wea.RDS')
  
  # Generate normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  outputFilePath <- "All_Norm_Counts_Final.csv"
  write.table(normalized_counts, file=outputFilePath, quote=FALSE, sep=",", row.names=TRUE, col.names=TRUE)
  
  # Generate results frame
  comparisons <- c("Coculture_SW12", "Coculture_VehA", "Coculture_SFM", "Coculture_VehB", "HBEC_SW12", "HBEC_Veh")
  comparison_num <- length(comparisons)/2
  n <- 1
  for (i in 1:comparison_num) {
    diff_ex <- results(dds, contrast = c("condition", comparisons[n], comparisons[n+1])) # Control in last position
    summary(diff_ex)
    outputFilePath <- paste0(comparisons[n],"_log2FC.csv")
    write.table(diff_ex, file=outputFilePath, quote=FALSE, sep=",", row.names=TRUE, col.names=TRUE)
    
    #Reduce to DEGs
    diff_ex_sig <- diff_ex[ which(diff_ex$padj <= 0.05), ]
    outputFilePath <- paste0(comparisons[n],"_log2FC_DEGs.csv")
    write.table(diff_ex_sig, file=outputFilePath, quote=FALSE, sep=",", row.names=TRUE, col.names=TRUE)
    
    n <- n+2
  }
}

generateData()

# Create vector of differentially expressed genes (DEGs) 
PAH_groups <- c("Coculture_SW12", "Coculture_SFM", "HBEC_SW12")
DEGs <- c()
for (i in PAH_groups) {
  df <- read.csv(paste0(i,"_log2FC_DEGs.csv"))
  DEGs <- c(DEGs, rownames(df))
}
DEGs <- unique(DEGs)