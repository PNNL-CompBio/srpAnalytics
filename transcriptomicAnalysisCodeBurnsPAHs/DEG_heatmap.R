# Load required packages
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

# Set working directory to folder containing normalized count data
# Normalized count data can be obtained from NCBI GEO
working_directory <- "C:/Users/burnkyle/Box/Manuscript files/Chapter 2 Manuscript/RNA sequencing/Data files/Normalized count files/"
setwd(working_directory)

# Read in GEO normalized counts and create a dataframe from all files
sample_nums <- 1:24
ST1 <- read.csv("ST1 normalized counts .csv")
i <- 0
for (x in c(2:length(sample_nums))) {
  next_sample <- read.csv(paste("ST",x, " normalized counts .csv", sep = ""))
  if (i == 0) {
    Norm_counts <- full_join(ST1, next_sample, by = "Gene.Symbol")
    i <- 1
  }
  else {
    Norm_counts <- full_join(Norm_counts, next_sample, by = "Gene.Symbol")
  }
}

# Convert the "Gene.Symbol" into row names
Norm_counts <- column_to_rownames(Norm_counts, var = "Gene.Symbol")

# Rename columns for simplicity
colnames(Norm_counts) <- c(
  "ST1", "ST2", "ST3", "ST4", 
  "ST5", "ST6", "ST7", "ST8", 
  "ST9", "ST10", "ST11", "ST12", 
  "ST13", "ST14", "ST15", "ST16",
  "ST17", "ST18", "ST19", "ST20",
  "ST21", "ST22", "ST23", "ST24"
  )

# Create vectors for each treatment group containing sample names 
Veh_Coculture_SW12_sample_names <- c("ST1", "ST2", "ST3", "ST4")
PAH_Coculture_SW12_sample_names <- c("ST5", "ST6", "ST7", "ST8")
Veh_Coculture_SFM_sample_names <- c("ST9", "ST10", "ST11", "ST12")
PAH_Coculture_SFM_sample_names <- c("ST13", "ST14", "ST15", "ST16")
Veh_HBEC_alone_SW12_sample_names <- c("ST17", "ST18", "ST19", "ST20")
PAH_HBEC_alone_SW12_sample_names <- c("ST21", "ST22", "ST23", "ST24")

# Calculate mean normalized count values for each gene for each vehicle group 
Norm_counts$CC_SW12_veh_mean <- rowMeans(Norm_counts[, c(Veh_Coculture_SW12_sample_names)])
Norm_counts$CC_SFM_veh_mean <- rowMeans(Norm_counts[, c(Veh_Coculture_SFM_sample_names)])
Norm_counts$HBEC_alone_SW12_veh_mean <- rowMeans(Norm_counts[, c(Veh_HBEC_alone_SW12_sample_names)])

# Create vectors of all sample names for each experiment
sample_names_CC_Sw12 <- c(Veh_Coculture_SW12_sample_names, PAH_Coculture_SW12_sample_names)
sample_names_CC_SFM <- c(Veh_Coculture_SFM_sample_names, PAH_Coculture_SFM_sample_names)
sample_names_HBEC_alone_Sw12 <- c(Veh_HBEC_alone_SW12_sample_names, PAH_HBEC_alone_SW12_sample_names)

# Calculate the log2 fold-change values for each replicate from normalized count data
# Create new columns with log2FC in data frame
for (name in sample_names_CC_Sw12) {
  Norm_counts[paste0(name, ".log2FC")] <- log(Norm_counts[name]/Norm_counts$CC_SW12_veh_mean, 2)
}
for (name in sample_names_CC_SFM) {
  Norm_counts[paste0(name, ".log2FC")] <- log(Norm_counts[name]/Norm_counts$CC_SFM_veh_mean, 2)
}
for (name in sample_names_HBEC_alone_Sw12) {
  Norm_counts[paste0(name, ".log2FC")] <- log(Norm_counts[name]/Norm_counts$HBEC_alone_SW12_veh_mean, 2)
}

# Create data frame with only log2FC data
log2FC <- select(Norm_counts, contains("log2FC")) 

# Filter data to include only differentially expressed genes (DEG) where adjusted p-value is <= 0.05
# This list of statistically significant DEGs comes from DESeq2 analysis
DEGs <- c("ZC3H12A", "MIR6732", "H2BC18", "H2BC21", "LOC124904412", "MTX1", 
          "BTG2", "AURKAP1", "LOC124905978", "RFX8", "LOC105373605", "NR4A2", 
          "LOC112268433", "LOC124906212", "CSRNP1", "NFKBIZ", "CXCL2", 
          "EGR1", "LOC401218", "DUSP1", "LOC124901135", "LOC124901256", 
          "H2BC4", "H3C4", "LOC101929555", "TNFAIP3", "ERV3-1", "PPP1R3B", 
          "TRIB1", "PTCSC2", "KLF4", "SLC25A25", "KLF6", "MAP3K8", "LOC124902455", 
          "PPP1R3C", "DUSP5", "FOSL1", "LOC105369669", "H4C16", "NR4A1", 
          "NR4A1AS", "DDIT3", "LOC128125814", "MIR616", "LINC01619", "PRKAB1", 
          "LOC124903033", "RGCC", "LINC00641", "G2E3-AS1", "SRSF5", "FAM161B", 
          "EHD4-AS1", "LOC124903484", "ANPEP", "THOC6", "LOC124903638", 
          "MT1H", "CENPBD1P", "MIR22HG", "LOC124903923", "RASD1", "PRR15L", 
          "H3-3B", "UBALD2", "LOC105371981", "MBP", "KLF2", "UPF1", "ZFP36", 
          "PLEKHG2", "TULP2", "NDUFAF5", "ZFAS1", "SNORD12B", "RUNX1-AS1", 
          "DIP2A-IT1", "MAFF", "TSPYL2", "MIR22HG_1", "TRNV", "ATP8", "LOC100996583", 
          "ALPL", "YRDC", "EDN2", "PDZK1IP1", "PCSK9", "JUN", "EPS8L3", 
          "FCGR2A", "ZNF648", "RGS2", "IKBKE", "IL19", "TRIB2", "RHOB", 
          "DUSP2", "NEURL3", "KYNU", "FLJ46875", "LOC105373714", "CCL20", 
          "IRAK2", "LOC105377642", "VIPR1", "SLC6A20", "TNNC1", "ABI3BP", 
          "TF", "CXCL8", "CXCL1", "CXCL3", "PARM1", "LOC107986289", "TNIP3", 
          "RNF150", "LOC124900826", "ENPP6", "MYL12BP2", "MTNR1A", "SLC9A3", 
          "SLC9A3-OT1", "ADRB2", "CD83", "H1-2", "IER3", "LINC03040", "LOC107986632", 
          "NCOA7", "SGK1", "SOD2", "IL6", "SLC26A4-AS1", "SLC26A4", "AHCYL2", 
          "LOC107986847", "LOC105379385", "TCIM", "PTP4A3", "GNE", "PRRX2", 
          "DEPP1", "IFIT2", "R3HCC1L", "LRRC55", "HYOU1-AS1", "LOC105369812", 
          "SOCS2-AS1", "SOCS2", "SLC5A8", "NFKBIA", "FOS", "TNFAIP2", "C2CD4A", 
          "C2CD4B", "CILP", "ZG16B", "IL17C", "DOC2B", "HS3ST3B1", "CSF3", 
          "MED24", "SOCS3", "C1QTNF1", "SERPINB7", "EBI3", "LDLR", "JUNB", 
          "IER2", "FCGBP", "CEACAM7", "RELB", "FUT2", "LOC105447645", "TP53INP2", 
          "SNAI1", "TNFRSF6B", "LOC105377139", "LSS", "ANKRD62P1-PARP4P3", 
          "LOC124905087", "PIM3", "MIR6821", "IL17REL", "HNRNPDLP5", "DEFB4A_1", 
          "ZG16B_1", "GCGR_1", "FCGBP_1", "NSDHL_1", "DEFB4B_2", "SLC5A8_1", 
          "DOC2B_1", "CEACAM7_1", "LSS_1", "IER3_1", "IER3_2", "IER3_3", 
          "IER3_4", "IER3_5", "CDC42EP5_3", "CYP4B1", "CAPN9", "LOC107985359", 
          "LOC124904546", "CYP1B1", "SYNPR-AS1", "LOC105377940", "PRKN", 
          "FZD9", "CYSRT1", "SYT8", "TNNT3", "INS-IGF2", "IGF2", "CYB561A3", 
          "CLCF1", "WTAPP1", "MMP1", "ADPRHL1", "CGNL1", "CYP1A1", "ALDH3A1", 
          "TMEM97", "TBC1D3D", "ADGRE2", "BICRA", "MFNG", "TNNT3_1")
log2FC_DEGs <- log2FC[DEGs,]

# Remove data on vehicle-treated samples and leave only data on PAH-treated samples
rep_vals_only <- log2FC_DEGs |>
  select(contains(c(PAH_Coculture_SW12_sample_names, 
           PAH_Coculture_SFM_sample_names, 
           PAH_HBEC_alone_SW12_sample_names)))

# Replace infinite values in data frame so that pheatmap() can function
rep_vals_only[is.infinite(as.matrix(rep_vals_only))] <- NA

### Create color scale
n <- 101  # odd so the middle color is exactly white
cols <- colorRampPalette(c("#0E4C92", "white", "#C21807"))(n)

# Min/max of 4 to better visualize data
minv <- -4
maxv <- 4

# Split breaks so that 0 lands exactly at the midpoint
n_low  <- floor(n/2) + 1   
n_high <- n - n_low + 1    

breaks <- c(seq(minv, 0, length.out = n_low),
        seq(0, maxv, length.out = n_high))  

# Remove one duplicated 0 to keep breaks length = n
breaks <- unique(breaks)

#Choose the legend ticks you want
legend_ticks <- c(minv, 0, maxv) #, 8, 12)

# Provide matching labels (you can customize the strings if desired)
legend_labs <- as.character(legend_ticks)

# Create custom colors for legend
light_blue <- rgb(148, 203, 236, maxColorValue = 255)
dark_blue <- rgb(046, 037, 133, maxColorValue = 255)

# Assign colors to groups
ann_colors <- list(
  Treatment = c("SW12" = "#606060", "SFM" = "#D4D4D4"),
  Model = c("HBEC alone" = light_blue, "Co-culture" = dark_blue))

ann_col <- data.frame(
  Treatment = factor(c(rep("SW12", 4), rep("SFM", 4), rep("SW12", 4))),
  Model = factor(c(rep("Co-culture", 8), rep("HBEC alone", 4))))

rownames(ann_col) <- colnames(rep_vals_only)

# Create heatmap
heatmap <- pheatmap(
  rep_vals_only,
  cluster_rows = T,
  cluster_cols = T,
  clustering_distance_cols = 'euclidean',
  clustering_distance_rows = 'euclidean',
  clustering_method = 'ward.D',
  border_color = "black",
  number_color = "black",
  col = cols,
  breaks = breaks,
  cellwidth = 15,
  cellheight = 1.5,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  annotation_names_col = F,
  show_colnames = F,
  show_rownames = F,
  fontsize_row = 4,
  fontsize_col = 5,
  legend = T,
  legend_labels = legend_labs,
  legend_breaks = legend_ticks#_in_range
)

# Save heatmap as .png
output_file <- "DEG_heatmap.png"

ggsave(
  output_file,
  plot = heatmap,
  width = 8,
  height = 8,
  units = "in",
  dpi = 1100
)