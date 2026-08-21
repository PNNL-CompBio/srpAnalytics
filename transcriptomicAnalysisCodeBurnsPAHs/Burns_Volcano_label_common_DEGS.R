library(ggplot2)
library(ggrepel)

generateVolcanoPlots <- function(
    filePath        = "Resub_FC_qvalue.csv",
    comparisons     = c("SW12", "SFM", "HBEC_Alone"),
    shared_comparisons = c("SW12", "HBEC_Alone"),  # <-- which comps define the "shared" set
    shared_label    = NULL,      # legend text; NULL = auto-generated
    shared_colour   = "#009E1F",
    pval_thresh     = 0.05,
    n_labels_padj   = 10,        # top N by significance
    n_labels_fc     = 10,        # top N by |log2FC| (among significant)
    fc_label_thresh = 1,         # only FC-labelled if |log2FC| >= this
    max_shared_labels = 25,      # cap on how many shared genes get names
    shared_same_direction = FALSE,   # TRUE = must move the same way in both comps
    highlight_in_all_plots = TRUE,   # also colour shared genes on the SFM plot
    outdir = "."
){
  data1 <- read.csv(filePath, header = TRUE, stringsAsFactors = FALSE)
  
  # ---- display label, robust to missing 'gene'/'product' columns ----
  has_gene_col    <- "gene"    %in% colnames(data1)
  has_product_col <- "product" %in% colnames(data1)
  if (has_gene_col) {
    is_missing_gene <- is.na(data1$gene) | data1$gene == "NA" | data1$gene == ""
    if (has_product_col) {
      data1$display_label <- ifelse(is_missing_gene,
                                    paste0(data1$Gene, " (", data1$product, ")"),
                                    data1$gene)
    } else {
      data1$display_label <- ifelse(is_missing_gene, data1$Gene, data1$gene)
    }
  } else {
    data1$display_label <- data1$Gene
  }
  if (!has_product_col) data1$product <- NA_character_
  
  # ---- keep only comparisons whose columns exist ----
  ok <- vapply(comparisons, function(cp)
    all(c(paste0("log2FoldChange_", cp), paste0("padj_", cp)) %in% colnames(data1)),
    logical(1))
  if (any(!ok)) warning("Dropping comparisons with missing columns: ",
                        paste(comparisons[!ok], collapse = ", "))
  comparisons <- comparisons[ok]
  shared_comparisons <- intersect(shared_comparisons, comparisons)
  if (length(shared_comparisons) < 2)
    stop("Need at least 2 valid comparisons in 'shared_comparisons'.")
  
  # ------------------------------------------------------------------
  # significance / direction matrices
  # ------------------------------------------------------------------
  sig_mat <- vapply(comparisons, function(cp) {
    padj <- data1[[paste0("padj_", cp)]]
    !is.na(padj) & padj < pval_thresh
  }, logical(nrow(data1)))
  dir_mat <- vapply(comparisons, function(cp) {
    fc <- data1[[paste0("log2FoldChange_", cp)]]
    ifelse(is.na(fc), NA_real_, sign(fc))
  }, numeric(nrow(data1)))
  colnames(sig_mat) <- colnames(dir_mat) <- comparisons
  
  # ---- diagnostic: DEG counts + all pairwise overlaps ----
  cat("\n--- Significant genes (padj <", pval_thresh, ") per comparison ---\n")
  print(colSums(sig_mat))
  if (length(comparisons) >= 2) {
    cat("\n--- Pairwise overlaps ---\n")
    pairs <- combn(comparisons, 2, simplify = FALSE)
    for (pr in pairs) {
      n_ov <- sum(sig_mat[, pr[1]] & sig_mat[, pr[2]])
      n_sd <- sum(sig_mat[, pr[1]] & sig_mat[, pr[2]] &
                    !is.na(dir_mat[, pr[1]]) & !is.na(dir_mat[, pr[2]]) &
                    dir_mat[, pr[1]] == dir_mat[, pr[2]])
      cat(sprintf("%-12s & %-12s : %5d shared (%d same direction)\n",
                  pr[1], pr[2], n_ov, n_sd))
    }
    cat(sprintf("All %d comparisons        : %5d shared\n",
                length(comparisons), sum(rowSums(sig_mat) == length(comparisons))))
  }
  
  # ------------------------------------------------------------------
  # define the shared set from 'shared_comparisons'
  # ------------------------------------------------------------------
  shared_idx <- which(rowSums(sig_mat[, shared_comparisons, drop = FALSE]) ==
                        length(shared_comparisons))
  if (shared_same_direction && length(shared_idx) > 0) {
    keep <- apply(dir_mat[shared_idx, shared_comparisons, drop = FALSE], 1,
                  function(x) !any(is.na(x)) && length(unique(x)) == 1)
    shared_idx <- shared_idx[keep]
  }
  shared_genes <- data1$Gene[shared_idx]
  
  if (is.null(shared_label))
    shared_label <- paste0("Shared: ", paste(shared_comparisons, collapse = " & "))
  
  cat("\n==============================================\n")
  cat("Shared set =", shared_label, "\n")
  cat("Genes:", length(shared_genes),
      if (shared_same_direction) "(same direction enforced)" else "", "\n")
  cat("==============================================\n")
  
  # ---- shared gene table (all comparisons' stats, sorted by mean padj of the pair) ----
  if (length(shared_idx) > 0) {
    shared_tbl <- data.frame(
      Gene    = data1$Gene[shared_idx],
      Label   = data1$display_label[shared_idx],
      Product = data1$product[shared_idx],
      stringsAsFactors = FALSE
    )
    for (cp in comparisons) {
      shared_tbl[[paste0("log2FC_", cp)]] <- round(data1[[paste0("log2FoldChange_", cp)]][shared_idx], 3)
      shared_tbl[[paste0("padj_",  cp)]]  <- signif(data1[[paste0("padj_", cp)]][shared_idx], 3)
    }
    rank_padj <- rowMeans(as.matrix(shared_tbl[, paste0("padj_", shared_comparisons), drop = FALSE]))
    mean_absFC <- rowMeans(abs(as.matrix(shared_tbl[, paste0("log2FC_", shared_comparisons), drop = FALSE])))
    shared_tbl$mean_padj_shared  <- signif(rank_padj, 3)
    shared_tbl$mean_absFC_shared <- round(mean_absFC, 3)
    shared_tbl <- shared_tbl[order(rank_padj), ]
    fn <- file.path(outdir, sprintf("Volcano_SharedGenes_%s.csv",
                                    paste(shared_comparisons, collapse = "_")))
    write.csv(shared_tbl, fn, row.names = FALSE)
    cat("Wrote:", fn, "\n")
    print(utils::head(shared_tbl[, c("Label", "mean_padj_shared", "mean_absFC_shared")], 25),
          row.names = FALSE)
    
    # genes chosen for labelling: most significant across the shared pair
    shared_to_label <- shared_tbl$Gene[seq_len(min(max_shared_labels, nrow(shared_tbl)))]
  } else {
    shared_to_label <- character(0)
  }
  
  # ------------------------------------------------------------------
  # one volcano per comparison
  # ------------------------------------------------------------------
  for (comp in comparisons) {
    fc_col   <- paste0("log2FoldChange_", comp)
    padj_col <- paste0("padj_", comp)
    
    df <- data.frame(
      Gene          = data1$Gene,
      display_label = data1$display_label,
      product       = data1$product,
      log2FC        = data1[[fc_col]],
      padj          = data1[[padj_col]],
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$log2FC) & !is.na(df$padj), ]
    
    df$Significance <- "Not Significant"
    df$Significance[df$padj < pval_thresh & df$log2FC > 0] <- "Up"
    df$Significance[df$padj < pval_thresh & df$log2FC < 0] <- "Down"
    
    show_shared <- highlight_in_all_plots || comp %in% shared_comparisons
    df$Shared <- df$Gene %in% shared_genes
    if (show_shared)
      df$Significance[df$Shared & df$Significance != "Not Significant"] <- shared_label
    
    df$Significance <- factor(df$Significance,
                              levels = c(shared_label, "Up", "Down", "Not Significant"))
    df$negLog10Padj <- -log10(pmax(df$padj, .Machine$double.xmin))
    
    # ---- pick labels ----
    sig_rows <- which(df$Significance != "Not Significant")
    lab_rows <- integer(0)
    if (show_shared)
      lab_rows <- c(lab_rows, which(df$Gene %in% shared_to_label &
                                      df$Significance == shared_label))
    if (length(sig_rows) > 0) {
      lab_rows <- c(lab_rows,
                    sig_rows[order(df$padj[sig_rows])][seq_len(min(n_labels_padj, length(sig_rows)))])
      big_fc <- sig_rows[abs(df$log2FC[sig_rows]) >= fc_label_thresh]
      if (length(big_fc) > 0)
        lab_rows <- c(lab_rows,
                      big_fc[order(-abs(df$log2FC[big_fc]))][seq_len(min(n_labels_fc, length(big_fc)))])
    }
    lab_rows <- sort(unique(lab_rows))
    
    df$label <- ""
    df$label[lab_rows] <- df$display_label[lab_rows]
    
    if (length(lab_rows) > 0) {
      labeled_info <- data.frame(
        Gene         = df$Gene[lab_rows],
        Label        = df$display_label[lab_rows],
        Product      = df$product[lab_rows],
        log2FC       = round(df$log2FC[lab_rows], 2),
        padj         = signif(df$padj[lab_rows], 3),
        Significance = df$Significance[lab_rows],
        InSharedSet  = df$Shared[lab_rows]
      )
      labeled_info <- labeled_info[order(labeled_info$padj), ]
      cat("\n=== Labeled genes for:", comp, "===\n")
      print(labeled_info, row.names = FALSE)
      csv_name <- file.path(outdir, sprintf("Volcano_LabeledGenes_%s.csv", comp))
      write.csv(labeled_info, csv_name, row.names = FALSE)
      cat("Wrote:", csv_name, "\n")
    }
    
    # ---- plot ----
    df_bg     <- df[df$Significance != shared_label, ]
    df_shared <- df[df$Significance == shared_label, ]
    
    pal <- c("Up" = "blue", "Down" = "red", "Not Significant" = "grey75")
    pal[shared_label] <- shared_colour
    
    sub <- paste0("padj < ", pval_thresh, " (no FC cutoff)")
    if (show_shared)
      sub <- paste0(sub, "; ", shared_colour_txt <- shared_label,
                    " (n = ", nrow(df_shared), ")")
    
    p <- ggplot(df, aes(x = log2FC, y = negLog10Padj)) +
      geom_point(data = df_bg, aes(colour = Significance), alpha = 0.5, size = 2) +
      geom_point(data = df_shared, aes(colour = Significance),
                 size = 2.8, alpha = 0.95, shape = 16) +
      scale_colour_manual(name = NULL, values = pal, drop = FALSE) +
      geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", colour = "black") +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
      geom_text_repel(aes(label = label, colour = Significance),
                      size = 3.2, fontface = "italic",
                      max.overlaps = Inf, show.legend = FALSE,
                      box.padding   = 0.6,    # was 0.35
                      point.padding = 0.6,    # was 0.25 - key one
                      force         = 8,      # was 3
                      force_pull    = 0.1,    # let labels drift far from their point
                      max.iter      = 1e5,
                      max.time      = 2,
                      segment.size = 0.25, segment.color = "grey50",
                      min.segment.length = 0, seed = 42) +
      labs(title = paste0("Volcano Plot - ", comp),
           subtitle = sub,
           x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
      theme(text = element_text(size = 18),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank(),
            legend.key = element_rect(fill = "white", colour = NA),
            legend.position = "top")
    
    ggsave(file.path(outdir, sprintf("Volcano_%s_SharedHighlighted.png", comp)),
           plot = p, width = 10, height = 8, dpi = 300)
  }
  
  invisible(shared_genes)
}