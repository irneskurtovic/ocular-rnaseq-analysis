
library(DESeq2)

tissues <- c("cornea", "limbus", "sclera")

for (tissue_name in tissues) {

  cat("\n============================\n")
  cat("Tissue:", tissue_name, "\n")
  cat("============================\n")

  # Load the saved DESeq2 model
  dds <- readRDS(
    paste0("results/dds_", tissue_name, ".rds")
  )

  # Recreate the differential expression results
  res <- results(
    dds,
    contrast = c("condition", "SARS-CoV-2", "mock"),
    alpha = 0.05
  )

  # Summarize statistical results
  cat("Total genes tested:", sum(!is.na(res$pvalue)), "\n")
  cat("Genes with NA adjusted p-values:",
      sum(is.na(res$padj)), "\n")
  cat("Genes with adjusted p-value < 0.05:",
      sum(res$padj < 0.05, na.rm = TRUE), "\n")

  cat("Genes with positive log2 fold change:",
      sum(res$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("Genes with negative log2 fold change:",
      sum(res$log2FoldChange < 0, na.rm = TRUE), "\n")

  # Save an MA plot
  png(
    paste0("results/figures/MA_", tissue_name, ".png"),
    width = 1200,
    height = 900,
    res = 150
  )

  plotMA(
    res,
    ylim = c(-5, 5),
    main = paste("MA plot:", tissue_name)
  )

  dev.off()
}

cat("\nDifferential expression QC completed.\n")
