# 07_top_genes.R
# Identify significantly differentially expressed genes

tissues <- c("cornea", "limbus", "sclera")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

top_genes_all <- data.frame()
significant_summary <- data.frame()

for (tissue in tissues) {

  input_file <- paste0(
    "results/tables/DE_",
    tissue,
    "_SARS-CoV-2_vs_mock_annotated.csv"
  )

  de <- read.csv(input_file, check.names = FALSE)

  # Keep genes with adjusted p-value below 0.05
  significant <- de[
    !is.na(de$padj) & de$padj < 0.05,
  ]

  # Sort by adjusted p-value
  significant <- significant[
    order(significant$padj),
  ]

  # Add tissue label
  significant$tissue <- tissue

  # Identify direction of expression change
  significant$direction <- ifelse(
    significant$log2FoldChange > 0,
    "Upregulated",
    ifelse(
      significant$log2FoldChange < 0,
      "Downregulated",
      "No change"
    )
  )

  # Save all significant genes for this tissue
  write.csv(
    significant,
    paste0(
      "results/tables/significant_genes_",
      tissue,
      ".csv"
    ),
    row.names = FALSE
  )

  # Collect top 20 genes for the combined table
  top20 <- head(significant, 20)
  top_genes_all <- rbind(top_genes_all, top20)

  # Summarize results
  significant_summary <- rbind(
    significant_summary,
    data.frame(
      tissue = tissue,
      significant_genes = nrow(significant),
      upregulated = sum(
        significant$log2FoldChange > 0,
        na.rm = TRUE
      ),
      downregulated = sum(
        significant$log2FoldChange < 0,
        na.rm = TRUE
      )
    )
  )
}

# Save combined results
write.csv(
  top_genes_all,
  "results/tables/top20_genes_per_tissue.csv",
  row.names = FALSE
)

write.csv(
  significant_summary,
  "results/tables/significant_gene_summary.csv",
  row.names = FALSE
)

print(significant_summary)

cat("\nTop genes per tissue:\n")

for (tissue in tissues) {
  cat("\n---", tissue, "---\n")

  print(
    head(
      top_genes_all[top_genes_all$tissue == tissue,
                    c("gene_id", "gene_symbol", "gene_name",
                      "log2FoldChange", "padj", "direction")],
      10
    )
  )
}