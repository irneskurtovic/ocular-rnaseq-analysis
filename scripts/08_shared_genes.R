# 08_shared_genes.R
# Compare significant genes across ocular tissues

tissues <- c("cornea", "limbus", "sclera")

gene_lists <- list()

for (tissue in tissues) {

  file_path <- paste0(
    "results/tables/DE_",
    tissue,
    "_SARS-CoV-2_vs_mock_annotated.csv"
  )

  de <- read.csv(file_path, check.names = FALSE)

  # Keep statistically significant genes
  significant <- de[
    !is.na(de$padj) & de$padj < 0.05,
  ]

  # Store gene IDs and their direction of change
  gene_lists[[tissue]] <- significant[
    ,
    c("gene_id", "gene_symbol", "log2FoldChange", "padj")
  ]

  names(gene_lists[[tissue]]) <- c(
    "gene_id", "gene_symbol",
    paste0("log2FC_", tissue),
    paste0("padj_", tissue)
  )
}

# Find genes significant in all three tissues
shared_all <- Reduce(
  function(x, y) merge(x, y, by = "gene_id"),
  gene_lists
)

# Check whether the direction agrees across all three tissues
shared_all$consistent_direction <-
  (
    sign(shared_all$log2FC_cornea) ==
      sign(shared_all$log2FC_limbus)
  ) &
  (
    sign(shared_all$log2FC_limbus) ==
      sign(shared_all$log2FC_sclera)
  )

# Save genes significant in all three tissues
write.csv(
  shared_all,
  "results/tables/genes_shared_all_three_tissues.csv",
  row.names = FALSE
)

# Count genes shared by each pair of tissues
pairwise_summary <- data.frame()

pairs <- list(
  c("cornea", "limbus"),
  c("cornea", "sclera"),
  c("limbus", "sclera")
)

for (pair in pairs) {

  genes_a <- gene_lists[[pair[1]]]$gene_id
  genes_b <- gene_lists[[pair[2]]]$gene_id

  common_ids <- intersect(genes_a, genes_b)

  pairwise_summary <- rbind(
    pairwise_summary,
    data.frame(
      tissue_1 = pair[1],
      tissue_2 = pair[2],
      shared_significant_genes = length(common_ids)
    )
  )
}

# Summarize significant genes in each tissue
tissue_summary <- data.frame(
  tissue = tissues,
  significant_genes = sapply(gene_lists, nrow)
)

write.csv(
  pairwise_summary,
  "results/tables/pairwise_shared_gene_summary.csv",
  row.names = FALSE
)

write.csv(
  tissue_summary,
  "results/tables/tissue_significant_gene_summary.csv",
  row.names = FALSE
)

cat("\nSignificant genes per tissue:\n")
print(tissue_summary)

cat("\nGenes shared by pairs of tissues:\n")
print(pairwise_summary)

cat("\nGenes significant in all three tissues:",
    nrow(shared_all), "\n")

cat("\nDirection of change among genes shared by all three:\n")
print(table(shared_all$consistent_direction))

cat("\nTop shared genes:\n")

top_shared <- shared_all[
  order(shared_all$padj_sclera),
  c(
    "gene_id",
    "gene_symbol",
    "log2FC_cornea",
    "log2FC_limbus",
    "log2FC_sclera",
    "consistent_direction"
  )
]

print(head(top_shared, 15))

  # Classify shared genes by direction across all three tissues

  shared_all$direction_pattern <- ifelse(
    shared_all$log2FC_cornea > 0 &
      shared_all$log2FC_limbus > 0 &
      shared_all$log2FC_sclera > 0,
    "Upregulated_all_three",
    ifelse(
      shared_all$log2FC_cornea < 0 &
        shared_all$log2FC_limbus < 0 &
        shared_all$log2FC_sclera < 0,
      "Downregulated_all_three",
      "Mixed_direction"
    )
  )

  direction_summary <- as.data.frame(
    table(shared_all$direction_pattern)
  )

  colnames(direction_summary) <- c(
    "direction_pattern",
    "number_of_genes"
  )

  print(direction_summary)

  write.csv(
    shared_all[shared_all$direction_pattern == "Upregulated_all_three", ],
    "results/tables/shared_genes_upregulated_all_three.csv",
    row.names = FALSE
  )

  write.csv(
    shared_all[shared_all$direction_pattern == "Downregulated_all_three", ],
    "results/tables/shared_genes_downregulated_all_three.csv",
    row.names = FALSE
  )

  write.csv(
    shared_all[shared_all$direction_pattern == "Mixed_direction", ],
    "results/tables/shared_genes_mixed_direction.csv",
    row.names = FALSE
  )

  write.csv(
    direction_summary,
    "results/tables/shared_gene_direction_summary.csv",
    row.names = FALSE
  )

  cat(
    "\nProvjera ukupnog broja gena:",
    sum(direction_summary$number_of_genes),
    "\n"
  )

  # Check GO enrichment input files
  up_file <- "results/tables/shared_genes_upregulated_all_three.csv"
  down_file <- "results/tables/shared_genes_downregulated_all_three.csv"

  file.exists(up_file)
  file.exists(down_file)