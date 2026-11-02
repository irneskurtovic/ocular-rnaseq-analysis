
# 1. Load DESeq2
library(DESeq2)

# 2. Load raw counts and metadata
counts <- read.delim(
  "data/raw/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/metadata/sample_metadata.csv",
  stringsAsFactors = FALSE
)

# 3. Match metadata to count matrix
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(counts), , drop = FALSE]

stopifnot(identical(rownames(metadata), colnames(counts)))

metadata$tissue <- factor(
  metadata$tissue,
  levels = c("cornea", "limbus", "sclera")
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("mock", "SARS-CoV-2")
)

# 4. Confirm sample numbers
cat("Samples per tissue and condition:\n")
print(xtabs(~ tissue + condition, data = metadata))

cat("\nCount matrix dimensions:\n")
print(dim(counts))

cat("\nMetadata and count matrix match successfully.\n")


# 5. Differential expression analysis
# Separate SARS-CoV-2 vs mock comparison for each tissue

tissues <- c("cornea", "limbus", "sclera")

for (current_tissue in tissues) {

  cat("\n====================================\n")
  cat("Analyzing tissue:", current_tissue, "\n")
  cat("====================================\n")

  # Select samples from the current tissue
  sample_idx <- metadata$tissue == current_tissue

  counts_tissue <- counts[, sample_idx, drop = FALSE]
  metadata_tissue <- metadata[sample_idx, , drop = FALSE]

  # Match metadata rows to count-matrix columns
  rownames(metadata_tissue) <- metadata_tissue$sample_id
  metadata_tissue <- metadata_tissue[
    colnames(counts_tissue), ,
    drop = FALSE
  ]

  stopifnot(
    identical(rownames(metadata_tissue), colnames(counts_tissue))
  )

  # Set mock as the reference condition
  metadata_tissue$condition <- factor(
    metadata_tissue$condition,
    levels = c("mock", "SARS-CoV-2")
  )

  stopifnot(all(table(metadata_tissue$condition) == 3))

  # Filter genes with very low counts
  keep <- rowSums(counts_tissue) >= 10
  counts_filtered <- counts_tissue[keep, , drop = FALSE]

  cat("Genes before filtering:", nrow(counts_tissue), "\n")
  cat("Genes after filtering:", nrow(counts_filtered), "\n")

  # Create DESeq2 dataset
  dds_tissue <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = metadata_tissue,
    design = ~ condition
  )

  # Run DESeq2
  dds_tissue <- DESeq(dds_tissue)

  # Compare SARS-CoV-2 against mock
  res <- results(
    dds_tissue,
    contrast = c("condition", "SARS-CoV-2", "mock"),
    alpha = 0.05
  )

  # Convert results to a data frame
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  # Sort by adjusted p-value
  res_df <- res_df[
    order(res_df$padj, na.last = TRUE),
  ]

  # Save full results
  output_file <- paste0(
    "results/tables/DE_",
    current_tissue,
    "_SARS-CoV-2_vs_mock.csv"
  )

  write.csv(res_df, output_file, row.names = FALSE)

  # Save fitted DESeq2 object
  saveRDS(
    dds_tissue,
    paste0("results/dds_", current_tissue, ".rds")
  )

  # Summarize results
  n_tested <- sum(!is.na(res_df$pvalue))
  n_significant <- sum(res_df$padj < 0.05, na.rm = TRUE)

  cat("Genes with a valid p-value:", n_tested, "\n")
  cat("Genes with adjusted p-value < 0.05:", n_significant, "\n")
  cat("Saved:", output_file, "\n")
}

cat("\nAll three tissue-specific analyses completed.\n")
