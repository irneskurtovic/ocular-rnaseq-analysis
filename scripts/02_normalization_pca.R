
# 1. Load packages
library(DESeq2)

# 2. Load raw counts and sample metadata
counts <- read.delim(
  "data/raw/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/metadata/sample_metadata.csv",
  stringsAsFactors = FALSE
)

# 3. Match metadata to the count matrix
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(counts), , drop = FALSE]

stopifnot(identical(rownames(metadata), colnames(counts)))
stopifnot(all(c("tissue", "condition") %in% colnames(metadata)))

metadata$tissue <- factor(metadata$tissue)
metadata$condition <- factor(
  metadata$condition,
  levels = c("mock", "SARS-CoV-2")
)

# 4. Create a DESeq2 object for exploratory analysis
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData = metadata,
  design = ~ 1
)

# 5. Estimate size factors for normalization
dds <- estimateSizeFactors(dds)

cat("Estimated size factors:\n")
print(sizeFactors(dds))

# 6. Calculate normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# 7. Apply variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# 8. Save results
write.csv(
  data.frame(
    sample_id = colnames(dds),
    library_size = colSums(counts(dds)),
    size_factor = sizeFactors(dds)
  ),
  "results/tables/normalization_factors.csv",
  row.names = FALSE
)

saveRDS(dds, "results/dds_normalized.rds")
saveRDS(vsd, "results/vst_transformed.rds")

cat("\nNormalization and VST completed successfully.\n")

library(ggplot2)
# 9. Prepare PCA data
pca_data <- plotPCA(
  vsd,
  intgroup = c("tissue", "condition"),
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar")
)

# 10. Create PCA plot
pca_plot <- plotPCA(
  vsd,
  intgroup = c("tissue", "condition")
)

print(pca_plot)

# 11. Save PCA plot
ggsave(
  "results/figures/03_pca_plot.png",
  plot = pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# 12. Save PCA coordinates
write.csv(
  pca_data,
  "results/tables/pca_coordinates.csv",
  row.names = FALSE
)

cat("\nPCA plot and coordinates saved successfully.\n")
