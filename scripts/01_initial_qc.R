
# 1. Load the raw count matrix and sample metadata

counts <- read.delim(
  "data/raw/GSE164073_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/metadata/sample_metadata.csv",
  stringsAsFactors = FALSE
)

# 2. Validate the relationship between counts and metadata

stopifnot(!anyDuplicated(colnames(counts)))
stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyDuplicated(metadata$sample_id))
stopifnot(setequal(colnames(counts), metadata$sample_id))

# Put metadata in exactly the same order as the count matrix
metadata <- metadata[
  match(colnames(counts), metadata$sample_id),
]

stopifnot(identical(metadata$sample_id, colnames(counts)))

# 3. Calculate basic sample-level QC metrics

sample_qc <- data.frame(
  sample_id = colnames(counts),
  tissue = metadata$tissue,
  condition = metadata$condition,
  library_size = colSums(counts),
  detected_genes = colSums(counts > 0)
)

# 4. Inspect the results

print(sample_qc)

cat("\nSummary of library sizes:\n")
print(summary(sample_qc$library_size))

cat("\nSummary of detected genes:\n")
print(summary(sample_qc$detected_genes))

# 5. Save the QC table

write.csv(
  sample_qc,
  "results/tables/initial_sample_qc.csv",
  row.names = FALSE
)

# 6. Plot library sizes

png(
  "results/figures/01_library_sizes.png",
  width = 1400,
  height = 900,
  res = 150
)

par(mar = c(9, 5, 4, 2))

barplot(
  sample_qc$library_size,
  names.arg = sample_qc$sample_id,
  las = 2,
  cex.names = 0.65,
  ylab = "Total gene counts",
  main = "RNA-seq library size by sample"
)

dev.off()

cat("\nQC table and library-size plot saved successfully.\n")
