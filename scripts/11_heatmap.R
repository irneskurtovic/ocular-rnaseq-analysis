
# 11_heatmap.R
library(DESeq2)
library(pheatmap)

# 1. Load VST-transformed expression data
vst_data <- readRDS("results/vst_transformed.rds")
vst_mat <- assay(vst_data)

# 2. Load genes significant in all three tissues
shared_genes <- read.csv(
  "results/tables/genes_shared_all_three_tissues.csv"
)

# 3. Rank shared genes by their strongest adjusted significance
shared_genes$best_padj <- pmin(
  shared_genes$padj_cornea,
  shared_genes$padj_limbus,
  shared_genes$padj_sclera,
  na.rm = TRUE
)

shared_genes <- shared_genes[
  is.finite(shared_genes$best_padj),
]
shared_genes <- shared_genes[
  order(shared_genes$best_padj),
]

# 4. Select up to 50 genes present in the VST matrix
# Convert gene IDs to character before matching
shared_genes$gene_id <- as.character(shared_genes$gene_id)

selected <- head(
  shared_genes$gene_id[
    shared_genes$gene_id %in% rownames(vst_mat)
  ],
  50
)

if (length(selected) < 2) {
  stop("Too few shared genes matched the VST matrix.")
}

heatmap_mat <- vst_mat[selected, , drop = FALSE]

# 5. Standardize expression for each gene
heatmap_mat <- t(scale(t(heatmap_mat)))
heatmap_mat <- heatmap_mat[
  apply(heatmap_mat, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

# 6. Create sample annotations
sample_info <- as.data.frame(colData(vst_data))

annotation_col <- data.frame(
  Tissue = factor(sample_info$tissue),
  Condition = factor(sample_info$condition)
)
rownames(annotation_col) <- colnames(heatmap_mat)

# 7. Use gene symbols as row labels when available
gene_labels <- shared_genes$gene_symbol[
  match(rownames(heatmap_mat), shared_genes$gene_id)
]

valid_labels <- !is.na(gene_labels) & gene_labels != ""
gene_labels[!valid_labels] <- rownames(heatmap_mat)[!valid_labels]

# 8. Save heatmap
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pheatmap(
  heatmap_mat,
  labels_row = make.unique(gene_labels),
  annotation_col = annotation_col,
  clustering_method = "complete",
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  main = "Shared Differentially Expressed Genes Across Ocular Tissues",
  filename = "results/figures/12_shared_genes_heatmap.png",
  width = 12,
  height = 10
)

message("Heatmap saved successfully.")
