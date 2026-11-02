
# 09_go_enrichment.R
# Gene Ontology enrichment analysis

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Load genes upregulated in all three tissues
up_genes <- read.csv(
  "results/tables/shared_genes_upregulated_all_three.csv"
)

# Load genes downregulated in all three tissues
down_genes <- read.csv(
  "results/tables/shared_genes_downregulated_all_three.csv"
)

# Load differential expression results for all tissues
cornea <- read.csv(
  "results/tables/DE_cornea_SARS-CoV-2_vs_mock_annotated.csv"
)

limbus <- read.csv(
  "results/tables/DE_limbus_SARS-CoV-2_vs_mock_annotated.csv"
)

sclera <- read.csv(
  "results/tables/DE_sclera_SARS-CoV-2_vs_mock_annotated.csv"
)

# Inspect column names
names(up_genes)
names(down_genes)
names(cornea)
names(limbus)
names(sclera)

# Check dimensions
dim(up_genes)
dim(down_genes)
dim(cornea)
dim(limbus)
dim(sclera)


# Define the common gene universe across all three tissues

common_tested_genes <- Reduce(
  intersect,
  list(
    cornea$gene_id,
    limbus$gene_id,
    sclera$gene_id
  )
)

# Convert gene IDs to character vectors
common_tested_genes <- as.character(common_tested_genes)

# Prepare foreground gene lists
up_gene_ids <- unique(as.character(up_genes$gene_id))
down_gene_ids <- unique(as.character(down_genes$gene_id))

# Check the number of genes
cat("Genes tested in all three tissues:",
    length(common_tested_genes), "\n")

cat("Upregulated genes:",
    length(up_gene_ids), "\n")

cat("Downregulated genes:",
    length(down_gene_ids), "\n")

# Check that all foreground genes are in the universe
cat("Upregulated genes in universe:",
    sum(up_gene_ids %in% common_tested_genes), "\n")

cat("Downregulated genes in universe:",
    sum(down_gene_ids %in% common_tested_genes), "\n")


# Gene Ontology enrichment: Biological Process

# Upregulated genes
ego_up <- enrichGO(
  gene          = up_gene_ids,
  universe      = common_tested_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

# Downregulated genes
ego_down <- enrichGO(
  gene          = down_gene_ids,
  universe      = common_tested_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

# Inspect results
cat("GO terms for upregulated genes:",
    nrow(as.data.frame(ego_up)), "\n")

cat("GO terms for downregulated genes:",
    nrow(as.data.frame(ego_down)), "\n")

head(as.data.frame(ego_up))
head(as.data.frame(ego_down))

# Save GO enrichment results

write.csv(
  as.data.frame(ego_up),
  "results/tables/GO_BP_upregulated_all_three.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down),
  "results/tables/GO_BP_downregulated_all_three.csv",
  row.names = FALSE
)

# Load plotting package
library(enrichplot)

# Plot enriched biological processes
pdf("results/figures/09_GO_upregulated.pdf",
    width = 10, height = 7)

print(dotplot(ego_up, showCategory = 15) +
        ggplot2::ggtitle("GO Biological Process - Upregulated Genes"))

dev.off()

pdf("results/figures/10_GO_downregulated.pdf",
    width = 10, height = 7)

print(dotplot(ego_down, showCategory = 15) +
        ggplot2::ggtitle("GO Biological Process - Downregulated Genes"))

dev.off()
# Save GO analysis summary

go_summary <- data.frame(
  gene_group = c("Upregulated", "Downregulated"),
  input_genes = c(length(up_gene_ids), length(down_gene_ids)),
  genes_in_GO_test = c(
    sum(up_gene_ids %in% common_tested_genes),
    sum(down_gene_ids %in% common_tested_genes)
  ),
  background_genes = rep(length(common_tested_genes), 2),
  significant_GO_terms = c(
    nrow(as.data.frame(ego_up)),
    nrow(as.data.frame(ego_down))
  )
)

write.csv(
  go_summary,
  "results/tables/GO_analysis_summary.csv",
  row.names = FALSE
)

print(go_summary)
