
# 10_volcano_plots.R
# Create volcano plots for each tissue

library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

tissues <- c("cornea", "limbus", "sclera")

for (tissue in tissues) {

  file_path <- paste0(
    "results/tables/DE_", tissue,
    "_SARS-CoV-2_vs_mock_annotated.csv"
  )

  results <- read.csv(file_path, check.names = FALSE)

  # Prepare data for plotting
  results$plot_padj <- results$padj
  results$plot_padj[is.na(results$plot_padj)] <- 1
  results$plot_padj <- pmax(results$plot_padj, 1e-300)

  results$significance <- "Not significant"
  results$significance[
    results$padj < 0.05 & results$log2FoldChange > 0
  ] <- "Upregulated"
  results$significance[
    results$padj < 0.05 & results$log2FoldChange < 0
  ] <- "Downregulated"

  p <- ggplot(
    results,
    aes(x = log2FoldChange, y = -log10(plot_padj),
        color = significance)
  ) +
    geom_point(alpha = 0.65, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40") +
    scale_color_manual(values = c(
      "Downregulated" = "#2878B5",
      "Not significant" = "grey75",
      "Upregulated" = "#D9534F"
    )) +
    labs(
      title = paste("Volcano Plot:", tools::toTitleCase(tissue)),
      subtitle = "SARS-CoV-2 exposure vs mock control",
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10]("Adjusted p-value")),
      color = "Expression"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(
    filename = paste0("results/figures/11_volcano_", tissue, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

message("Volcano plots saved successfully.")
