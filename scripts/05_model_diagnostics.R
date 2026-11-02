
library(DESeq2)

tissues <- c("cornea", "limbus", "sclera")

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

diagnostic_summary <- data.frame()

for (tissue_name in tissues) {

  cat("\n============================\n")
  cat("Tissue:", tissue_name, "\n")
  cat("============================\n")

  dds <- readRDS(
    paste0("results/dds_", tissue_name, ".rds")
  )

  res <- results(
    dds,
    contrast = c("condition", "SARS-CoV-2", "mock"),
    alpha = 0.05
  )

  # 1. Pregled statistickih rezultata
  cat("Broj uzoraka:", ncol(dds), "\n")
  cat("Procijenjeni size factors:\n")
  print(sizeFactors(dds))

  cat("Broj NA p-values:", sum(is.na(res$pvalue)), "\n")
  cat("Broj NA adjusted p-values:", sum(is.na(res$padj)), "\n")
  cat("Broj znacajnih gena:", sum(res$padj < 0.05, na.rm = TRUE), "\n")

  # 2. Graf procjena disperzije
  png(
    paste0("results/figures/dispersion_", tissue_name, ".png"),
    width = 1200,
    height = 900,
    res = 150
  )

  plotDispEsts(dds, main = paste("Dispersion plot:", tissue_name))
  dev.off()

  # 3. Cookova distanca po uzorku
  cooks <- assays(dds)[["cooks"]]

  if (!is.null(cooks)) {

    png(
      paste0("results/figures/cooks_", tissue_name, ".png"),
      width = 1200,
      height = 900,
      res = 150
    )

    boxplot(
      log10(cooks + 1e-8),
      names = colnames(dds),
      las = 2,
      ylab = "log10(Cook's distance)",
      main = paste("Cook's distance:", tissue_name)
    )

    dev.off()
  }
  # 3b. Provjera Cook's distance praga
  cooks <- assays(dds)[["cooks"]]

  if (!is.null(cooks)) {
    cooks_cutoff <- mcols(dds)$maxCooks

    if (is.null(cooks_cutoff)) {
      cat("maxCooks prag nije dostupan u ovom objektu.\n")
    } else {
      max_cooks_per_gene <- apply(
        cooks,
        1,
        max,
        na.rm = TRUE
      )

      genes_above_cutoff <- sum(
        max_cooks_per_gene > cooks_cutoff,
        na.rm = TRUE
      )

      cat(
        "Broj gena s Cook's distance iznad maxCooks praga:",
        genes_above_cutoff,
        "\n"
      )
    }
  }

  # 4. Sažetak za poređenje tkiva
  diagnostic_summary <- rbind(
    diagnostic_summary,
    data.frame(
      tissue = tissue_name,
      samples = ncol(dds),
      genes_tested = sum(!is.na(res$pvalue)),
      na_pvalue = sum(is.na(res$pvalue)),
      na_padj = sum(is.na(res$padj)),
      significant_genes = sum(res$padj < 0.05, na.rm = TRUE)
    )
  )
}

write.csv(
  diagnostic_summary,
  "results/tables/model_diagnostic_summary.csv",
  row.names = FALSE
)

print(diagnostic_summary)

cat("\nModel diagnostics completed.\n")
