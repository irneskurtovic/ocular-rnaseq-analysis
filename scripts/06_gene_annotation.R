
library(AnnotationDbi)
library(org.Hs.eg.db)

# Tkiva koja analiziramo
tissues <- c("cornea", "limbus", "sclera")

# Folder za rezultate
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Funkcija za anotaciju jedne DE tabele
annotate_de_table <- function(file_path, output_path) {

  # Ucitavanje originalnih DE rezultata
  de <- read.csv(file_path, check.names = FALSE)

  # Pretvaranje ID-jeva u tekst radi sigurnog mapiranja
  gene_ids <- unique(as.character(de$gene_id))

  # Mapiranje Entrez ID -> simbol i puni naziv gena
  annotation <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENTREZID",
    columns = c("SYMBOL", "GENENAME")
  )

  # Uklanjanje duplikata mapiranja za isti ID
  annotation <- annotation[
    !duplicated(annotation$ENTREZID),
  ]

  # Dodavanje anotacije uz ocuvanje originalnog redoslijeda
  de$gene_id <- as.character(de$gene_id)

  de$gene_symbol <- annotation$SYMBOL[
    match(de$gene_id, annotation$ENTREZID)
  ]

  de$gene_name <- annotation$GENENAME[
    match(de$gene_id, annotation$ENTREZID)
  ]

  # Broj uspjesno anotiranih ID-jeva
  mapped <- !is.na(de$gene_symbol) &
    de$gene_symbol != ""

  cat("\nFajl:", basename(file_path), "\n")
  cat("Ukupno redova:", nrow(de), "\n")
  cat("Anotirani geni:", sum(mapped), "\n")
  cat("Bez simbola:", sum(!mapped), "\n")

  # Spremanje u novi fajl
  write.csv(de, output_path, row.names = FALSE)

  return(data.frame(
    file = basename(file_path),
    total_rows = nrow(de),
    annotated = sum(mapped),
    unmapped = sum(!mapped)
  ))
}

# Obrada sva tri tkiva
annotation_summary <- data.frame()

for (tissue in tissues) {

  input_file <- paste0(
    "results/tables/DE_",
    tissue,
    "_SARS-CoV-2_vs_mock.csv"
  )

  output_file <- paste0(
    "results/tables/DE_",
    tissue,
    "_SARS-CoV-2_vs_mock_annotated.csv"
  )

  current_summary <- annotate_de_table(
    input_file,
    output_file
  )

  annotation_summary <- rbind(
    annotation_summary,
    current_summary
  )
}

# Cuvanje sazetka
write.csv(
  annotation_summary,
  "results/tables/gene_annotation_summary.csv",
  row.names = FALSE
)

print(annotation_summary)

cat("\nGene annotation completed.\n")

# Validate annotated DE tables against original results

tissues <- c("cornea", "limbus", "sclera")

stat_columns <- c(
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

validation_summary <- data.frame()

for (tissue in tissues) {

  original_file <- paste0(
    "results/tables/DE_", tissue, "_SARS-CoV-2_vs_mock.csv"
  )

  annotated_file <- paste0(
    "results/tables/DE_", tissue, "_SARS-CoV-2_vs_mock_annotated.csv"
  )

  original <- read.csv(original_file, check.names = FALSE)
  annotated <- read.csv(annotated_file, check.names = FALSE)

  same_gene_ids <- identical(
    as.character(original$gene_id),
    as.character(annotated$gene_id)
  )

  same_statistics <- all(vapply(stat_columns, function(column) {
    isTRUE(all.equal(
      original[[column]],
      annotated[[column]],
      tolerance = 1e-12
    ))
  }, logical(1)))

  validation_summary <- rbind(
    validation_summary,
    data.frame(
      tissue = tissue,
      same_number_of_rows = nrow(original) == nrow(annotated),
      same_gene_ids = same_gene_ids,
      same_statistics = same_statistics,
      annotated_genes = sum(
        !is.na(annotated$gene_symbol) &
          annotated$gene_symbol != ""
      )
    )
  )
}

print(validation_summary)

write.csv(
  validation_summary,
  "results/tables/annotation_validation.csv",
  row.names = FALSE
)
