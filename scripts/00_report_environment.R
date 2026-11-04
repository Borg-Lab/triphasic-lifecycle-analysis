#!/usr/bin/env Rscript

required_packages <- c(
  "AnnotationDbi", "BiocGenerics", "BiocManager", "ComplexHeatmap", "DESeq2",
  "GO.db", "GenomeInfoDb", "GenomicRanges", "IRanges", "S4Vectors",
  "SummarizedExperiment", "WGCNA", "ape", "circlize", "data.table",
  "dendextend", "dplyr", "edgeR", "ggbio", "ggplot2", "ggrastr", "ggridges",
  "gridExtra", "myTAI", "patchwork", "pheatmap", "purrr", "readr", "rlang",
  "rtracklayer", "scales", "stringr", "tibble", "tidyr", "tidyverse", "yaml"
)

output_path <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  commandArgs(trailingOnly = TRUE)[[1L]]
} else {
  "docs/session_info.txt"
}

installed <- installed.packages()
package_table <- data.frame(
  package = required_packages,
  installed = required_packages %in% rownames(installed),
  version = vapply(
    required_packages,
    function(package) {
      if (package %in% rownames(installed)) installed[package, "Version"] else NA_character_
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
report <- capture.output({
  cat("Triphasic lifecycle analysis environment\n")
  cat("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n\n")
  cat("Required package versions\n")
  write.table(package_table, sep = "\t", row.names = FALSE, quote = FALSE, na = "MISSING")
  cat("\nR session information\n")
  print(sessionInfo())
})
writeLines(report, output_path)
message("Wrote environment report: ", output_path)
