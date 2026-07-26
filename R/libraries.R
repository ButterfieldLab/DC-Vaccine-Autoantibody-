# ============================================================================
# libraries.R
# ----------------------------------------------------------------------------
# Package loader for the PhIP-seq / REAP autoantibody manuscript pipeline.
#
# This file was reconstructed from the original (lost) source file
# `autoantibody_libraries.R`. It loads every package referenced anywhere in
# the manuscript pipeline. Packages are loaded defensively: a missing package
# produces a warning listing what to install rather than aborting the session,
# so scripts that only need a subset of packages still run.
#
# Install helper (run once):
#   source("R/libraries.R"); install_autoantibody_deps()
# ============================================================================

# ---- CRAN packages ---------------------------------------------------------
.cran_pkgs <- c(
  "tidyverse", "dplyr", "tidyr", "readr", "readxl", "writexl", "tibble",
  "stringr", "stringi", "purrr", "rlang", "magrittr", "data.table",
  "ggplot2", "ggpubr", "ggrepel", "ggtext", "ggthemes", "ggprism", "ggbreak",
  "ggvenn", "ggalt", "VennDiagram", "patchwork", "aplot", "cowplot",
  "gridExtra", "hrbrthemes", "scales", "RColorBrewer", "viridis", "circlize",
  "pheatmap", "corrr", "broom", "rstatix", "PMCMRplus", "janitor", "vtable",
  "cluster", "vegan", "lsa", "Rtsne", "umap", "randomForest", "digest",
  "devtools", "Rcpp", "extrafont", "tm", "msigdbr", "gprofiler2", "flextable"
)

# ---- Bioconductor packages -------------------------------------------------
.bioc_pkgs <- c(
  "limma", "ComplexHeatmap", "SummarizedExperiment", "GSEABase",
  "clusterProfiler", "enrichplot", "DOSE", "fgsea", "org.Hs.eg.db",
  "pathview", "rWikiPathways", "HPAanalyze", "EnhancedVolcano", "Glimma",
  "oncoprint", "airway"
)

.all_pkgs <- c(.cran_pkgs, .bioc_pkgs)

# ---- Loader ----------------------------------------------------------------
#' Load all pipeline packages, warning (not erroring) on any that are missing.
load_autoantibody_libs <- function(quiet = TRUE) {
  missing <- character(0)
  for (p in .all_pkgs) {
    ok <- suppressWarnings(suppressPackageStartupMessages(
      requireNamespace(p, quietly = TRUE)
    ))
    if (ok) {
      suppressWarnings(suppressPackageStartupMessages(
        library(p, character.only = TRUE)
      ))
    } else {
      missing <- c(missing, p)
    }
  }
  if (length(missing)) {
    warning(
      "The following packages are not installed and were skipped:\n  ",
      paste(missing, collapse = ", "),
      "\nRun install_autoantibody_deps() to install them.",
      call. = FALSE
    )
  }
  invisible(missing)
}

# ---- Installer -------------------------------------------------------------
#' Install any missing CRAN / Bioconductor dependencies.
install_autoantibody_deps <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  inst <- rownames(installed.packages())
  need_cran <- setdiff(.cran_pkgs, inst)
  need_bioc <- setdiff(.bioc_pkgs, inst)
  if (length(need_cran)) install.packages(need_cran)
  if (length(need_bioc)) BiocManager::install(need_bioc, update = FALSE, ask = FALSE)
  invisible(list(cran = need_cran, bioc = need_bioc))
}

# Load on source() (matches the original behaviour of autoantibody_libraries.R)
load_autoantibody_libs()
