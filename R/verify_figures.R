# ============================================================================
# verify_figures.R
# ----------------------------------------------------------------------------
# Reproducibility check: regenerates a representative set of manuscript figures
# using the EXACT plotting code from the analysis scripts, driven by the saved
# processed data objects. Confirms that (a) the reconstructed package stack and
# colour/theme tools load, and (b) the figure code renders against the data.
#
# Outputs PNGs to results/figures/. Run from the project root:
#   Rscript R/verify_figures.R
# ============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(magrittr); library(forcats)
  library(ggprism); library(tibble)
}))

find_root <- function() {
  env <- Sys.getenv("AUTOANTIBODY_ROOT", unset = NA)
  if (!is.na(env) && dir.exists(env)) return(normalizePath(env))
  d <- normalizePath(getwd())
  for (i in 1:6) {
    if (file.exists(file.path(d, "R", "setup.R"))) return(d)
    if (dirname(d) == d) break; d <- dirname(d)
  }
  normalizePath(getwd())
}
root  <- find_root()
rdata <- Sys.getenv("AUTOANTIBODY_RDATA", unset = file.path(root, "data", ".RData"))
if (!file.exists(rdata))
  stop("Saved workspace not found: ", rdata,
       "\nSet AUTOANTIBODY_RDATA or place it at data/.RData")
figdir <- file.path(root, "results", "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

source(file.path(root, "R", "color_tools.R"))   # colours + IM_theme/QC_theme
e <- new.env(); load(rdata, envir = e); g <- function(x) get(x, envir = e)
PL <- g("PHIPseq_list")

ok <- character(0); fail <- character(0)
render <- function(id, file, expr) {
  p <- tryCatch(force(expr), error = function(err) {
    message("  x ", id, " : ", conditionMessage(err)); fail[[length(fail)+1]] <<- id; NULL })
  if (is.null(p)) return(invisible())
  tryCatch({
    ggsave(file.path(figdir, file), p, width = 8, height = 6, dpi = 150)
    message("  ok ", id, " -> ", file); ok[[length(ok)+1]] <<- id
  }, error = function(err) { message("  x save ", id, " : ", conditionMessage(err)); fail[[length(fail)+1]] <<- id })
}

# ---- Fig 1A : PhIP-seq expression distribution (QC / background) -----------
Row_sums_rpk <- PL$Row_sums_rpk
render("1A histogram", "fig_1A_expression_distribution.png",
  Row_sums_rpk %>%
    mutate(Zero_counts = ifelse(All_0_counts < 105, "Low", "High")) %>%
    ggplot(aes(x = log2(Row_sums_All))) +
    geom_histogram(bins = 250, aes(fill = Zero_counts)) +
    scale_fill_manual(values = c("red", "black")) +
    theme_prism() +
    labs(title = "Distribution of PhIP-seq Expression Values",
         x = "Log2(rpk)", y = "count") +
    theme(aspect.ratio = 1/1) +
    facet_wrap(~Background, scales = "free_x"))

# ---- Fig 1B (1BB-4) : DEGpatterns cluster expression profiles --------------
clust.tab <- PL$clusters_vacc.cycle_outcome
clust.tab$cluster <- as.factor(clust.tab$cluster)
clust.tab %<>% mutate(Vac_cycle = gsub("day_", "", Vac_cycle),
                      Vac_cycle = gsub("89_", "", Vac_cycle))
clust.tab$Vac_cycle <- factor(clust.tab$Vac_cycle, levels = c("0", "43", "101"))
cluster_labels <- setNames(paste("Cluster", unique(clust.tab$cluster)),
                           unique(clust.tab$cluster))
render("1BB-4 cluster profiles", "fig_1B_cluster_profiles.png",
  ggplot(clust.tab, aes(Vac_cycle, value, color = Outcome)) +
    geom_line(data = clust.tab %>% filter(Outcome == "Good"),
              aes(group = genes), color = "red", size = 0.1, alpha = 0.5) +
    geom_line(data = clust.tab %>% filter(Outcome == "Bad"),
              aes(group = genes), color = "grey45", size = 0.1, alpha = 0.5) +
    scale_color_manual(values = c(colors_LHB$Out_median_lines)) +
    theme_prism() +
    geom_smooth(aes(group = Outcome), method = "loess") +
    labs(title = "Auto-antibody Expression Profiles",
         x = "Vaccine cycle", y = "Scaled(rpk)") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 14, face = "bold"),
          legend.position = "right") +
    facet_wrap(~cluster, labeller = labeller(cluster = cluster_labels)))

# ---- Fig 1B (1BB-5) : pathway enrichment dot plot (enrichGO cluster 1) -----
enrichGO_res <- g("result_list")$enrichGO_res
render("1BB-5 pathway dotplot", "fig_1B_pathway_enrichment_dotplot.png",
  enrichGO_res %>%
    arrange(desc(Count)) %>% filter(pvalue <= 0.05) %>% filter(Count >= 3) %>%
    head(30) %>%
    ggplot(aes(x = -log2(pvalue), y = fct_reorder(Description, -log2(pvalue)))) +
    geom_point(aes(size = Count, color = as.factor(Count))) +
    theme_bw() +
    labs(title = "Pathway Enrichment", x = "-log2(pvalue)", y = "Pathway") +
    scale_size_continuous(range = c(0, 9)) +
    theme(legend.position = "right"))

# ---- Fig (boxplot) : single-gene RPK over vaccine cycle --------------------
gene_int <- "PSMA5"
render("single-gene boxplot", "fig_boxplot_single_gene.png",
  PL$R3_rpk_long %>%
    filter(Outcome != "HD") %>%
    dplyr::filter(Gene %in% gene_int) %>%
    ggplot(aes(x = Vac_cycle, y = log(Expression))) +
    geom_boxplot() +
    geom_line(aes(group = Subject_ID), color = "Grey25", alpha = 0.3) +
    geom_point(aes(col = Response)) +
    scale_color_manual(values = c(colors_LHB$Response)) +
    theme_bw() + theme(aspect.ratio = 1) +
    ggtitle(paste0("Expression of ", gene_int, " gene")) +
    facet_wrap(~Outcome + Gene))

# ---- Fig Tissue : subcellular-location classification of cluster genes -----
render("tissue gene localization", "fig_tissue_gene_localization.png",
  g("abc") %>% count(cluster, Category) %>%
    ggplot(aes(x = as.factor(cluster), y = n, fill = Category)) +
    geom_col() +
    scale_fill_manual(values = col) +
    theme_prism() +
    labs(title = "Subcellular localisation of cluster auto-antigens",
         x = "Cluster", y = "Gene count") +
    theme(legend.position = "right"))

message("\n==== verify_figures summary ====")
message("Rendered OK (", length(ok), "): ", paste(ok, collapse = " | "))
if (length(fail)) message("FAILED (", length(fail), "): ", paste(fail, collapse = " | "))
