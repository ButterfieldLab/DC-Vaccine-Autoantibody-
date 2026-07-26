# ============================================================================
# build_figure_source_data.R
# ----------------------------------------------------------------------------
# Builds `results/tables/Figure_Source_Data.xlsx`: one worksheet per manuscript
# figure (or figure group), each holding the underlying data for that figure.
#
# Data are extracted from the saved analysis workspace (.RData) and the
# processed object list (PHIPseq_list.Rda). Figures whose plotting objects are
# computed on-the-fly inside the analysis scripts (and were not saved) are
# listed in the `00_Figure_Index` sheet with the script + object needed to
# regenerate them.
#
# Usage:
#   Rscript R/build_figure_source_data.R  /path/to/.RData  /path/to/output.xlsx
# Defaults point at the original project workspace.
# ============================================================================

suppressWarnings(suppressPackageStartupMessages({
  library(writexl); library(dplyr); library(tibble)
}))

args <- commandArgs(trailingOnly = TRUE)

# Locate the project root (env override -> walk up from wd -> wd).
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
root <- find_root()

# Saved workspace: arg 1 -> AUTOANTIBODY_RDATA -> data/.RData under the project.
rdata_path <- if (length(args) >= 1) args[[1]] else
  Sys.getenv("AUTOANTIBODY_RDATA", unset = file.path(root, "data", ".RData"))
out_path <- if (length(args) >= 2) args[[2]] else
  file.path(root, "results", "tables", "Figure_Source_Data.xlsx")

if (!file.exists(rdata_path))
  stop("Saved workspace not found: ", rdata_path,
       "\nPass it as arg 1, or set AUTOANTIBODY_RDATA, or place it at data/.RData")

message("Loading workspace: ", rdata_path)
e <- new.env(); load(rdata_path, envir = e)
g  <- function(x) get(x, envir = e)
has <- function(x) exists(x, envir = e)

PL <- g("PHIPseq_list")

# ---- helper: coerce anything to a tidy data.frame -------------------------
as_df <- function(x) {
  if (is.data.frame(x)) return(as.data.frame(x))
  if (is.matrix(x)) return(tibble::rownames_to_column(as.data.frame(x), "rowname"))
  as.data.frame(x)
}

sheets <- list()
add <- function(name, df) {
  df <- as_df(df)
  sheets[[substr(name, 1, 31)]] <<- df
  message(sprintf("  + %-31s %d x %d", substr(name,1,31), nrow(df), ncol(df)))
}

# ---- Figure tabs (data recovered from the saved workspace) -----------------

# Fig 1A : QC / background gene selection (scatter + histogram) -> Row_sums_rpk
add("Fig1A_QC_background", PL$Row_sums_rpk)

# Fig 1B : DEGpatterns cluster expression profiles (all clusters, Good/Bad)
add("Fig1B_cluster_profiles", PL$clusters_vacc.cycle_outcome)

# Fig 1B : pathway enrichment dot plots, per cluster (combined across clusters)
if (has("combined_df")) add("Fig1B_pathways_by_cluster", g("combined_df"))

# Fig 1B : genes-in-pathways table-figure (the flextable saved as flextable.png)
if (has("df")) add("Fig1B_pathwayGenes_table", g("df"))

# Fig 1B : cluster-1 enrichment results (Reactome / GO / WikiPathways)
if (has("result_list")) {
  rl <- g("result_list")
  if (!is.null(rl$enrichPathway_res)) add("Fig1B_enrichReactome", rl$enrichPathway_res)
  if (!is.null(rl$enrichGO_res))      add("Fig1B_enrichGO", rl$enrichGO_res)
  if (!is.null(rl$enrichWP_res))      add("Fig1B_enrichWikiPathways", rl$enrichWP_res)
}

# Fig Tissue : selected cluster pathways feeding the oncoprint
if (has("CLUST.combined.df")) add("FigTIS_cluster_pathways", g("CLUST.combined.df"))

# Fig Tissue : gene -> cluster -> subcellular Category (oncoprint row annotation)
if (has("abc")) add("FigTIS_gene_localization", g("abc"))

# Fig Tissue : oncoprint matrix (cluster x gene presence)
if (has("aa_meta")) add("FigTIS_oncoprint_matrix", g("aa_meta"))

# Fig Tissue : cluster genes with subcellular location + Category
if (has("CLUST.tab.tissue")) add("FigTIS_cluster_gene_tissue", g("CLUST.tab.tissue"))

# Reference tables used to classify HPA terms into categories
if (has("cell.location.tbl")) add("Ref_subcellular_classes", g("cell.location.tbl"))
if (has("anatomical_categories")) add("Ref_HPA_anatomical_cats", g("anatomical_categories"))
if (has("cell_structure_categories")) add("Ref_HPA_cellstruct_cats", g("cell_structure_categories"))

# ---- Single-gene boxplot figures: per-gene RPK expression ------------------
# (1BA CAMKV, 1BB TDP2/PSMA5, 1BC MR1, 2A STMN4, 2B ZPR1, survival genes, etc.)
boxplot_genes <- c("CAMKV","TDP2","PSMA5","MR1","STMN4","ZPR1",
                   "KRT78","ARTN","RBM17","BAK1","VHL")
long <- PL$R3_rpk_long
keep_cols <- intersect(
  c("Gene","Sample","Subject_ID","Expression","Vac_cycle","Cancer",
    "Outcome","Outcome.1","Response","Best_Clinical_outcome",
    "Overall_Survival","PFS"),
  names(long))
box_df <- long[long$Gene %in% boxplot_genes, keep_cols, drop = FALSE]
box_df <- box_df[order(box_df$Gene, box_df$Subject_ID), ]
add("Fig_boxplot_genes_RPK", box_df)

# ---- 00_Figure_Index : master map of every manuscript figure ---------------
idx <- tibble::tribble(
  ~Figure, ~Script, ~Description, ~Data_object, ~In_workbook, ~Tab_or_note,
  "1A-1..4","01_data_processing","QC scatter + histogram, background-gene selection","Row_sums_rpk","YES","Fig1A_QC_background",
  "1BA-1","02_limma_baseline_melanoma_vs_HD","Volcano Melanoma vs HD baseline","toptab / fin.volcano.tab","NO","Run 02 (limma topTable)",
  "1BA-2/3","02_limma_baseline_melanoma_vs_HD","Boxplot CAMKV by Cancer / over cycle","R3_rpk_long (CAMKV)","YES","Fig_boxplot_genes_RPK",
  "1BA-4","02_limma_baseline_melanoma_vs_HD","Top-gene boxplot panels","gene_data (R3_rpk_long)","PARTIAL","Fig_boxplot_genes_RPK (subset)",
  "1BA-5","02_limma_baseline_melanoma_vs_HD","Melanoma vs HD row-sum scatter","my_df_plot","NO","Run 02",
  "1BB-1","03_limma_vaccine_induced_kinetics","Volcano baseline patient contrast","top_tables_list[[data_set]]","NO","Run 03",
  "1BB-2/3","03_limma_vaccine_induced_kinetics","Boxplot TDP2/PSMA5 over cycle","R3_rpk_long","YES","Fig_boxplot_genes_RPK",
  "1BB-4/9/18","03_limma_vaccine_induced_kinetics","DEGpatterns cluster expression profiles","clust.tab / clusters_vacc.cycle_outcome","YES","Fig1B_cluster_profiles",
  "1BB-5","03_limma_vaccine_induced_kinetics","enrichGO dot plot (cluster 1)","enrichGO_results[[1]]","YES","Fig1B_enrichGO",
  "1BB-6","03_limma_vaccine_induced_kinetics","enrichWP dot plot (cluster 1)","enrichWP_results[[1]]","YES","Fig1B_enrichWikiPathways",
  "1BB-7","03_limma_vaccine_induced_kinetics","Selected cluster pathways dot plot","combined_df","YES","Fig1B_pathways_by_cluster",
  "1BB-8","03_limma_vaccine_induced_kinetics","Genes-in-pathways table (flextable.png)","df","YES","Fig1B_pathwayGenes_table",
  "1BB-10..13","03_limma_vaccine_induced_kinetics","GSEA lollipop / fold-enrichment","ego_summary / out_gsea@result","NO","Run 03 (GSEA)",
  "1BB-14..18","03_limma_vaccine_induced_kinetics","Up/down counts + expression per cluster","abc_wide / abc / clust.tab","PARTIAL","Fig1B_cluster_profiles",
  "1BC-1","04_limma_baseline_variantC","Volcano Good vs Bad outcome","top_tables_list[[a]]","NO","Run 04",
  "1BC-2/3","04_limma_baseline_variantC","Boxplot MR1 by outcome / cycle","R3_rpk_long_wo_HDs / R3_rpk_long","YES","Fig_boxplot_genes_RPK",
  "1BC-4","04_limma_baseline_variantC","Top-down gene boxplot panels","gene_data","PARTIAL","Fig_boxplot_genes_RPK (subset)",
  "1BC-5/6","04_limma_baseline_variantC","Background histogram + Mel-vs-HD scatter","df_test","NO","Run 04",
  "2A-1","05_main_figures_partA","Regulated auto-Ab volcano","d89_raw_stats_filt","NO","Run 05",
  "2A-2 / 2B-5","05/06_main_figures","Sample-distance heatmap","sampleDistMatrix","NO","Run 05/06",
  "2A-3 / 2B-6","05/06_main_figures","PCA PC1 vs PC2","PC","NO","Run 05/06",
  "2A-4..8 / 2B-7","05/06_main_figures","Diversity: Richness/Clonality/Shannon","Phip.seq.div","NO","Run 05/06 (vegan)",
  "2A-6/11 / 2B-11/12","05/06_main_figures","Antibody-quantile abundance","Phip_seq_median","NO","Run 05/06",
  "2A-9/10 / 2B-8/9/10","05/06_main_figures","Clonotype counts by quantile","Phip_seq_final","NO","Run 05/06",
  "2A-12/13 / 2B-13/14","05/06_main_figures","Morisita-Horn / cosine similarity heatmaps","Phip_mx_mat / Phip_cos_mat","NO","Run 05/06 (lsa/divo)",
  "2A-14/15 / 2B-15/16","05/06_main_figures","Cosine index over cycle / at day0","cos_combined_df / cos_day_0_df","NO","Run 05/06",
  "2A-16","05_main_figures_partA","Volcano HD vs Melanoma (limma)","volcano.top.table","NO","Run 05",
  "2A-17/18 / 2B-4","05/06_main_figures","Boxplot STMN4 / ZPR1","R3_rpk_long","YES","Fig_boxplot_genes_RPK",
  "2A-19/22 / 2A-23/24","05_main_figures_partA","GSEA lollipop / dotplot / enrichment curve","my_gsea_list@result / gse","NO","Run 05 (fgsea/DOSE)",
  "2A-20 / 2B-22","05/06_main_figures","Heatmap of top DE / significant proteins","sub_df / d0 (matrix_phip)","NO","Run 05/06",
  "2A-21","05_main_figures_partA","Tissue-expression counts bar","tab_x","NO","Run 05 (HPA)",
  "2A-25/31 2B-17/23","05/06_main_figures","Forest plots linear regression Outcome","Plot_Outcome.1 / Plot_Outcome.2","NO","Run 05/06",
  "2A-26/32 2B-18/24","05/06_main_figures","Estimate vs -log2(padj) volcano","Adj.p_Modelunnest","NO","Run 05/06",
  "2A-27/28/34 2B-19/20/26","05/06_main_figures","Protein-of-interest point/line plots","Phip_seq_long_df","NO","Run 05/06",
  "2A-33 / 2B-25","05/06_main_figures","Richness vs best clinical outcome","Phip.seq.0h.div","NO","Run 05/06",
  "2A-35/36 2B-27/28/30","05/06_main_figures","Cox OS / PFS forest plots","Unnest_Model_Cox_OS / _PFS","NO","Run 05/06 (survival)",
  "2A-37 / 2B-29","05/06_main_figures","Kaplan-Meier PFS curve","res.PFS (Phip_fit)","NO","Run 05/06 (survminer)",
  "2A-38 / 2B-31","05/06_main_figures","Baseline Melanoma vs HD volcano","raw_statistics_df","NO","Run 05/06",
  "2B-1","06_main_figures_partB","Density before/after filtering","filt_Phip_rpk_long","NO","Run 06",
  "2B-2/3","06_main_figures_partB","Volcano + interactive (Glimma)","allDEresults","NO","Run 06",
  "TIS-1","07_tissue_expression_HPA","Oncoprint cluster genes x tissue/cell category","aa_meta / CLUST.combined.df","YES","FigTIS_oncoprint_matrix + FigTIS_cluster_pathways",
  "TIS-2","07_tissue_expression_HPA","Paired oncoprints Melanoma vs HD","mel_pts_mat / mel_HD_mat","NO","Run 07",
  "TIS-3/7","07_tissue_expression_HPA","Tissue-expression pattern heatmaps","my.df_mat / tissue_toptab_wide","NO","Run 07",
  "TIS-4/5","07_tissue_expression_HPA","Tissue-expression bar / tile heatmap","tissue_toptab","NO","Run 07",
  "TIS-6","07_tissue_expression_HPA","Patient heatmap of tissue-annotated genes","sub_df_mat (matrix_phip)","NO","Run 07",
  "TIS (ref)","07_tissue_expression_HPA","Gene subcellular classification / cluster tissue","cell.location.tbl / CLUST.tab.tissue / abc","YES","Ref_subcellular_classes / FigTIS_*"
)

# Put the index first
sheets <- c(list("00_Figure_Index" = as.data.frame(idx)), sheets)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(sheets, path = out_path)
message("\nWrote ", length(sheets), " sheets to:\n  ", out_path)
message("Sheets: ", paste(names(sheets), collapse = ", "))
