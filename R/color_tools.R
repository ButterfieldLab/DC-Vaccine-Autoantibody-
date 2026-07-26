# ============================================================================
# color_tools.R
# ----------------------------------------------------------------------------
# Colour palettes, ggplot themes and small text helpers for the manuscript
# figures. Reconstructed from the original (lost) `LHB_color_tools.R` using the
# exact objects recovered from the saved workspace (.RData): colours, the
# IM_theme / QC_theme definitions, jitter positions and the gene-wrapping
# helpers. Requires ggplot2 (loaded via R/libraries.R).
# ============================================================================

# ---- Named colour scales (used by scale_*_manual in the figures) -----------
colors_LHB <- list(
  Response = c(HD = "magenta", PR = "red", SD = "green",
               NED1 = "lightblue", NED2 = "blue", PD = "grey45"),
  Outcome  = c(HD = "magenta", Good = "red", Bad = "grey45"),
  Outcome.1 = c(HD = "magenta", Good = "red", Stable = "green", Bad = "grey45"),
  Vac_cycle = c(day_0 = "grey65", day_43 = "green2", day_89_101 = "orange1"),
  Cancer1  = c(HD = "#ADD8E6", Melanoma = "#FFB6C1"),
  Cancer2  = c(Healthy = "#ADD8E6", Cancer = "#FFB6C1"),
  Out_median_lines = c(Good = "red3", Bad = "grey10")
)

# ---- Discrete qualitative palettes -----------------------------------------
colors_dutch <- c(
  "#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67", "#F79F1F", "#A3CB38",
  "#1289A7", "#D980FA", "#B53471", "#EE5A24", "#009432", "#0652DD", "#9980FA",
  "#833471", "#EA2027", "#006266", "#1B1464", "#5758BB", "#6F1E51"
)

colors_spanish <- c(
  "#ff5252", "#ff793f", "#d1ccc0", "#ffb142", "#ffda79", "#2c2c54", "#474787",
  "#aaa69d", "#227093", "#218c74", "#b33939", "#cd6133", "#84817a", "#cc8e35",
  "#ccae62", "#40407a", "#706fd3", "#f7f1e3", "#34ace0", "#33d9b2", "#4DA167",
  "#947EB0", "#083D77", "#06707E"
)

cell_type_cols <- c(
  "#FFC312", "#C4E538", "#006266", "#FDA7DF", "#ED4C67", "#009432", "#0652DD",
  "#F79F1F", "#A3CB38", "#1289A7", "#D980FA", "#0652DD", "#ff5252", "#1B1464",
  "#9980FA", "#833471", "#B53471", "#5758BB", "#6F1E51", "#ccae62", "#84817a",
  "#ffda79", "#ff793f", "#d1ccc0", "#34ace0", "#33d9b2", "#4DA167", "#947EB0",
  "#083D77", "#aaa69d"
)

custom_colors <- list(
  discrete = c(colors_dutch, colors_spanish),
  colors_dutch = colors_dutch,
  colors_spanish = colors_spanish,
  cell_type_cols = cell_type_cols
)

# Subcellular-location colour key used in the tissue / localisation figures
col <- c(
  "Cytoplasm and Associated Structures" = "#1f78b4",
  "Nucleus and Nuclear Structures"      = "#33a02c",
  "Organelles"                          = "#e31a1c",
  "Membrane-associated Structures"      = "#ff7f00",
  "Cytoskeleton"                        = "#6a3d9a"
)

cluster_labels <- setNames(paste("Cluster", 1:12), as.character(1:12))

# ---- ggplot themes ---------------------------------------------------------
# IM_theme: manuscript "immune" figure theme (horizontal x labels).
IM_theme <- ggplot2::theme(
  axis.title   = ggplot2::element_text(family = "Arial", colour = "black", size = 18),
  axis.title.x = ggplot2::element_text(colour = "black", size = 20),
  axis.title.y = ggplot2::element_text(colour = "black", size = 20),
  axis.text    = ggplot2::element_text(colour = "black", size = 15),
  axis.text.x  = ggplot2::element_text(hjust = 1, vjust = 1, angle = 0),
  axis.ticks   = ggplot2::element_line(colour = "black"),
  axis.ticks.length = ggplot2::unit(0.4, "lines"),
  axis.line    = ggplot2::element_line(colour = "black"),
  legend.key   = ggplot2::element_blank(),
  legend.text  = ggplot2::element_text(family = "Arial", size = 15),
  legend.title = ggplot2::element_blank(),
  panel.background = ggplot2::element_rect(fill = "white", colour = "black"),
  panel.grid   = ggplot2::element_blank(),
  strip.text   = ggplot2::element_text(face = "bold", size = 12)
)

# QC_theme: identical to IM_theme but with 90-degree rotated x labels
# (used for the QC / many-category plots).
QC_theme <- IM_theme + ggplot2::theme(
  axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, angle = 90)
)

# ---- Jitter positions (fixed seeds for reproducible point spread) ----------
jitter_0 <- ggplot2::position_jitter(width = 0.0, height = 0, seed = 1)
jitter_1 <- ggplot2::position_jitter(width = 0.1, height = 0, seed = 1)
jitter_2 <- ggplot2::position_jitter(width = 0.2, height = 0, seed = 1)

# ---- Text helpers: wrap long gene strings across lines ---------------------
add_newline_after_5 <- function(gene_string) {
  genes <- unlist(strsplit(gene_string, " "))
  if (length(genes) > 3) {
    paste(sapply(seq(1, length(genes), by = 3), function(i) {
      paste(genes[i:min(i + 9, length(genes))], collapse = " ")
    }), collapse = "\n")
  } else {
    gene_string
  }
}

add_newline_after_10 <- function(gene_string) {
  genes <- unlist(strsplit(gene_string, " "))
  if (length(genes) > 10) {
    paste(sapply(seq(1, length(genes), by = 10), function(i) {
      paste(genes[i:min(i + 9, length(genes))], collapse = " ")
    }), collapse = "\n")
  } else {
    gene_string
  }
}
