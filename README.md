# Autoantibody PhIP-seq / REAP analysis — melanoma cancer vaccine cohort

Analysis code for a manuscript profiling the **autoantibody (auto-Ab) repertoire**
of melanoma patients receiving a cancer vaccine, using **PhIP-seq** (phage
immunoprecipitation sequencing) and **REAP** data, compared against healthy
donors (HD). The pipeline covers QC and background-gene selection, differential
auto-Ab enrichment (limma), vaccine-induced kinetics and expression-pattern
clustering, pathway enrichment, tissue / subcellular-localisation annotation
(Human Protein Atlas), diversity and similarity metrics, and outcome / survival
modelling.

> This repository was reorganised from a working analysis folder for public
> release. The raw data are **not** included (see [Data](#data)); the code,
> reconstructed dependencies, a figure-source-data workbook, and example
> regenerated figures are.

## Repository layout

```
autoantibody-phipseq-melanoma/
├── README.md
├── .gitignore
├── R/
│   ├── setup.R                  # entry point: paths + loads libraries & colours
│   ├── libraries.R              # package loader / installer (reconstructed)
│   ├── color_tools.R            # figure palettes + IM_theme/QC_theme (reconstructed)
│   ├── 01_data_processing.Rmd           # metadata + PhIP-seq RPK processing, background selection
│   ├── 02_limma_baseline_melanoma_vs_HD.Rmd
│   ├── 03_limma_vaccine_induced_kinetics.Rmd   # kinetics, DEGpatterns clusters, pathway enrichment
│   ├── 04_limma_baseline_variantC.Rmd          # Good- vs Bad-outcome baseline
│   ├── 05_main_figures_partA.Rmd        # diversity, similarity, GSEA, tissue, regression, survival
│   ├── 06_main_figures_partB.Rmd        # day-0 / patient-only parallel analyses
│   ├── 07_tissue_expression_HPA.Rmd     # HPA tissue / subcellular-localisation figures
│   ├── gene_tissue_classification_source.R     # helper sourced by 07
│   ├── build_figure_source_data.R       # builds results/tables/Figure_Source_Data.xlsx
│   └── verify_figures.R                  # regenerates representative figures (repro check)
├── data/                        # inputs live here (git-ignored; see data/README.md)
└── results/
    ├── figures/                 # regenerated example figures (PNG)
    └── tables/
        └── Figure_Source_Data.xlsx   # one tab per figure of underlying data
```

## Pipeline order

1. **`01_data_processing`** — cleans clinical/vaccine metadata, imports the
   PhIP-seq RPK matrix (fixing genes mis-read as dates), splits amplification
   rounds 2/3, removes all-zero antigens, selects background vs "clear" antigens
   (bottom-5% expression + zero-count filter), and saves `PHIPseq_list.Rda`
   (the processed object list re-used by every downstream script).
2. **`02` / `04`** — limma differential auto-Ab enrichment: Melanoma vs HD (02)
   and Good vs Bad outcome (04) at baseline. Volcano + single-gene boxplots.
3. **`03`** — vaccine-induced kinetics; DEGpatterns clustering of expression
   trajectories; per-cluster pathway enrichment (GO / WikiPathways / Reactome)
   and GSEA; genes-in-pathways table.
4. **`05` / `06`** — main figures: repertoire diversity (Richness / Clonality /
   Shannon), Morisita-Horn & cosine similarity heatmaps, GSEA, DE heatmaps,
   HPA tissue mapping, linear-regression & Cox OS/PFS forest plots, Kaplan-Meier.
5. **`07`** — HPA tissue and subcellular-localisation figures (oncoprints,
   tissue-pattern heatmaps, localisation bars/tiles).

## Requirements

- **R ≥ 4.4** (developed/verified on R 4.5.2).
- Packages are declared and loaded by `R/libraries.R` (CRAN + Bioconductor).
  Install everything with:

  ```r
  source("R/libraries.R")
  install_autoantibody_deps()
  ```

  Notes: `ggalt` is not available for R ≥ 4.5 (used only for two density plots;
  a `geom_density` fallback works). `airway` / `oncoprint` are not required by
  the figure code (the oncoprints use `ComplexHeatmap::oncoPrint`).

## Running

Every script starts with `source("R/setup.R")`, which locates the project root,
defines `data_path` / `results_path`, and loads the package stack + figure
colours/themes. Then either knit an `.Rmd` in RStudio, or from the project root:

```r
setwd("R"); source("setup.R")   # loads everything
# then run chunks of the desired NN_*.Rmd
```

Point the pipeline at your data either by placing files in `data/` or by
setting an environment variable before launching R:

```bash
export AUTOANTIBODY_DATA=/path/to/your/data
```

## Reconstructed dependencies

The original pipeline sourced two machine-specific files that were lost
(`autoantibody_libraries.R`, `LHB_color_tools.R`). Both have been rebuilt:

- **`R/libraries.R`** — from the full set of `library()` calls across the
  pipeline; loads defensively (missing packages warn, don't abort).
- **`R/color_tools.R`** — from the exact colour/theme objects recovered from the
  saved workspace: `colors_LHB`, `custom_colors`, the subcellular `col` key,
  `IM_theme` / `QC_theme` (verified: identical except x-axis label angle 0° vs
  90°), the jitter positions and the gene-string wrapping helpers.

All original hard-coded absolute paths were replaced by project-relative paths.

## Figure source data

`results/tables/Figure_Source_Data.xlsx` provides, **one worksheet per figure**,
the underlying data table for that figure, plus a `00_Figure_Index` master sheet
mapping all 40+ manuscript figures to their source object and whether the data
is included or must be regenerated by running the producing script. Rebuild it
with:

```bash
Rscript R/build_figure_source_data.R
```

Figures whose plotting object was computed on the fly inside a script (volcano
top-tables, diversity/similarity matrices, regression & Cox tables, on-the-fly
heatmap matrices) are flagged `NO`/`PARTIAL` in the index with the script to run;
figures backed by saved objects (QC, cluster profiles, pathway enrichment,
localisation, per-gene expression) are included directly.

## Reproducibility check

`R/verify_figures.R` regenerates representative figures (expression-distribution
histogram, DEGpatterns cluster profiles, pathway-enrichment dot plot, single-gene
boxplot, subcellular-localisation bars) with the exact script code and the saved
data, confirming the reconstructed environment renders figures end-to-end.
Example outputs are in `results/figures/`.

## Data

Not distributed here — several raw inputs exceed GitHub's 100 MB/file limit
(largest ~720 MB). See **[`data/README.md`](data/README.md)** for the full input
manifest, sizes, which script uses each, and how to obtain them (the Human
Protein Atlas TSVs are public).
