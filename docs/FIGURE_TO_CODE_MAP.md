# Figure → Source-code map

Maps every panel in `figures pdf.pdf` (6 main figures, 9 supplemental figures,
4 supplemental tables) to the code block that generates it.

**How to read line refs:** `v3:2865` = `REAP Butterfield v3.R` line 2865.
Files referenced:
- `v3` = `REAP Butterfield v3.R` — **the definitive REAP figure script** (superset; carries the author's own `FIG…` comments). *Currently in the old project folder, not yet copied into this repo — see note at bottom.*
- `v2` = `R/REAP Butterfield v2.R` — polished subset of the REAP panels.
- `v1` = `R/REAP Butterfield .R` — older: tSNE clinical panels, REAP hit counts, immune-autoantibody pheatmap.
- `NB 0x` = the PhIP-seq R-notebook pipeline in `R/` (`01`–`07`).

> ⚠️ **The `FIG1A` / `FIG1C` comments inside the REAP code are out of date.**
> They reflect an earlier panel arrangement. The mapping below is against the
> **final** panel letters in `figures pdf.pdf`.

---

## Which file produces which figure (overview)

| Analysis | Figures | Source |
|---|---|---|
| **REAP autoantibody screen** | Fig 1, Fig 2A/B/D, Fig 3A/C, Supp 1, Supp 6, Supp 9, Supp Table 2 | `v3` (primary), `v2`, `v1` |
| **REAP cytokine assay** | Fig 2C | `v3` only |
| **PhIP-seq epitope-spreading** | Fig 3B/D/E | `v3` only |
| **PhIP-seq outcome pipeline** | Fig 5, Supp 3/4/5/8 | `NB 05`/`NB 06` (+ `NB 02`–`04`) |
| **PhIP-seq fragment/epitope** | Fig 6 | PhIP-seq peptide-level (see notes) |
| **DC CyTOF / metabolism correlations** | Fig 4, Supp 7 | *Not in the provided files* |

---

## Main figures

### Figure 1 — Autoantibodies are increased in melanoma patients  *(REAP)*
| Panel | What it shows | Code |
|---|---|---|
| **1A** | REAP-score heatmap, HD vs Melanoma, proteins grouped by functional category, columns split by prior checkpoint | `v3:508` / `v2:499` (`heatmap_baseline_melanoma_vs_hd.png`) — ComplexHeatmap on `heat_map_data` |
| **1B** | Violin plots "No. of Reactivities (score ≥1)", HD vs Melanoma, per pathway category | `v3:709` / `v2:698` (`pathway_hits_graph.png`) — *code comment says "FIG1C"* |
| **1C** | Individual-protein REAP dot plots: IL15RA, DEFB118-N, APOH, FGFBP3, ASIP, APOC3 | `v3:928` (titles at `v3:782–917`) / `v2:917` (`individual_protein_hd_vs_melanoma_graph.png`) |

### Figure 2 — Autoantibody profiling & vaccine clinical efficacy  *(REAP)*
| Panel | What it shows | Code |
|---|---|---|
| **2A** | Baseline heatmaps stacked "Up in NRs" over "Up in Rs" (HD/R/NR) | `v3:2206` / `v2:2193` (`heatmap_baseline_melanoma_R_vs_NR.png`, **live save**) |
| **2B** | CXCL9 & CXCR3 baseline REAP score, HD/NR/R, patient-labeled | Responder per-protein jitter (same pattern as `responder_checkpoint_graph`, `v3:~2320–2372`); **no dedicated `ggsave` name** — block filters `Protein %in% c("CXCL9","CXCR3")` |
| **2C** | IL-6 (pg/mL) over Day 0/43/89, NR vs R | `v3:4381` (`cytokine_graph.png`, **live save**; cytokine section `v3:4300–4417`) — **v3 only** |
| **2D** | Anti-CTLA-4 (NRs) & Anti-PD-1 (Rs) REAP over time | `v3:2373` (`checkpoints_over_time_graph.png`) — **v3 only** |

### Figure 3 — REAP detects B-cell epitope spreading  *(REAP + PhIP-seq RPK)*
| Panel | What it shows | Code |
|---|---|---|
| **3A** | Vaccine proteins RPK over time: MAGE-A6, MART-1 (MLANA), Tyrosinase (TYR) | `v3:2821` / `v2:2744` (`vaccine_responses.png`) |
| **3B** | Non-vaccine antigens: gp100, NY-ESO-1, MLPH | `v3:2865` (`gp100_nyeso_mlph_graph.png`) — **v3 only** |
| **3C** | Related MAGE: MAGEA10, MAGEA4, MAGEB1, MAGEB3, MAGEC1 | `v3:2935` / `v2:2810` (`mage_graphs.png`) |
| **3D** | Cancer-testis-antigen de-novo autoAb count, Day 43/89, NR vs R | `v3:3280` (`cta_hits_graphs_ifn.png`) — **v3 only** |
| **3E** | All de-novo autoantibodies, Day 43/89, NR vs R | `v3:3424` (`denovo_hits_graphs_ifn.png`) — **v3 only** |

### Figure 4 — cDC2 frequencies & glycolytic DC metabolism  *(NOT REAP)*
Spearman-correlation bar panels (auto-Ab vs peripheral-population frequencies /
metabolic mDC profiles). **No generating code in v1/v2/v3** — this is a separate
DC CyTOF / metabolism correlation analysis not among the files provided.

### Figure 5 — Autoantibody signatures associate with clinical outcome  *(PhIP-seq pipeline)*
| Panel | What it shows | Code (this repo) |
|---|---|---|
| **5A** | Auto-Ab-targeted protein categories, NR vs R (bar) | `NB 05` limma + category annotation |
| **5B** | NR / R auto-Ab pathway enrichment (dot) | `NB 05` GSEA/enrich sections |
| **5C** | Cox-regression OS / PFS forest plots (Hazard Ratio) | `NB 05` `Unnest_Model_Cox_OS` / `_PFS` (figs 2A-35/36) |
| **5D** | Auto-Ab categories with increased HR (OS/PFS bar) | `NB 05` (derived from Cox results) |
| **5E–G** | Cosine-similarity index over vaccine cycle (by outcome) | `NB 05` `cos_combined_df` (fig 2A-14) |
| **5H** | PCA of auto-Ab signature | `NB 05` `PC` (fig 2A-3) |
| **5I** | GSVA scores, Good vs Bad outcome | `NB 05`/`NB 06` GSVA block |

### Figure 6 — Autoantibody specificities distinguish response  *(PhIP-seq fragment-level)*
Fragment-resolution KIR2DL5B heatmap + per-protein RPK timecourses (ADAM12, CD27,
EGF, NCAM1, APP, MUC16, RNF169) + AlphaFold structure models. Plot style matches
the `v3` induced-protein timecourses (`v3:2660–3020`) but at **peptide/fragment**
resolution (uses the peptide-level PhIP-seq matrix, not the gene-collapsed one).
Structure ribbons are external (AlphaFold/PyMOL), not R.

---

## Supplemental figures & tables

| Item | Title | Source |
|---|---|---|
| **Supp Fig 1** | REAP scores per patient/timepoint (jitter + mean±SD, faceted PT_1…PT_35, HD1…HD5) | `v1:479` (`stat_summary` mean±SD on `Paul_patient_tidy`, faceted by patient) |
| **Supp Fig 2** | (formerly supp fig 4) | *see page 9 — mixed panels; confirm intent* |
| **Supp Fig 3** | RPK scores among all patients/timepoints (heatmap) | `NB` PhIP-seq (peptide/gene RPK heatmap) |
| **Supp Fig 4** | PhIP-seq melanoma unique autoAb repertoire (DEGs + pathway enrichment) | `NB 02`/`NB 05` (limma DEGs, enrichment) |
| **Supp Fig 5** | Increased in Melanoma vs HD, Log2(RPK) | `NB 05` baseline volcano / category counts |
| **Supp Fig 6** | Humoral response to other MAGE A/B/C proteins | `v3:2935` (`mage_graphs.png` family) |
| **Supp Fig 7** | Immune correlates of autoAb production (formerly fig 4) | DC-correlation analysis (not in these files) |
| **Supp Fig 8** | AutoAb concentrations & OS/PFS associations | `NB 05` Cox OS/PFS (+ REAP concentration) |
| **Supp Fig 9** | Machine learning (XGBoost feature importance, bootstrap ROC AUC=0.929, PCA, auto-Ab heatmap, top-feature Log2(RPK) boxplots) | **PhIP-seq** XGBoost ML script (uses Log2(RPK) auto-Ab signature) — *not the REAP randomForest/tSNE at `v3:1112`/`v1` tSNE, which is a separate unused exploration* |
| **Supp Table 1** | Dates of checkpoint administration | metadata table (not plotted) |
| **Supp Table 2** | Induced autoantibodies by REAP at Day 89 | `v3:3443` `write_xlsx("Autoatnbiody hits at day 89.xlsx")` |
| **Supp Table 3** | Proteins in functional categories & pathway analyses | `v3` annotation section (`v3:~398`) / `annotation_df` |
| **Supp Table 4** | AA sequences unique to Rs/NRs | peptide-sequence extraction (PhIP-seq library) |

---

## Notes & recommendations

1. **Add `REAP Butterfield v3.R` to the repo.** It is the definitive REAP
   figure script (generates Fig 1, 2A/B/C/D, 3A–E, Supp 1/6/9, Supp Table 2).
   `v2` reproduces only a subset and `v1` only the tSNE/ML/pheatmap panels.
2. **Almost every `ggsave`/`png()` line is commented out** in all three REAP
   files — figures were exported by hand. Only a few are live
   (`heatmap_baseline_melanoma_R_vs_NR.png`, `baseline_responder_volcano.png`,
   `pt_10_graph.png` in v2; `cytokine_graph.png` in v3).
3. **Stale panel comments:** the `FIG1A`/`FIG1C`/`FIG1D` markers in the REAP code
   predate the final layout. Trust this document's panel letters, not the code
   comments.
4. **Figures 4 & Supp 7** (DC CyTOF/metabolism correlations) have **no source
   in the REAP files** — provide that script separately if it should be mapped.
5. **Figures 5, 6 and Supp 3/4/5/8** come from the **PhIP-seq pipeline**
   (`R/01`–`R/07`), not the REAP scripts.
