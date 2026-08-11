# Build the multi-tab Supplemental Tables workbook (Table-S-master.xlsx)
#
# Assembles the manuscript's key data tables (each already a flat CSV in the
# repo) into one Excel workbook, one tab per table, with a Table-of-Contents
# tab up front. Re-run to regenerate reproducibly after any upstream table
# changes.
#
# Scope: gene-level and summary-level tables that read sensibly in Excel.
# Per-cell tables (barcode-level annotations / QC filters, ~10^4-10^5 rows) are
# intentionally EXCLUDED here -- those belong in the Zenodo data deposit, not a
# supplemental-tables workbook.
#
# Writer: WriteXLS (Perl-backed; committed in renv). Sheet names <= 31 chars.
#
# Output: supplementary-tables/Table-S-master.xlsx  (+ Table-S-contents.csv)

# ---- Prepare Env ----
require(conflicted)
require(tidyverse)
require(glue)
require(cli)
require(WriteXLS)
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, .quiet = TRUE)

setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

out.dir <- "./supplementary-tables"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

read_one <- function(path) {
  stopifnot(file.exists(path))
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  # drop a leading unnamed row-index column if present (e.g. QC table)
  if (names(df)[1] %in% c("", "...1")) df <- df[, -1, drop = FALSE]
  as.data.frame(df)
}

# combine a set of same-family CSVs, tagging provenance with `id`
read_bind <- function(paths, id_from = basename) {
  ids <- id_from(paths)
  purrr::map2(paths, ids, function(p, i) {
    dplyr::mutate(read_one(p), .source = i, .before = 1)
  }) |> dplyr::bind_rows() |> as.data.frame()
}


# ============================================================
# Signature gene lists -> one tidy tab (program, species, set, gene)
# ============================================================
sig_files <- list.files("./references/signatures", pattern = "\\.csv$", full.names = TRUE)
sig_files <- sig_files[!grepl("C3-signature-overlap", sig_files)]  # kept as its own tab
signatures_tidy <- purrr::map_dfr(sig_files, function(p) {
  bn  <- sub("\\.csv$", "", basename(p))
  sp  <- ifelse(grepl("_hs$", bn), "human", "mouse")
  set <- dplyr::case_when(grepl("_refined", bn) ~ "refined",
                          grepl("_original", bn) ~ "original",
                          TRUE ~ "original")
  prog <- bn |> sub("_(hs|mm)$", "", x = _) |> sub("_(refined|original)", "", x = _)
  g <- read_one(p)
  tibble::tibble(program = prog, species = sp, set = set,
                 gene = as.character(g[["gene"]]))
}) |> dplyr::arrange(program, species, set, gene) |> as.data.frame()


# ============================================================
# Curated sheet manifest (order preserved). Each entry -> one tab.
# ============================================================
sheets <- list()
meta   <- list()   # parallel list of TOC metadata

add <- function(sheet, title, ref, df) {
  sheets[[sheet]] <<- df
  meta[[sheet]]   <<- tibble::tibble(Sheet = sheet, Title = title,
                                     N_rows = nrow(df))
}

add("S00_Sample_manifest",
    "Full Flex scRNA-seq sample and library manifest (23 libraries incl. excluded/pilot replicate; mouse, sex, survival, irradiation, hybridization batch, QC; footnotes i-vi)",
    "STAR", read_one("./metadata/sample-table-full-n23.csv"))

add("S00_VisiumHD_manifest",
    "Visium HD spatial sample manifest (4 publication samples: healthy control, 2 primary, 1 post-radiotherapy recurrent; 8/16 um bins)",
    "STAR (Visium HD)", read_one("./metadata/sample-table-visium-hd.csv"))

add("S01_QC_per_sample", "Per-sample QC summary (raw Flex object, 16 samples)",
    "-", read_one("./results/03-qc/03-1-qc-samples-seurat-raw-s16-c92737.csv"))

# ---- Custom RCAS transgene probe sequences (Flex + Visium HD panels) ----
# Oligo sequences for the custom RCAS transgene probes referenced in STAR
# (Key Resources / Oligonucleotides). Built from the committed probe CSVs so the
# table is reproducible from the repo. Flex panel = full designed set; the
# Visium HD panel uses the reduced retained subset.
.probe_cols  <- c("probe_id", "sequence", "aliases", "included", "spliced")
.flex_probes <- readr::read_csv("./references/rcas-tva-probes/10x-flex-custom-probes.csv",
                                col_names = .probe_cols, show_col_types = FALSE, progress = FALSE)
.hd_probes   <- readr::read_csv("./references/rcas-tva-probes/10x-visium-hd-custom-probes.csv",
                                col_names = .probe_cols, show_col_types = FALSE, progress = FALSE)
custom_probes <- .flex_probes |>
  dplyr::transmute(
    Transgene      = sub("_(\\d+|U|nHA)$", "", sub("^CUSTOMPROBE_", "", probe_id)),
    Probe_ID       = probe_id,
    Probe_alias    = vapply(strsplit(aliases, "\\|"),
                            function(x) { tg <- x[grepl("^Tg-", x)]; if (length(tg)) tg[1] else x[1] },
                            character(1)),
    Sequence       = sequence,
    Flex_panel     = "yes",
    VisiumHD_panel = ifelse(probe_id %in% .hd_probes$probe_id, "yes", "no")
  ) |>
  as.data.frame()
add("Custom_probes", "Custom RCAS transgene probe sequences (Flex and Visium HD panels)",
    "STAR (Oligonucleotides)", custom_probes)

add("S02_Markers_Level4", "Cluster markers, all cell types (Level_4, FindAllMarkers)",
    "-", read_one("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Level_4.csv"))

add("S03_Markers_Neoplastic", "Neoplastic subcluster markers (Level_4, FindAllMarkers)",
    "-", read_one("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Neoplastic-Level_4.csv"))

add("S04_Markers_NonNeoplastic", "Non-neoplastic cell-type markers (Level_3, FindAllMarkers)",
    "-", read_one("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-NonNeoplastic-Level_3.csv"))

add("S05_Astro_subtype_markers", "Astrocyte subtype markers (TE / NT / R, FindAllMarkers)",
    "-", read_one("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv"))

add("S06_AstroR_Tumor_vs_Healthy", "Astrocyte R: tumour vs healthy DEG",
    "-", read_one("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstroR-TvH.csv"))

add("S07_NeoplACR_vs_AstroR", "Malignant Neopl-ACR vs reactive Astrocyte R DEG (both directions)",
    "-", read_one("./results/07-6-neopl-acr-vs-astrocyte-r-deg/07-6-FindMarkers-NeoplACR-vs-AstrocyteR.csv"))

add("S08_Signature_gene_lists", "Program signature gene lists (mouse + human, original/refined)",
    "-", signatures_tidy)

add("S09_Signature_refinement", "ACR / NACR signature overlap and refinement membership",
    "-", read_one("./references/signatures/C3-signature-overlap-supplementary-table.csv"))

add("S10_Composition_stats", "Per-condition composition: pairwise Wilcoxon (BH-corrected)",
    "-", read_bind(list.files("./manuscript-figures/figure-2",
                                   pattern = "Figure-S2-composition-stats-.*\\.csv$",
                                   full.names = TRUE)))

add("S11_ExtVal_correlations", "External-cohort signature correlations (Soni/GBmap/Sussman/Suter/Nomura)",
    "-", read_bind(list.files("./results/31-external-validation",
                                   pattern = "\\.csv$", full.names = TRUE)))

add("S12_CellChat_signaling", "CellChat secreted-signalling sender/receiver weights",
    "-", read_bind(list.files("./manuscript-figures/figure-6",
                                   pattern = "CellChat_.*Signaling_(sender|receiver)\\.csv$",
                                   full.names = TRUE)))

add("S13_BulkIR_DESeq2", "Human astrocyte irradiation bulk RNA-seq: DESeq2 IR vs non-IR",
    "-", read_one("./results/07-5-astrocyte-ir-bulkseq/DESeq2-IR-vs-nonIR-DMSO.csv"))

add("S14_BulkIR_GSEA", "Astrocyte R program GSEA in irradiation-ranked bulk RNA-seq",
    "-", read_one("./results/07-5-astrocyte-ir-bulkseq/GSEA-AstrocyteR-program-in-IR.csv"))

# ---- Astrocyte R recurrence DEG (belongs with the S06/S07 DEG lists) ----
add("S17_AstroR_Recur_vs_Prim", "Astrocyte R: post-radiotherapy recurrent vs primary DEG",
    "Fig 3I", read_one("./manuscript-figures/figure-3/Figure-3-Volcano_AstroR_TRvsTP-sourcedata.csv"))

# ---- Visium HD composition-niche cell-type matrix (Fig 7F) ----
# Per-niche mean RCTD weight for every cell type, plus the summed malignant (Neopl-*) fraction
# and the tumor/interface/normal class that follows from it. Documents the 12 niche definitions
# and the tumor-vs-normal call. NB: re-run 16-3 first so the source CSV carries the final niche
# names (name niches by dominant component, not the rare Neopl-ACR/ECM enrichment).
niche_wide <- read_one("./manuscript-figures/figure-7/Figure-7-niches-celltype-heatmap-sourcedata.csv") |>
  dplyr::select(niche_name, cell_type, mw) |>
  tidyr::pivot_wider(names_from = cell_type, values_from = mw)
.neopl <- grep("^Neopl", names(niche_wide), value = TRUE)
niche_wide <- niche_wide |>
  dplyr::mutate(Neopl_fraction = round(rowSums(dplyr::across(dplyr::all_of(.neopl))), 3),
                Class = dplyr::case_when(Neopl_fraction > 0.5  ~ "Tumor",
                                         Neopl_fraction > 0.15 ~ "Interface",
                                         TRUE ~ "Normal"),
                .after = niche_name) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))) |>
  dplyr::arrange(dplyr::desc(Neopl_fraction)) |>
  as.data.frame()
add("S18_Niche_composition",
    "Visium HD composition-niche cell-type matrix (mean RCTD weight per niche; malignant fraction + tumor/normal class; Fig 7F)",
    "-", niche_wide)

# NOTE: enrichment / activity / stats results (AUCell, fgsea, CellChat rankNet, GFAP correlation)
# are FIGURE SOURCE DATA shipped per panel (`*-sourcedata.csv`), NOT manuscript supplemental
# tables — intentionally NOT added here. Manuscript supp tables = gene lists (markers, DEGs,
# signatures), oligos, sample/QC, plus the niche-composition matrix (S18, added by request).
# If reviewers need the other analysis results collated, build a separate reviewers-only workbook.


# ============================================================
# Single Table S1: the whole workbook IS "Table S1"; every sheet is a tab within
# it, NOT a separately numbered Table. Strip the "S0x_" ordering prefixes so tab
# labels read as descriptive sheets and are not mistaken for Tables S1-S18.
# ============================================================
strip_num <- function(x) sub("^S\\d+_", "", x)
names(sheets) <- strip_num(names(sheets))
meta <- lapply(meta, function(m) { m$Sheet <- strip_num(m$Sheet); m })

# ============================================================
# Table of contents (first tab) + write workbook
# ============================================================
toc <- dplyr::bind_rows(meta) |>
  dplyr::mutate(Tab = dplyr::row_number(), .before = 1) |>
  as.data.frame()

all_sheets <- c(list(TOC = toc), sheets)

readr::write_csv(toc, file.path(out.dir, "Table-S-contents.csv"))

xlsx_path <- file.path(out.dir, "Table-S-master.xlsx")
WriteXLS::WriteXLS(
  x            = all_sheets,
  ExcelFileName = xlsx_path,
  SheetNames   = names(all_sheets),
  row.names    = FALSE,
  AdjWidth     = TRUE,
  BoldHeaderRow = TRUE,
  FreezeRow    = 1,
  na           = ""
)

cli::cli_alert_success("Wrote {xlsx_path} ({length(all_sheets)} tabs incl. TOC)")
print(toc)
