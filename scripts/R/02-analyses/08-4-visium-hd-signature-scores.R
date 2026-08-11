# 08-4 — Visium HD per-bin signature scores (harmonized panel)
# ============================================================================
# Scores the harmonized signature panel (our 6 neoplastic programs + 14 non-
# neoplastic + Neftel 6 / Richards 2 / Nomura 10 reference states, each tagged
# with source + compartment; mouse lists) per 16um bin with AddModuleScore
# (seed = 42), plus HALLMARK_HYPOXIA. Attaches Tg-tumor membership and the
# 08-3 neoplastic fraction. Consumed by the Figure 7 correlation script (16-4).
#
# INPUTS:
#   - Per-sample 16um Seurat objects (Zenodo): data/processed/seurat/visium-hd
#     (env VISIUM_HD_OBJ_ROOT)
#   - Tissue masks (tracked): references/visium-hd-masks (env VISIUM_HD_MASK_DIR)
#   - Signature panel: scripts/R/03-manuscript-figures/_signature-panel.R
#     (reads references/signatures + references/genesets)
#   - 08-3 RCTD weights (run 08-3 first; optional, adds neopl_frac)
# OUTPUT (results/08-visium-hd-derived/):
#   signature-scores-<sample>-16um.csv.gz + signature-panel-meta.csv
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(tibble); library(cli)
})
set.seed(42)
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-2-mask-utils.R"))
source(file.path(fig_dir, "_mask-outlines.R"))        # tg_membership
source(file.path(fig_dir, "_signature-panel.R"))      # build_panel_full
out_dir   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
OBJ_ROOT  <- Sys.getenv("VISIUM_HD_OBJ_ROOT",
  unset = file.path(proj_root, "data/processed/seurat/visium-hd"))
mask_dir  <- Sys.getenv("VISIUM_HD_MASK_DIR",
  unset = file.path(proj_root, "references/visium-hd-masks/masks-workspace"))
BIN_SIZE  <- 16L
SAMPLES   <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")

# ---- harmonized panel — BOTH original + REFINED Neopl-ACR/Astrocyte-R targets,
# so the HD corrplots match the Heiland analyses AND de-confound the
# ACR<->reactive-astro overlap (shared genes removed from the refined lists).
panel <- build_panel_full("mouse", acr_targets = "both")
# append HALLMARK_HYPOXIA (hypoxia co-enrichment is called out in the Fig 7 text
# but is not part of build_panel_full); source='Hallmark', compartment='neoplastic'.
hyp_genes <- msigdbr::msigdbr(species = "Mus musculus", category = "H") |>
  dplyr::filter(gs_name == "HALLMARK_HYPOXIA") |> dplyr::pull(gene_symbol) |> unique()
panel$sigs[["HALLMARK_HYPOXIA"]] <- hyp_genes
panel$meta <- dplyr::bind_rows(panel$meta,
  tibble::tibble(label = "HALLMARK_HYPOXIA", source = "Hallmark", compartment = "neoplastic"))
write_csv(panel$meta, file.path(out_dir, "signature-panel-meta.csv"))
cli::cli_alert_info("panel: {nrow(panel$meta)} sigs ({sum(panel$meta$compartment=='neoplastic')} neoplastic); ACR/Astro-R = BOTH original+refined")

# ---- score each sample (AddModuleScore, seed 42) ---------------------------
for (sid in SAMPLES) {
  cli::cli_h2("{sid}")
  o <- readRDS(file.path(OBJ_ROOT, sprintf("seurat_visium-hd_%s_%dum_v1.1.rds", sid, BIN_SIZE)))
  assay <- grep("Spatial", Assays(o), value = TRUE)[1]; DefaultAssay(o) <- assay
  o <- NormalizeData(o, verbose = FALSE); o <- apply_tissue_mask(o, sid, mask_dir, bin_size = BIN_SIZE)
  feats <- lapply(panel$sigs, function(g) intersect(g, rownames(o)))
  ok <- vapply(feats, length, integer(1)) >= 5; feats <- feats[ok]
  o <- AddModuleScore(o, features = feats, name = "SIG_", assay = assay, seed = 42)
  sc <- o@meta.data[, paste0("SIG_", seq_along(feats)), drop = FALSE]; colnames(sc) <- names(feats)
  co <- GetTissueCoordinates(o); cc <- intersect(c("x","y","imagerow","imagecol"), colnames(co))[1:2]
  df <- data.frame(sample_id = sid, barcode = colnames(o), x = co[[cc[1]]], y = co[[cc[2]]], sc, check.names = FALSE)
  df$in_tg <- tg_membership(sid, df$x, df$y)
  w <- file.path(out_dir, sprintf("rctd-weights-%s-%dum.csv.gz", sid, BIN_SIZE))
  if (file.exists(w)) { wd <- read_csv(w, show_col_types = FALSE); df <- left_join(df, wd[, c("barcode","neopl_frac")], by = "barcode") }
  write_csv(df, file.path(out_dir, sprintf("signature-scores-%s-%dum.csv.gz", sid, BIN_SIZE)))
  cli::cli_alert_success("{sid}: {nrow(df)} bins x {length(feats)} sigs | Tg-tumor = {sum(df$in_tg)}")
  rm(o); gc(FALSE)
}
cli::cli_alert_success("08-4 panel scores -> signature-scores-*.csv.gz (+ signature-panel-meta.csv)")
