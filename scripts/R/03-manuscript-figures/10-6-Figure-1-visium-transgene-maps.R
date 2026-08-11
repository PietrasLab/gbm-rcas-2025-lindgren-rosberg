# 10-6 — Figure 1 (panel TBD; likely 1E/1F): Visium HD Olig2 + RCAS transgene maps
# ============================================================================
# Spatial maps of Olig2 + the RCAS transgenes (Tg-hPDGFB-nHA driver, Tg-RFP reporter)
# on the 4 publication Visium HD samples (HB_hm01 healthy control, TP_1083, TP_987,
# TR_03d_990), in the canonical spatial style: true 16um square bins, fig1d continuous
# fill, and the Tg / tissue / exclusion mask outlines. HB_hm01 = the transgene-negative
# control (no Tg overlay). Supersedes the transgene portion of the pilot 10-2-Figure-1-
# visium.R (which covered only the 2 pilots).
#
# INPUTS (rerun order):
#   - Tg overlay + QuPath annotation masks: run 06-2-visium-hd-tg-mask.R first
#     (stages into results/06-2-tg-masks/, gitignored).
#   - Per-sample 16um Seurat objects (Zenodo): data/processed/seurat/visium-hd
#     (env VISIUM_HD_OBJ_ROOT). Bin tissue mask + capture parquet (tracked):
#     references/visium-hd-masks/masks-workspace (env VISIUM_HD_MASK_DIR).
# OUTPUT: manuscript-figures/figure-1/ (combined 3-row PDF + per-feature PDFs + sourcedata)
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(RANN); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-2-mask-utils.R"))          # apply_tissue_mask
source(file.path(fig_dir, "_fig1d-spatial-style.R"))                 # fig1d_fill, theme_fig1d
source(file.path(fig_dir, "_mask-outlines.R"))                       # mask_outlines, spatial_bin_grid (reads results/06-2-tg-masks)
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-1"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

OBJ_ROOT  <- Sys.getenv("VISIUM_HD_OBJ_ROOT",
  unset = file.path(proj_root, "data/processed/seurat/visium-hd"))
MASK_DIR  <- Sys.getenv("VISIUM_HD_MASK_DIR",
  unset = file.path(proj_root, "references/visium-hd-masks/masks-workspace"))
BIN_SIZE  <- 16L
SAMP_ORD  <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")
GENES     <- c("Olig2", "Tg-hPDGFB-nHA", "Tg-RFP")
GENE_LAB  <- c("Olig2" = "Olig2 (oligo-lineage)", "Tg-hPDGFB-nHA" = "hPDGFB (RCAS driver)", "Tg-RFP" = "RFP (RCAS reporter)")

# ---- fetch the 3 features per bin (normalized 'data' layer) -----------------
df <- lapply(SAMP_ORD, function(sid) {
  cli::cli_alert("fetch {sid}")
  o <- readRDS(file.path(OBJ_ROOT, sprintf("seurat_visium-hd_%s_%dum_v1.1.rds", sid, BIN_SIZE)))
  assay <- grep("Spatial", Assays(o), value = TRUE)[1]; DefaultAssay(o) <- assay
  o <- NormalizeData(o, verbose = FALSE); o <- apply_tissue_mask(o, sid, MASK_DIR, bin_size = BIN_SIZE)
  co <- GetTissueCoordinates(o); cc <- intersect(c("x","y","imagerow","imagecol"), colnames(co))[1:2]
  present <- intersect(GENES, rownames(o))
  ex <- if (length(present)) FetchData(o, vars = present, layer = "data") else NULL
  d <- tibble(sample_id = sid, barcode = colnames(o), x = co[[cc[1]]], y = co[[cc[2]]])
  for (g in GENES) d[[g]] <- if (!is.null(ex) && g %in% colnames(ex)) ex[[g]] else NA_real_
  rm(o); gc(FALSE); d
}) |> bind_rows()
df$sample_id <- factor(df$sample_id, SAMP_ORD)
write_csv(df, file.path(plot_dir, "Figure-1-visium-transgene-olig2-sourcedata.csv.gz"))

outlines <- bind_rows(lapply(SAMP_ORD, mask_outlines)); outlines$sample_id <- factor(outlines$sample_id, SAMP_ORD)

# ---- one canonical row per feature (fig1d fill, per-gene 1-99% clamp) --------
rows <- lapply(GENES, function(g) {
  lim <- quantile(df[[g]], c(0.01, 0.99), na.rm = TRUE)
  if (!is.finite(lim[2]) || lim[2] <= lim[1]) lim <- c(0, max(df[[g]], na.rm = TRUE))
  spatial_bin_grid(df, g, SAMP_ORD, outlines, title = GENE_LAB[[g]],
                   fill_scale = fig1d_fill(limits = as.numeric(lim), name = g))
})
combined <- wrap_plots(rows, ncol = 1)
ggsave(file.path(plot_dir, "Figure-1-visium-transgene-olig2-maps.pdf"), combined, width = 15, height = 13.5)
for (i in seq_along(GENES)) {
  if (GENES[i] == "Tg-RFP") next   # Tg-RFP is the Supp Fig S1A panel (saved below), not a main Figure 1 map
  ggsave(file.path(plot_dir, sprintf("Figure-1-visium-map-%s.pdf", gsub("[^A-Za-z0-9]+","_",GENES[i]))), rows[[i]], width = 15, height = 4.6)
}
# The RFP reporter row is the Supplementary Figure S1A panel (4-sample map, replaces the old
# 2-sample CUSTOMPROBE-RFP S1A). Tg-RFP lives ONLY here, not as a duplicate main Figure 1 map.
ggsave(file.path(plot_dir, "Figure-S1-visium-map-Tg_RFP.pdf"), rows[[match("Tg-RFP", GENES)]], width = 15, height = 4.6)
cli::cli_alert_success("10-6 transgene + Olig2 maps -> {plot_dir}")
