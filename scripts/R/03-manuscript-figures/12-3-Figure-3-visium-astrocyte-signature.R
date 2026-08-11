# 12-3 — Figure 3C (updated): Visium HD astrocyte-subtype signature maps (4 samples)
# ============================================================================
# Updated Figure 3C. The original 3C (12-1-Figure-3.R) scored the astrocyte-subtype gene
# SIGNATURES (top-10 markers per subtype from 07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv)
# with AddModuleScore on the 16um Visium HD bins, for the 2 pilot samples (H + TP_1083). This
# rebuild extends it to the 4 publication samples and flips the layout to rows = subtype:
#   rows = Astrocyte TE, Astrocyte NT, Astrocyte R  (Pan-Astrocyte dropped to limit the panel)
#   cols = HB_hm01 (healthy), TP_1083, TP_987 (primary), TR_03d_990 (post-RT recurrent)
# Signature-based (method-independent of RCTD vs TransferAnchors), 16um (Spatial.016um), gradient
# capped at 0.75. Canonical spatial style (true square bins, fig1d fill, Tg/tissue outlines) with a
# rotated left-hand subtype label per row.
#
# INPUTS (rerun order):
#   - Tg overlay + QuPath annotation masks: run 06-2-visium-hd-tg-mask.R first (stages the outlines
#     into results/06-2-tg-masks/, gitignored) - same dependency as 10-6.
#   - Per-sample 16um Seurat objects (Zenodo): data/processed/seurat/visium-hd (env
#     VISIUM_HD_OBJ_ROOT); bin tissue mask (tracked): references/visium-hd-masks/masks-workspace (env VISIUM_HD_MASK_DIR).
#   - Astrocyte-subtype markers: results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv.
# OUTPUT: manuscript-figures/figure-3/ (combined 3-row PDF + per-subtype PDFs + sourcedata)
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(RANN); library(grid); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-2-mask-utils.R"))          # apply_tissue_mask
source(file.path(fig_dir, "_fig1d-spatial-style.R"))                 # fig1d_fill, theme_fig1d
source(file.path(fig_dir, "_mask-outlines.R"))                       # mask_outlines, spatial_bin_grid (reads results/06-2-tg-masks)
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-3"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

OBJ_ROOT  <- Sys.getenv("VISIUM_HD_OBJ_ROOT",
  unset = file.path(proj_root, "data/processed/seurat/visium-hd"))
MASK_DIR  <- Sys.getenv("VISIUM_HD_MASK_DIR",
  unset = file.path(proj_root, "references/visium-hd-masks/masks-workspace"))
BIN <- 16L; TOPN <- 10L; CAP <- 0.75
SAMP_ORD  <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")
ASTRO     <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R")    # rows (Pan dropped)

# ---- astrocyte-subtype signatures: top-10 markers per subtype (as 12-1) ------
deg <- read_csv(file.path(proj_root, "results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv"), show_col_types = FALSE)
sig <- deg |> filter(avg_log2FC > 0) |> group_by(cluster) |> arrange(p_val_adj, .by_group = TRUE) |>
  slice_head(n = TOPN) |> summarise(genes = list(gene), .groups = "drop")
sig.list <- setNames(sig$genes, sig$cluster)[ASTRO]
cli::cli_alert_info("signatures: {paste(sprintf('%s(%d)', names(sig.list), lengths(sig.list)), collapse=', ')}")

# ---- score each sample's 16um bins with AddModuleScore -----------------------
df <- lapply(SAMP_ORD, function(sid) {
  cli::cli_alert("score {sid}")
  o <- readRDS(file.path(OBJ_ROOT, sprintf("seurat_visium-hd_%s_%dum_v1.1.rds", sid, BIN)))
  assay <- grep("Spatial", Assays(o), value = TRUE)[1]; DefaultAssay(o) <- assay
  o <- NormalizeData(o, verbose = FALSE); o <- apply_tissue_mask(o, sid, MASK_DIR, bin_size = BIN)
  feats <- lapply(sig.list, function(g) intersect(g, rownames(o)))
  o <- AddModuleScore(o, features = feats, name = "AC_", assay = assay, seed = 42)
  sc <- o@meta.data[, paste0("AC_", seq_along(ASTRO)), drop = FALSE]; colnames(sc) <- ASTRO
  co <- GetTissueCoordinates(o); cc <- intersect(c("x","y","imagerow","imagecol"), colnames(co))[1:2]
  d <- data.frame(sample_id = sid, barcode = colnames(o), x = co[[cc[1]]], y = co[[cc[2]]], sc, check.names = FALSE)
  rm(o); gc(FALSE); d
}) |> bind_rows()
df$sample_id <- factor(df$sample_id, SAMP_ORD)
write_csv(df, file.path(plot_dir, "Figure-3-ModuleScores-VisiumHD-AstroTENTR-4sample-sourcedata.csv.gz"))
outlines <- bind_rows(lapply(SAMP_ORD, mask_outlines)); outlines$sample_id <- factor(outlines$sample_id, SAMP_ORD)

# one canonical row per astrocyte subtype (fig1d fill, fixed 0..CAP so rows are comparable).
# each row carries a rotated left-hand label (the subtype) instead of a top title.
row_lab <- function(g) {
  p   <- spatial_bin_grid(df, g, SAMP_ORD, outlines, title = NULL,
                          fill_scale = fig1d_fill(limits = c(0, CAP), name = g))
  lab <- patchwork::wrap_elements(grid::textGrob(g, rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 13)))
  patchwork::wrap_plots(lab, patchwork::wrap_elements(full = p), widths = c(1, 32), nrow = 1)
}
rows <- lapply(ASTRO, row_lab)
combined <- wrap_plots(rows, ncol = 1)
ggsave(file.path(plot_dir, "Figure-3-ModuleScores-VisiumHD-AstroTENTR-4sample.pdf"), combined, width = 15.5, height = 13.5)
for (i in seq_along(ASTRO)) ggsave(
  file.path(plot_dir, sprintf("Figure-3-ModuleScores-VisiumHD-Astro-%s-4sample.pdf", gsub("[^A-Za-z0-9]+","_",ASTRO[i]))),
  rows[[i]], width = 15.5, height = 4.6)
cli::cli_alert_success("12-3 astrocyte signature panel (AddModuleScore 16um, cap={CAP}) -> {plot_dir}")
