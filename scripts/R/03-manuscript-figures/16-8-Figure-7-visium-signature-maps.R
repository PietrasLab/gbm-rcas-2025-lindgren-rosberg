# 16-8 — Figure 7 supplement: Visium HD spatial maps of the GBM-state signatures
# ============================================================================
# The our-model companion to MF5 (Heiland UKF260 maps): the same GBM-state programs
# scored per 16um bin on the 4 publication Visium HD samples and shown as spatial
# maps, STYLED EXACTLY AS MF4 (Fig 3C / script 12-3): row = signature (rotated
# left-hand label), cols = HB_hm01 / TP_1083 / TP_987 / TR_03d_990, fig1d fill +
# Tg/tissue outlines. INDIVIDUAL per-signature PDFs (no composite). Signature-based
# (independent of RCTD); illustrative (small spatial cohort). Signature sources =
# Fig 4 / MF5 (mouse Richards/Neftel RDS + refined _mm CSVs).
#
# INPUTS (as 12-3): 06-2-visium-hd-tg-mask.R staged outlines (results/06-2-tg-masks);
#   per-sample 16um Seurat objects (Zenodo; env VISIUM_HD_OBJ_ROOT); bin tissue mask
#   (tracked: references/visium-hd-masks/masks-workspace; env VISIUM_HD_MASK_DIR).
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(ggplot2); library(patchwork)
  library(RANN); library(grid); library(msigdbr); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-2-mask-utils.R"))     # apply_tissue_mask
source(file.path(fig_dir, "_fig1d-spatial-style.R"))            # fig1d_fill
source(file.path(fig_dir, "_mask-outlines.R"))                  # mask_outlines, spatial_bin_grid
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
OBJ_ROOT  <- Sys.getenv("VISIUM_HD_OBJ_ROOT",
  unset = file.path(proj_root, "data/processed/seurat/visium-hd"))
MASK_DIR  <- Sys.getenv("VISIUM_HD_MASK_DIR",
  unset = file.path(proj_root, "references/visium-hd-masks/masks-workspace"))
sig_dir   <- file.path(proj_root, "references/signatures")
BIN <- 16L
SAMP_ORD  <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")

# ---- MF5 GBM-state signature set (mouse), same sources as Fig 4 / MF5 --------
gs_mm <- readRDS(file.path(proj_root, "references/genesets/Richards_NatCancer_2021_GeneSets_Mouse.rds"))
sig.list <- list(
  "Neopl-ACR"                = read_csv(file.path(sig_dir, "neopl_ACR_refined_mm.csv"), show_col_types = FALSE)$gene,
  "Neftel MES1"              = gs_mm[["Neftel_Cell_2019_MES1"]],
  "Richards Injury-Response" = gs_mm[["InHouse_BulkRNAseq_2019_InjuryResponseGSC"]],
  "Richards Developmental"   = gs_mm[["InHouse_BulkRNAseq_2019_DevelopmentalGSC"]],
  "HALLMARK Hypoxia"         = msigdbr(species = "Mus musculus", category = "H") |>
                                 dplyr::filter(gs_name == "HALLMARK_HYPOXIA") |> dplyr::pull(gene_symbol) |> unique(),
  "Astrocyte R"              = read_csv(file.path(sig_dir, "astrocyte_R_refined_mm.csv"), show_col_types = FALSE)$gene)
SIGS <- names(sig.list)
cli::cli_alert_info("signatures: {paste(sprintf('%s(%d)', SIGS, lengths(sig.list)), collapse=', ')}")

# ---- score each sample's 16um bins with AddModuleScore (as 12-3) ------------
df <- lapply(SAMP_ORD, function(sid) {
  cli::cli_alert("score {sid}")
  o <- readRDS(file.path(OBJ_ROOT, sprintf("seurat_visium-hd_%s_%dum_v1.1.rds", sid, BIN)))
  assay <- grep("Spatial", Assays(o), value = TRUE)[1]; DefaultAssay(o) <- assay
  o <- NormalizeData(o, verbose = FALSE); o <- apply_tissue_mask(o, sid, MASK_DIR, bin_size = BIN)
  feats <- lapply(sig.list, function(g) intersect(g, rownames(o)))
  o <- AddModuleScore(o, features = feats, name = "M_", assay = assay, seed = 42)
  sc <- o@meta.data[, paste0("M_", seq_along(SIGS)), drop = FALSE]; colnames(sc) <- SIGS
  co <- GetTissueCoordinates(o); cc <- intersect(c("x","y","imagerow","imagecol"), colnames(co))[1:2]
  d <- data.frame(sample_id = sid, barcode = colnames(o), x = co[[cc[1]]], y = co[[cc[2]]], sc, check.names = FALSE)
  rm(o); gc(FALSE); d
}) |> bind_rows()
df$sample_id <- factor(df$sample_id, SAMP_ORD)
write_csv(df, file.path(plot_dir, "Figure-7-visium-signature-maps-sourcedata.csv.gz"))
outlines <- bind_rows(lapply(SAMP_ORD, mask_outlines)); outlines$sample_id <- factor(outlines$sample_id, SAMP_ORD)

# ---- MF4 style: one row per signature (rotated label); per-signature 1-99% fill ----
row_lab <- function(g) {
  lim <- quantile(df[[g]], c(0.01, 0.99), na.rm = TRUE)
  d <- df; d[[g]] <- pmin(pmax(d[[g]], lim[1]), lim[2])
  p   <- spatial_bin_grid(d, g, SAMP_ORD, outlines, title = NULL,
                          fill_scale = fig1d_fill(limits = lim, name = g))
  lab <- patchwork::wrap_elements(grid::textGrob(g, rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 12)))
  patchwork::wrap_plots(lab, patchwork::wrap_elements(full = p), widths = c(1, 32), nrow = 1)
}
# INDIVIDUAL per-signature PDFs (no composite)
for (g in SIGS) ggsave(
  file.path(plot_dir, sprintf("Figure-7-visium-signature-map-%s.pdf", gsub("[^A-Za-z0-9]+", "_", g))),
  row_lab(g), width = 15.5, height = 4.6)
cli::cli_alert_success("16-8 Visium HD signature maps (individual, MF4-style, {length(SIGS)} sigs) -> {plot_dir}")
