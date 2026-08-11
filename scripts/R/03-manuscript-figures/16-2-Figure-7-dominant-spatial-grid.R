# 16-2 — Figure 7: dominant (argmax) RCTD cell-type spatial grid (16um, 4 samples)
# ============================================================================
# Each 16um bin -> its single dominant (argmax) RCTD cell type (14-type schema), rendered
# in the canonical spatial style (true square bins + Tg/tissue/exclusion outlines) with the
# predefined .color_pal Level_4ACM identity colours. Neuron is paled to near-background so the
# tumour/astrocyte/immune structure is legible. Reads the per-bin RCTD weight CSVs from
# results/08-visium-hd-derived (env VISIUM_HD_DERIVED_DIR); RCTD deconvolution is
# upstream (not re-run here).
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-0-source-functions.R"))    # .color_pal
source(file.path(fig_dir, "_fig1d-spatial-style.R"))
source(file.path(fig_dir, "_mask-outlines.R"))                       # mask_outlines, spatial_bin_grid
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
DERIVED   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))

# ---- 14-type -> predefined .color_pal identity colours ----------------------
pal4 <- .color_pal[["Level_4ACM"]]
rep_map <- c("Neopl-ACR"="Neopl-ACR","Neopl-ECM"="Neopl-ECM","Neopl-OPC"="Neopl-OPC",
             "Neopl-COP"="Neopl-COP","Neopl-Bulk"="Neopl-Bulk",
             "Astrocyte R"="Astrocyte R","Astrocyte TE-NT"="Astrocyte TE/NT","OPC-COP-OLG"="OPC/COP/OLG",
             "Neuron"="Neural","Ependymal"="Ependymal","Choroid"="Choroid","Myeloid"="Myeloid",
             "Immune-other"="Immune-other","Vascular"="Endothelial")
CT <- names(rep_map)
CT_PAL <- setNames(unname(pal4[rep_map]), CT); CT_PAL[is.na(CT_PAL)] <- "grey40"

SAMP_ORD <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")
.cond <- function(s) ifelse(grepl("^HB", s), "Healthy", ifelse(grepl("^TP", s), "Primary", "PostRT"))

# ---- load all per-bin RCTD weights + argmax dominant type -------------------
raw <- lapply(SAMP_ORD, function(s)
  read_csv(file.path(DERIVED, sprintf("rctd-weights-%s-16um.csv.gz", s)), show_col_types = FALSE)) |> bind_rows()
raw$dominant  <- factor(CT[max.col(as.matrix(raw[, CT]), ties.method = "first")], CT)
raw$sample_id <- factor(raw$sample_id, SAMP_ORD)

outlines <- bind_rows(lapply(SAMP_ORD, mask_outlines)); outlines$sample_id <- factor(outlines$sample_id, SAMP_ORD)

# FIGURE DEFAULT = Neuron paled: neuronal parenchyma dominates the argmax
# and swamps the tumour/astrocyte/immune structure, so it is set near-background pale blue.
CT_PAL_fig <- CT_PAL; CT_PAL_fig[["Neuron"]] <- "#DCE6F4"
p_grid <- spatial_bin_grid(raw, "dominant", SAMP_ORD, outlines,
             title = "Dominant RCTD cell type per bin - Neuron paled (blue=Tg tumor, grey=tissue, dashed=excluded)",
             fill_scale = scale_fill_manual(values = CT_PAL_fig, name = "dominant cell type", drop = FALSE))
ggsave(file.path(plot_dir, "Figure-7-dominant-spatial-grid.pdf"), p_grid, width = 15, height = 4.6)
write_csv(raw[, c("sample_id","barcode","x","y","dominant")],
          file.path(plot_dir, "Figure-7-dominant-spatial-grid-sourcedata.csv.gz"))
cli::cli_alert_success("16-2 dominant RCTD spatial grid -> {plot_dir}")
