# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(glue)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
  library(dplyr)
})
conflicts_prefer(dplyr::filter, .quiet = TRUE)

# Canonical manuscript palette
source("./scripts/R/00-0-source-functions.R")

# 10-3 — Supplementary Figure 1 (Figure 1 supplement): per-sample split UMAP
# (revision). Shows each sample spans the shared embedding,
# i.e. clustering is not sample- or downsampling-driven.
# Output in the matching main-figure folder (figure-1/) as Figure-S1-* per repo convention.

obj_path <- "data/processed/seurat/seurat_flex_filtered_v1.0.rds"
plot.dir <- "./manuscript-figures/figure-1"
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

sample_levels <- c("BRAIN_01","BRAIN_04a","BRAIN_04b","BRAIN_05","BRAIN_06",
                   "TP_01b","TP_04","TP_05","TP_06","TP_08",
                   "TR_01","TR_02","TR_03","TR_04","TR_06","TR_07")

cli_alert_info("Loading object {obj_path}")
obj <- readRDS(obj_path)
obj$sample_id   <- factor(obj$sample_id, levels = sample_levels)
obj$sample_type <- factor(obj$sample_type, levels = c("Healthy","Primary","Recurrent"))

sample_colors    <- .color_pal[["sample_id"]][sample_levels]   # canonical
condition_colors <- .color_pal[["sample_type"]]

theme_umap <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          axis.text = element_blank(), axis.ticks.length = unit(0, "pt"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = base_size))
}

# ---- one highlighted mini-UMAP per sample (grey background) --------------------
make_sample_panel <- function(sid) {
  col <- sample_colors[[sid]]
  # small vivid points, with the point layer rasterised
  # via ggrastr so the 16-panel supp PDF stays small (~MB, not 57 MB vector).
  p <- DimPlot(obj,
          cells.highlight = WhichCells(obj, expression = sample_id == sid),
          cols.highlight  = col, cols = "grey85",
          pt.size = 0.08, sizes.highlight = 0.12, raster = FALSE, label = FALSE) +
    labs(title = sid, x = NULL, y = NULL) +
    theme_umap(8) +
    theme(legend.position = "none", plot.title = element_text(color = col, face = "bold"))
  rasterise(p, dpi = 300)
}

# robust row label (vertical text strip) so condition grouping is explicit
row_label <- function(txt, col) {
  ggplot() + annotate("text", x = 0, y = 0, label = txt, color = col,
                      fontface = "bold", size = 4, angle = 90) +
    theme_void()
}
make_row <- function(idx, txt) {
  panels <- lapply(sample_levels[idx], make_sample_panel)
  wrap_plots(c(list(row_label(txt, condition_colors[[txt]])), panels), nrow = 1,
             widths = c(0.12, rep(1, length(panels))))
}

cli_alert_info("Building 16 per-sample panels...")
p_split <- make_row(1:5,   "Healthy") / make_row(6:10, "Primary") / make_row(11:16, "Recurrent")

ggsave(file.path(plot.dir, "Figure-S1-umap-per-sample.pdf"), p_split, width = 17, height = 11)
cli_alert_success("Supp Fig 1 (per-sample UMAP) -> {plot.dir}/Figure-S1-umap-per-sample.pdf")

# Source data: shared embedding + per-sample / condition grouping (panels are highlight
# overlays of the same UMAP), one row per cell.
readr::write_csv(
  cbind(
    Embeddings(obj, "umap_AllCells")[, 1:2],
    obj@meta.data[, c("sample_id", "sample_type"), drop = FALSE]
  ),
  file.path(plot.dir, "Figure-S1-umap-per-sample-sourcedata.csv")
)
