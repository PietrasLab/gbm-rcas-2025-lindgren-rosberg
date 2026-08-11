# 11-5 — Supp Fig 2 (Figure 2 supplement): per-cell-type QC on uncapped UMI
# ============================================================================
# Shows the cell types are not QC artefacts and the conclusions are not shaped
# by UMI-downsampling.
# Level_3 (broad cell-type) companion of the Level_4 per-cluster QC panel
# (13-2 Panel A, Figure S4). Same computation, coarser grouping (Level_3).
#
# Panel:
#   A  Per-cell-type QC (nCount/nFeature/percent.mt) on RAW, UNCAPPED UMI across
# all Level_3 broad cell types.
# NB reads the RAW (non-downsampled) object only for uncapped counts, matched to
#    the published Level_3 labels — no re-clustering needed (that is 11-3).
# ============================================================================

suppressMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(tidyr)
  library(readr); library(cli)
})
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
source("./scripts/R/00-0-source-functions.R")   # canonical .color_pal

plot.dir <- "./manuscript-figures/figure-2"
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

# ---- objects ---------------------------------------------------------------
cli::cli_alert_info("Loading UMI-downsampled object (for Level_3 labels)")
obj <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")

cli::cli_alert_info("Loading raw (non-downsampled) object for uncapped QC")
raw <- readRDS("./data/processed/seurat/seurat_flex_raw_v1.0.rds")
if (length(SeuratObject::Layers(raw, search = "counts")) > 1)
  raw <- SeuratObject::JoinLayers(raw)

# ---- match raw cells to published Level_3 labels ---------------------------
# published : {sample}_{barcode}-1     raw : {barcode}-1_{mergeIndex}
pub_core <- mapply(function(cn, s) sub(paste0("^", s, "_"), "", cn),
                   colnames(obj), as.character(obj$orig.ident), USE.NAMES = FALSE)
ann      <- setNames(as.character(obj$Level_3), paste(as.character(obj$orig.ident), pub_core))
raw_key  <- paste(as.character(raw$orig.ident), sub("_\\d+$", "", colnames(raw)))
raw <- raw[, raw_key %in% names(ann)]
raw$Level_3 <- unname(ann[raw_key[raw_key %in% names(ann)]])
raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^mt-")  # mito genes present pre-filter
cli::cli_alert_info("Uncapped-QC cells (raw ∩ published): {ncol(raw)}")

# ---- A — per-cell-type QC on uncapped UMI (Level_3) -------------------------
lvl_order <- intersect(names(.color_pal[["Level_3"]]), unique(raw$Level_3))
qc_all <- raw@meta.data |>
  select(Level_3, nCount_RNA, nFeature_RNA, percent.mt) |>
  mutate(Level_3 = factor(Level_3, levels = lvl_order))
qc_long <- tidyr::pivot_longer(qc_all, c(nCount_RNA, nFeature_RNA, percent.mt),
                               names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric,
                         levels = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                         labels = c("UMI count (uncapped)", "Genes detected", "Mitochondrial %")))
p_qc <- ggplot(qc_long, aes(Level_3, value, fill = Level_3)) +
  geom_violin(scale = "width", alpha = 0.85, show.legend = FALSE) +
  geom_boxplot(width = 0.07, outlier.shape = NA, fill = "white", alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = .color_pal[["Level_3"]], guide = "none") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(title = "Per-cell-type QC across all Level_3 states (RAW, uncapped UMI)",
       subtitle = "Non-UMI-downsampled counts; broad cell types ordered by compartment (neural → neoplastic)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot.dir, "Figure-S2-QC-uncapped-Level3-violins.pdf"), p_qc, width = 10, height = 10)
write_csv(qc_all, file.path(plot.dir, "Figure-S2-QC-uncapped-Level3-violins-sourcedata.csv"))

cli::cli_alert_success("Supp Fig 2 Level_3 uncapped QC -> {plot.dir}")
