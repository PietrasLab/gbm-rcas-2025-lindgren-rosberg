# 11-3 — Supp Fig 2 (Figure 2 supplement): downsampling robustness
# ============================================================================
# Clarifies downsampled vs raw counts and shows the conclusions are not shaped
# by UMI-downsampling.
# Panels:
#   A  Per-sample sequencing depth — RAW (uncapped) UMI/cell only:
#      motivates capping to a common depth.
#   B  Non-downsampled re-clustering (raw counts, same cells/filter/params) —
#      side-by-side UMAP (published UMI-downsampled vs non-downsampled, Level_4).
#   C  Cluster-concordance heatmap (new clusters vs published Level_4; ARI/NMI).
# NB the uncapped-UMI per-cluster QC (Supp Fig 4) lives in 13-2 — it needs
#    only the raw counts per published label, not this re-clustering.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(readr); library(tidyr); library(cli); library(patchwork)
})
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
source("./scripts/R/00-0-source-functions.R")   # canonical .color_pal
plot.dir <- "./manuscript-figures/figure-2"
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

COND <- c("Healthy", "Primary", "Recurrent")
CCOL <- .color_pal[["sample_type"]]   # canonical palette

# ---- A — per-sample RAW depth ----------------------------------------------
pub <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")   # UMI-downsampled
raw <- readRDS("./data/processed/seurat/seurat_flex_raw_v1.0.rds")        # uncapped

# match raw -> published cells by (sample + core 10x barcode)
pub_core <- mapply(function(cn, s) sub(paste0("^", s, "_"), "", cn),
                   colnames(pub), as.character(pub$orig.ident), USE.NAMES = FALSE)
pub_key  <- paste(as.character(pub$orig.ident), pub_core)
raw_key  <- paste(as.character(raw$orig.ident), sub("_\\d+$", "", colnames(raw)))
raw_nc   <- setNames(raw$nCount_RNA, raw_key)

df <- tibble(
  sample    = as.character(pub$orig.ident),
  condition = factor(pub$sample_type, levels = COND),
  raw       = as.numeric(raw_nc[pub_key])           # uncapped library size
) |> filter(!is.na(raw))

permed <- df |> group_by(sample, condition) |>
  summarise(median_raw = median(raw), n = n(), .groups = "drop") |> arrange(median_raw)
fold <- max(permed$median_raw) / min(permed$median_raw)
df$sample <- factor(df$sample, levels = permed$sample)

p_depth <- ggplot(df, aes(sample, raw, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.85, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, linewidth = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = CCOL, name = NULL) +
  coord_flip() +
  labs(title = "Per-sample sequencing depth (raw, uncapped)",
       subtitle = sprintf("Raw per-sample median UMI/cell spans %.1f-fold across samples — motivates UMI-downsampling to a common cap.", fold),
       x = NULL, y = "UMI per cell (log10)") +
  theme_bw(base_size = 9)
ggsave(file.path(plot.dir, "Figure-S2-downsampling-depth-by-sample.pdf"), p_depth, width = 9, height = 5)
write_csv(df, file.path(plot.dir, "Figure-S2-downsampling-depth-by-sample-sourcedata.csv"))
cli::cli_alert_success("Supp Fig 2 depth panel (raw only) -> {plot.dir}")

# ---- B/C — non-downsampled re-clustering + concordance -----------------------
# Re-derive clustering + QC from RAW (uncapped) counts on the SAME cells, gene
# filter (G1) and params (nfeat 2000 / 30 PCs / res 0.5) as the published object
# (04-1/04-2), to show the cell-state structure and the QC conclusions are not
# artefacts of UMI-downsampling.
NFEATURES <- 2000; DIMS_PCA <- 30; RES <- 0.5
gene.list <- readRDS("./references/genesets/gene-list.Rds")
ann_col <- "Level_4"
ann <- setNames(as.character(pub@meta.data[[ann_col]]), pub_key)   # key -> label

if (length(SeuratObject::Layers(raw, search = "counts")) > 1) raw <- SeuratObject::JoinLayers(raw)
keep <- which(raw_key %in% names(ann))
rr <- raw[, keep]
rr$pub_label <- factor(unname(ann[raw_key[keep]]))
cli::cli_alert_info("Re-clustering {ncol(rr)} raw cells (shared with published)")

# gene filter G1 (identical to 04-2) then re-run reductions + clustering on RAW counts
G1 <- .CreateSeuratFilter(
  min.cells.expressed = 25,
  features.remove = c(gene.list$mitochondria, gene.list$haemoglobin,
                      gene.list$rcas_flex_probes, gene.list$y_chr, "Cst3"))
rr <- rr |> .seuratFilterFoo(filter.object = G1)
rr <- rr |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(nfeatures = NFEATURES, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(npcs = DIMS_PCA, verbose = FALSE) |>
  RunUMAP(dims = 1:DIMS_PCA, reduction.name = "umap_nonds", verbose = FALSE) |>
  FindNeighbors(dims = 1:DIMS_PCA, verbose = FALSE) |>
  FindClusters(resolution = RES, algorithm = 1, verbose = FALSE)
rr$nonds_cluster <- rr$seurat_clusters

# concordance: new clusters vs published Level_4
new_chr <- as.character(rr$nonds_cluster); pub_chr <- as.character(rr$pub_label)
ari <- aricode::ARI(new_chr, pub_chr); nmi <- aricode::NMI(new_chr, pub_chr)
tab <- table(new = new_chr, pub = pub_chr)
cluster_major <- apply(tab, 1, function(r) colnames(tab)[which.max(r)])
acc <- mean(pub_chr == cluster_major[new_chr])

# --- Panel C: cluster-concordance heatmap ---
ct <- as.data.frame(table(new = new_chr, published = pub_chr)) |>
  dplyr::group_by(published) |> dplyr::mutate(frac = Freq / sum(Freq)) |> dplyr::ungroup()
pub_order <- sort(unique(pub_chr))
new_order <- names(cluster_major)[order(match(cluster_major, pub_order))]
ct$published <- factor(ct$published, levels = pub_order)
ct$new       <- factor(ct$new, levels = new_order)
p_heat <- ggplot(ct, aes(published, new, fill = frac)) + geom_tile() +
  scale_fill_viridis_c(name = "frac of\npublished type") +
  labs(title = sprintf("Cluster concordance — ARI %.2f | NMI %.2f | label acc %.0f%%",
                       ari, nmi, 100 * acc),
       x = sprintf("published %s", ann_col), y = "new non-ds cluster") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot.dir, "Figure-S2-downsampling-cluster-concordance.pdf"), p_heat, width = 10, height = 7)
readr::write_csv(ct, file.path(plot.dir, "Figure-S2-downsampling-cluster-concordance-sourcedata.csv"))

# --- Panel B: side-by-side UMAP (published UMI-downsampled vs non-downsampled) ---
pal <- .color_pal[["Level_4"]]
red_pub <- grep("umap", names(pub@reductions), value = TRUE, ignore.case = TRUE)[1]
p_pub <- DimPlot(pub, reduction = red_pub, group.by = ann_col, cols = pal,
                 label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Published (UMI-downsampled)")
p_lab <- DimPlot(rr, reduction = "umap_nonds", group.by = "pub_label", cols = pal,
                 label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Non-downsampled (raw counts)")
ggsave(file.path(plot.dir, "Figure-S2-downsampling-UMAP-raw-vs-downsampled.pdf"),
       p_pub + p_lab, width = 15, height = 7)
mk <- function(pp, panel, lab_col) {
  d <- pp$data; xy <- d[, 1:2]; names(xy) <- c("dim1", "dim2")
  data.frame(panel = panel, xy, label = as.character(d[[lab_col]]))
}
readr::write_csv(rbind(mk(p_pub, "published", ann_col), mk(p_lab, "non_downsampled", "pub_label")),
                 file.path(plot.dir, "Figure-S2-downsampling-UMAP-raw-vs-downsampled-sourcedata.csv"))

cli::cli_alert_success("Supp Fig 2 depth + non-ds re-clustering — ARI {round(ari,3)} / NMI {round(nmi,3)} / acc {round(100*acc)}%")
