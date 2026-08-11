# 13-2 — Supp Fig 4 (Figure 4 supplement): per-cluster QC + neoplastic hypoxia
# ============================================================================
# Shows that neoplastic clusters, incl. Neopl-Bulk, are not QC artefacts, and
# that hypoxia is distributed across neoplastic clusters (not a discrete
# post-RT hypoxic state).
# Panels:
#   A  Per-cluster QC (nCount/nFeature/percent.mt) on RAW, UNCAPPED UMI across
#      ALL Level_4 clusters — neoplastic clusters sit within the full-tissue QC
# range, well above Neopl-RNA-low.
# (Level_3 / broad-cell-type companion of this panel, in 11-5 -> Fig S2)
# B HALLMARK_HYPOXIA module score per neoplastic cluster (violin).
# C Hypoxia score on the neoplastic UMAP (+ cluster annotation).
# NB Panel A reads the RAW (non-downsampled) object only for uncapped counts; it
# does NOT need the non-downsampled re-clustering (that is 11-3).
# ============================================================================

suppressMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(msigdbr); library(readr); library(glue); library(cli)
})
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
source("./scripts/R/00-0-source-functions.R")   # canonical .color_pal

plot.dir <- "./manuscript-figures/figure-4"
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

NEOPL_ORDER <- c("Neopl-Bulk","Neopl-CC-I","Neopl-CC-II","Neopl-CC-III",
                 "Neopl-OPC","Neopl-COP","Neopl-NC","Neopl-ACR",
                 "Neopl-ECM","Neopl-RNA-low")

# ---- object ----------------------------------------------------------------
cli::cli_alert_info("Loading UMI-downsampled, QC-filtered Seurat object")
obj <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
obj <- NormalizeData(obj, verbose = FALSE)

neopl <- subset(obj, subset = Level_4 %in% NEOPL_ORDER)
neopl$Level_4 <- factor(neopl$Level_4, levels = NEOPL_ORDER)
cli::cli_alert_info("Neoplastic cells: {ncol(neopl)}")

# ============================================================================
# A — per-cluster QC on RAW (uncapped) UMI, ALL Level_4 clusters
# ============================================================================
# The published object is UMI-downsampled (nCount capped at 10k), so its nUMI is
# flat-topped and uninformative. Pull the TRUE library size from the raw object,
# matched to the published Level_4 labels by (sample + core 10x barcode):
#   published : {sample}_{barcode}-1     raw : {barcode}-1_{mergeIndex}
# No re-clustering is needed — only per-cell counts under the published labels.
cli::cli_alert_info("Loading raw (non-downsampled) object for uncapped QC")
raw <- readRDS("./data/processed/seurat/seurat_flex_raw_v1.0.rds")
if (length(SeuratObject::Layers(raw, search = "counts")) > 1)
  raw <- SeuratObject::JoinLayers(raw)

pub_core <- mapply(function(cn, s) sub(paste0("^", s, "_"), "", cn),
                   colnames(obj), as.character(obj$orig.ident), USE.NAMES = FALSE)
ann      <- setNames(as.character(obj$Level_4), paste(as.character(obj$orig.ident), pub_core))
raw_key  <- paste(as.character(raw$orig.ident), sub("_\\d+$", "", colnames(raw)))
raw <- raw[, raw_key %in% names(ann)]
raw$Level_4 <- unname(ann[raw_key[raw_key %in% names(ann)]])
raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^mt-")  # mito genes present pre-filter
cli::cli_alert_info("Uncapped-QC cells (raw ∩ published): {ncol(raw)}")

lvl_order <- intersect(names(.color_pal[["Level_4"]]), unique(raw$Level_4))
qc_all <- raw@meta.data |>
  select(Level_4, nCount_RNA, nFeature_RNA, percent.mt) |>
  mutate(Level_4 = factor(Level_4, levels = lvl_order))
qc_long <- tidyr::pivot_longer(qc_all, c(nCount_RNA, nFeature_RNA, percent.mt),
                               names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric,
                         levels = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                         labels = c("UMI count (uncapped)", "Genes detected", "Mitochondrial %")))
p_qc <- ggplot(qc_long, aes(Level_4, value, fill = Level_4)) +
  geom_violin(scale = "width", alpha = 0.85, show.legend = FALSE) +
  geom_boxplot(width = 0.07, outlier.shape = NA, fill = "white", alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = .color_pal[["Level_4"]], guide = "none") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(title = "Per-cluster QC across all Level_4 states (RAW, uncapped UMI)",
       subtitle = "Non-UMI-downsampled counts; every annotated cluster ordered by compartment (non-neoplastic → neoplastic)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot.dir, "Figure-S4-QC-uncapped-all-Level4-violins.pdf"), p_qc, width = 14, height = 10)
write_csv(qc_all, file.path(plot.dir, "Figure-S4-QC-uncapped-all-Level4-violins-sourcedata.csv"))
rm(raw); gc()

# ============================================================================
# B + C — HALLMARK_HYPOXIA module score
# ============================================================================
hypoxia_genes <- msigdbr(species = "Mus musculus", category = "H") |>
  filter(gs_name == "HALLMARK_HYPOXIA") |> pull(gene_symbol) |> unique()
hypoxia_found <- intersect(hypoxia_genes, rownames(obj))
cli::cli_alert_info("HALLMARK_HYPOXIA genes found: {length(hypoxia_found)}/{length(hypoxia_genes)}")

obj <- AddModuleScore(obj, features = list(hypoxia_found), name = "Hypoxia_score")
neopl$Hypoxia_score <- obj$Hypoxia_score1[colnames(neopl)]

p_hypoxia_vln <- ggplot(neopl@meta.data,
                        aes(x = Level_4, y = Hypoxia_score, fill = Level_4)) +
  geom_violin(scale = "width", alpha = 0.85, show.legend = FALSE) +
  geom_boxplot(width = 0.07, outlier.shape = NA, fill = "white",
               alpha = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = .color_pal[["Level_4"]], guide = "none") +
  labs(x = NULL, y = "HALLMARK_HYPOXIA module score",
       title = "Hypoxia program per neoplastic cluster") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(plot.dir, "Figure-S4-hypoxia-neoplastic-violin.pdf"),
       p_hypoxia_vln, width = 10, height = 5)
write_csv(neopl@meta.data[, c("Level_4", "Hypoxia_score")],
          file.path(plot.dir, "Figure-S4-hypoxia-neoplastic-violin-sourcedata.csv"))

obj_neopl <- subset(obj, subset = Level_4 %in% NEOPL_ORDER)
obj_neopl$Hypoxia_score <- neopl$Hypoxia_score[colnames(obj_neopl)]
p_annot <- DimPlot(obj_neopl, group.by = "Level_4", cols = .color_pal[["Level_4"]], pt.size = 0.5, label = TRUE,
                   label.size = 3, repel = TRUE) +
  labs(title = "Neoplastic cluster annotation", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
p_umap <- FeaturePlot(obj_neopl, features = "Hypoxia_score", pt.size = 0.5, order = TRUE) +
  scale_color_gradientn(colours = c("grey90", "#FEE08B", "#D9534F")) +
  labs(title = "HALLMARK_HYPOXIA score — neoplastic cells", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 11)
ggsave(file.path(plot.dir, "Figure-S4-hypoxia-UMAP.pdf"), p_annot | p_umap, width = 14, height = 6)
# source data: plotted UMAP coords + Level_4 + hypoxia score per neoplastic cell
sd_umap <- p_annot$data
sd_umap$Hypoxia_score <- obj_neopl$Hypoxia_score[rownames(sd_umap)]
write_csv(sd_umap, file.path(plot.dir, "Figure-S4-hypoxia-UMAP-sourcedata.csv"))

# summaries (hypoxia only; per-cluster QC medians now reported by 11-3)
neopl@meta.data |> group_by(Level_4) |>
  summarise(n_cells = n(),
            mean_hypoxia = round(mean(Hypoxia_score), 4),
            pct_hypoxia_pos = round(mean(Hypoxia_score > 0) * 100, 1), .groups = "drop") |>
  write_csv(file.path(plot.dir, "Figure-S4-hypoxia-summary.csv"))

cli::cli_alert_success("Supp Fig 4 (neoplastic hypoxia) -> {plot.dir}")
