
# ---- Prepare Env ----
require(conflicted)
require(SeuratExtend)
require(DelayedMatrixStats)
require(tidyverse)
require(qs2)
require(cli)
require(glue)
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
# library(scCustomize)
# require(peRcebe)
require(CellMetaVerse)
require(kableExtra)
require(scCustomize)
require(scDblFinder)
require(BiocParallel)
library(reticulate)
library(dplyr)
library(patchwork)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::first)
conflicts_prefer(matrixStats::colRanks)
set.seed(169)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table

plot.dir <- glue("./manuscript-figures/figure-4")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)


# ---- Prepare objects ----

## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features with added metadata")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)
# add level annotations

seurat.object <-  seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData(verbose = T)
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)


# ---- Load Gene Signature Lists ----
richards_sigs <- readRDS(file.path("./references/genesets/Richards_NatCancer_2021_GeneSets_Mouse.rds"))


# ---- Figure 4 UMAP: NEoplastic subclsuters ------
fig.name <- "Figure-4-umap-snn-subclusters-neoplastic"
annot <- "Level_4"
reduction.name <- "umap_AllCells"

plot.object <- seurat.object

Idents(plot.object) <- annot
table(plot.object@meta.data[[annot]])

levels(plot.object) <- base::union(levels(plot.object), "other")
non_neoplastic <- WhichCells(plot.object, expression = Level_1 != "Neoplastic")
Idents(plot.object, cells = non_neoplastic) <- "other"
Idents(plot.object) <- droplevels(Idents(plot.object))
table(Idents(plot.object))

p <- DimPlot(pt.size = 1, alpha = 1,shuffle = T,
  plot.object,
  reduction = reduction.name,
  # group.by = annot,
  cols = .color_pal[["Level_4"]],
  label = TRUE,
  raster = FALSE
)
p <- p + ggtitle(glue("{fig.name}"))

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 13, height = 9)
# DimPlot uses Idents(plot.object) (non-neoplastic recoded to "other")
readr::write_csv(
  cbind(Embeddings(plot.object, reduction.name)[, 1:2],
        data.frame(ident = as.character(Idents(plot.object)),
                   row.names = colnames(plot.object))),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ---- Neoplastic RCAS Signature -----
# Signature From Fig 2
deg_df <- readr::read_csv("./results/03-qc/03-6-rcas-findallmarkers.csv")
deg_df <- .FactorizeMdata(deg_df)
#str(deg_df$Level_1)
#table(deg_df$Level_1)
topN <- 25
topGenes <- deg_df %>%
  arrange(p_val_adj) %>%
  slice_head(n=topN)
print(topGenes, n=topN)
table(topGenes$cluster)

sig.list <- list(neoplastic_rcas = topGenes$gene)
names(sig.list)
#module.scores <- lapply(deg.df.list, function(deg.df) {
seurat.sub <- seurat.object
# seurat.sub <- seurat.object %>% subset(SNN_Clusters_broad %in% c("Mixed Clusters","Neoplastic Clusters"))
object.tmp <- AddModuleScore(
  object = seurat.sub,
  features = sig.list,
  # ctrl = deg.df$ctrl,
  # assay = "Spatial.016um",
  name = glue("ModuleScore__")
)

colnames(object.tmp@meta.data)
module.score.df <- data.frame(object.tmp@meta.data[, grep(glue("ModuleScore__"), colnames(object.tmp@meta.data) )])
module_names <- paste0(glue("ModuleScore__"), names(sig.list))
module_names <- gsub(" |-", "_", module_names)
colnames(module.score.df) <- module_names
# saveRDS(module.score.df, file = file.path(plot.dir, "ModuleScores-Neoplastic.rds"))

for(col.i in module_names){
  object.tmp[[col.i]] <- module.score.df[, col.i]
}

p <- Seurat::VlnPlot(
  #split.by = "sample_type",
  pt.size = 0,
  object.tmp,
  cols = .color_pal[["Level_4"]],
  features = "ModuleScore__neoplastic_rcas",
  group.by = "Level_4"
)
print(p)

fig.name <- "Figure-4-rcas-Neoplastic-ModuleScore-Level4"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 9, height = 5)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))

colnames(seurat.object@meta.data)
p <- Seurat::VlnPlot(
  #split.by = "sample_type",
  pt.size = 0,
  seurat.object %>% subset(Level_1 == "Neoplastic"),
  cols = .color_pal[["Level_4"]],
  features = "nFeature_RNA",
  group.by = "Level_4"
)
print(p)




# ---- Neoplastic Cluster Signatures -----
annot <- "Level_4"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = nrow(seurat.subset)) %>%  # all genes (v5: Inf trips 1:Inf on large objects)
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Neoplastic-Level_4.csv")
deg_df <- deg_df %>%
  dplyr::rename(Level_4 = cluster) %>%
  .FactorizeMdata()
deg_df$Level_4


## ----  Get topGenes (for figure) ----
topN <- 50
topGenes <- deg_df %>%
  arrange(!!sym(annot)) %>%
  group_by(!!sym(annot)) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=50)
table(topGenes$Level_4)

# turn into list
sig.list <- split(topGenes, f = topGenes$Level_4)

sig.list <- lapply(sig.list, function(df){
  df$gene
})
# remove any signatures with less than 5 genes
sig.list <- sig.list[sapply(sig.list, length) >= 5]

names(sig.list)
#module.scores <- lapply(deg.df.list, function(deg.df) {
seurat.sub <- seurat.object %>% subset(Level_1 == "Neoplastic")
# seurat.sub <- seurat.object %>% subset(SNN_Clusters_broad %in% c("Mixed Clusters","Neoplastic Clusters"))
object.tmp <- AddModuleScore(
  object = seurat.sub,
  features = sig.list,
  # ctrl = deg.df$ctrl,
  # assay = "Spatial.016um",
  name = glue("ModuleScore__")
)

colnames(object.tmp@meta.data)
module.score.df <- data.frame(object.tmp@meta.data[, grep(glue("ModuleScore__"), colnames(object.tmp@meta.data) )])
module_names <- paste0(glue("ModuleScore__"), names(sig.list))
module_names <- gsub(" |-", "_", module_names)
colnames(module.score.df) <- module_names
# saveRDS(module.score.df, file = file.path(plot.dir, "ModuleScores-Neoplastic.rds"))

for(col.i in module_names){
  object.tmp[[col.i]] <- module.score.df[, col.i]
}

p <- Seurat::VlnPlot(
  #split.by = "sample_type",
  pt.size = 0,
  object.tmp,
  cols = .color_pal[["Level_4"]],
  features = module_names,
  group.by = "Level_4"
)
print(p)




# ---- CorrPlot CellTypes Level4 ------

## ---- Level 4, All cell types  ----
annot <- "Level_4"
plot.object <- seurat.object
plot.object <- subset(plot.object, subset = Level_1 != "Ambiguous")
Idents(plot.object) <- annot
table(Idents(plot.object))

plot.object <- plot.object %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

# p <- SCpubr::do_CorrelationPlot(sample = plot.object)
p <- SCpubr::do_CorrelationHeatmap(sample = plot.object)
p

#
fig.name <- glue("Figure-4-CorrPlot-CellTypes-Level4")
pdf.name <- file.path(plot.dir, paste0(fig.name, ".pdf"))
ggsave(pdf.name, p, width = 10, height = 10)
# do_CorrelationHeatmap $data holds the long-format pairwise correlations (tiles)
readr::write_csv(p$data, file.path(plot.dir, paste0(fig.name, "-sourcedata.csv")))


# Get the subtype labels from x axis - to add annotations
plot_data <- ggplot_build(p)
# OR (in case it's stored as a factor with levels)
x_levels <- plot_data$layout$panel_params[[1]]$x$breaks
meta_data <- data.frame(Level_4 = x_levels)
p <- .plot_annotation_grid(
  plot_annotations = c("Level_4"),
  sort_levels = as.character(x_levels),
  meta_data = meta_data,
  reverse = F,
  color_palettes = .color_pal
)

pdf.name <- file.path(plot.dir, paste0(fig.name, "-AnnotationBar.pdf"))
ggsave(pdf.name, p, width = 10, height = 5)
# Annotation bar plots the ordered Level_4 category strip (x-axis of the corr heatmap)
readr::write_csv(meta_data, file.path(plot.dir, paste0(fig.name, "-AnnotationBar-sourcedata.csv")))



# ---- Heatmaps Neoplastic Markers  -----
annot <- "Level_4"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = nrow(seurat.subset)) %>%  # all genes (v5: Inf trips 1:Inf on large objects)
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Neoplastic-Level_4.csv")
deg_df <- deg_df %>%
  dplyr::rename(Level_4 = cluster) %>%
  .FactorizeMdata()
deg_df$Level_4


## ----  Get topGenes (for figure) ----
topN <- 50
topGenes <- deg_df %>%
  arrange(!!sym(annot)) %>%
  group_by(!!sym(annot)) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=50)
table(topGenes$Level_4)

topGenes.csv <- file.path(plot.dir, glue("Figure-4-Table-TopGenes-Neoplastic.csv"))
readr::write_csv(topGenes, topGenes.csv)

fig.name <- "Figure-4-TopGenes-SeuratDoHeatmap-Neoplastic-Level4"
dev.new()
p <- Seurat::DoHeatmap(
  group.by = "Level_4",
  group.colors = .color_pal[[annot]],
  draw.lines = TRUE,
  lines.width = 50,
  seurat.subset,
  features = topGenes$gene,
  label = FALSE
)
p
str(p)

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
p2 <- p + ggtitle(glue("{fig.name}")) + NoLegend()
p2 <- p + ggtitle(glue("{fig.name}")) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")[2:8]))
ggsave(pdf.name, p2, width = 24, height = 5)
png.name <- file.path(plot.dir, glue("{fig.name}.png"))
ggsave(png.name, p2, width = 24, height =6)
# Source data: per-cluster (Identity) mean of the plotted scaled expression (Feature x cluster).
# The full per-cell DoHeatmap matrix is ~13M rows / ~900 MB — not repo-suitable.
readr::write_csv(
  p$data |> dplyr::filter(!is.na(Identity)) |>
    dplyr::group_by(Feature, Identity) |>
    dplyr::summarise(mean_scaled_expression = mean(Expression), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Identity, values_from = mean_scaled_expression),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


# Get the subtype labels from x axis - to add annotations
plot_data <- ggplot_build(p)
# OR (in case it's stored as a factor with levels)
# x_levels <- plot_data$layout$panel_params[[1]]$x$breaks

p3 <- .plot_annotation_grid(
  plot_annotations = c("Level_4","sample_type"),
  #sort_levels = as.character(x_levels),
  meta_data = seurat.subset@meta.data,
  reverse = F,
  color_palettes = .color_pal
)
p3
pdf.name <- file.path(plot.dir, paste0(fig.name, "-AnnotationBar.pdf"))
ggsave(pdf.name, p3, width = 10, height = 5)
readr::write_csv(
  tibble::rownames_to_column(
    seurat.subset@meta.data[, c("Level_4", "sample_type"), drop = FALSE], "cell_name"),
  file.path(plot.dir, paste0(fig.name, "-AnnotationBar-sourcedata.csv"))
)





## ---- Neoplastic vs non-neoplastic Signature -----
genelist.neoplastic <- list()
topN <- 5
deg_df <- readr::read_csv("./results//07-1-deg-findmarkers/07-1-FindAllMarkers-Level_1.csv")
deg_df <- .FactorizeMdata(deg_df)
str(deg_df$Level_1)
table(deg_df$Level_1)
topGenes <- deg_df %>%
  dplyr::filter(Level_1 == "Neoplastic") %>%
  #arrange(Level_1) %>%
  group_by(Level_1) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=topN)
table(topGenes$cluster)
topGenes <- .FactorizeMdata(topGenes)
genelist.neoplastic[["neopl"]] <- topGenes$gene

## ---- Neoplastic SNN -----
topN <- 3
deg_df <- readr::read_csv("./results//07-1-deg-findmarkers/07-1-FindAllMarkers-Neoplastic-Level_4.csv")
deg_df <- deg_df %>%
  dplyr::rename(Level_4 = cluster) %>%
  .FactorizeMdata()
deg_df$Level_4
deg_df <- .FactorizeMdata(deg_df)
topGenes <- deg_df %>%
  arrange(Level_4) %>%
  group_by(Level_4) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
genelist.neoplastic[["snn"]] <- topGenes$gene


p <- SCpubr::do_DotPlot(cluster = F,
  diverging.direction = 1,
  diverging.palette = "Greens",
  zscore.data = T,
  group.by = "Level_4", flip = F,
  seurat.object %>% subset(Level_1 == "Neoplastic"),
  features = unlist(genelist.neoplastic)
  )

p <- p +
  scale_color_gradientn(
    colours = c("#FFF0F5", "#F48FB1", "#C2185B", "#7B1E4E"),
    name = "Z-Scored Avg. Exp."
  ) +
  scale_fill_gradientn(
    colours = c("#FFF0F5", "#F48FB1", "#C2185B", "#7B1E4E"),
    name = "Z-Scored Avg. Exp."
  ) +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(expand = expansion(add = 0.5))
p
fig.name <- file.path("Figure-4-DotPlot-FindMarkers-Neoplastic-Level_4-top3")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 14, height = 5.2)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))



# ---- Figure 4D: Barplot ----
annot <- "Level_4"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

fig.name <- glue("Figure-4-Barplots-Neoplastic-Level_4")
Idents(seurat.subset) <- annot
table(Idents(seurat.subset))
mdata.plot <- tibble(seurat.subset@meta.data)

p1 <- .barplot_stacked(
  my.pal = .color_pal[[annot]],
  plot_df = mdata.plot,
  group.var = "sample_type",
  color.var = annot,
  scaled.y = T,
  scaled.y.wihtin.color = F
)
p1
p1 <- p1 + ggtitle(glue("{fig.name}\"All cells"))
# combine plots side by side


pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1, width = 5, height = 7)
readr::write_csv(p1$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))



# ---- Figure Pathways Neoplastic ----
require(SeuratExtend)

annot <- "Level_4"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)
Idents(seurat.subset) <- annot
table(Idents(seurat.subset))


### ---- Neftel + Richards + Nomura GBM cell states  ----
## Richards: FULL developmental / injury-response lists (the short top-250 variant is not
## shown here). Provenance = Richards & Whitley 2021 Nat Cancer, their in-house bulk-RNAseq
## GSC developmental<->injury axis (Supp Table 6 == Supp Table 7 membership; proven identical).
## Displayed as Richards_{Developmental,InjuryResponse}_2021.
u <- grep("Neftel|InHouse_BulkRNAseq_2019", names(richards_sigs))
neftel_richards <- richards_sigs[u]

## Nomura 2025 GBM meta-programs collapsed to the paper's 9 cell states + cell cycle
## (Extended Data Fig. 3; RP/MIC/LQ technical MPs dropped, NEU/Stress each union two MPs).
nomura        <- readRDS("./references/genesets/Nomura_NatGenet_2025_GeneSets_Mouse.rds")
nom_map       <- attr(nomura, "mp_state_map")
nom_pfx       <- attr(nomura, "prefix")
nom_statelab  <- c("AC-like"="Nomura_AC","MES-like"="Nomura_MES","NPC-like"="Nomura_NPC",
                   "OPC-like"="Nomura_OPC","GPC-like"="Nomura_GPC","NEU-like"="Nomura_NEU",
                   "Hypoxia"="Nomura_Hypoxia","Stress"="Nomura_Stress",
                   "Cilia-like"="Nomura_Cilia","Cell cycle"="Nomura_CC")
nom_states <- nom_map %>%
  dplyr::filter(state != "Unassigned") %>%                 # drop RP / MIC / LQ
  dplyr::mutate(key = paste0(nom_pfx, mp))
nomura_state_sigs <- split(nom_states$key, nom_states$state) %>%
  lapply(function(k) unique(unlist(nomura[k], use.names = FALSE)))
names(nomura_state_sigs) <- nom_statelab[names(nomura_state_sigs)]

gene_sigs <- c(neftel_richards, nomura_state_sigs)
names(gene_sigs)

seurat.subset <- SeuratExtend::GeneSetAnalysis(seurat.subset, genesets = gene_sigs)
matr <- seurat.subset@misc$AUCell$genesets
str(matr)
rownames(matr) <- gsub("\\.", "", rownames(matr))
## canonical Richards display names (FULL dev/injury lists)
rownames(matr) <- gsub("InHouse_BulkRNAseq_2019_DevelopmentalGSC",  "Richards_Developmental_2021",  rownames(matr))
rownames(matr) <- gsub("InHouse_BulkRNAseq_2019_InjuryResponseGSC", "Richards_InjuryResponse_2021", rownames(matr))
# saveRDS(matr, file.path(results.dir, "07-2-GeneSetAnalysis-AUCell-NeftelRichads-Neoplastic.rds"))


my.method <- "zscore"
conflicts_prefer(matrixStats::count)


stats.df <-  CalcStats(
  matr, f = Idents(seurat.subset),
  #  "mean", "median", "zscore", "tscore", "p", or "logFC". Default: "zscore".
  method = my.method,
  order = "p", n = Inf)

brbg_colors <- rev(corrplot::COL2("BrBG"))
my.pal <- brbg_colors[25:175]
method = "spearman"
p <- ComplexHeatmap::Heatmap(cluster_columns = F,
  stats.df, name = method,
  col = my.pal,
  clustering_distance_rows = "euclidean",
  border = TRUE,
  row_split = 5,
  row_gap = grid::unit(5, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.rect(x = x, y = y, width = width, height = height,
      gp = grid::gpar(col = "black", fill = NA, lwd = 0.5))
  }
)
p
fig.name <- glue("Figure-4-GeneSetAnalysis-AUCell-NeftelRichads-Neoplastic-Heatmap")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 9, height = 8)
print(p)
dev.off()
# ComplexHeatmap of the AUCell z-score matrix (gene sets x Level_4)
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(stats.df), "geneset"),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)

p <- SeuratExtend::Heatmap(
  CalcStats(
    matr, f = Idents(seurat.subset),
    #  "mean", "median", "zscore", "tscore", "p", or "logFC". Default: "zscore".
    method = my.method,
    order = "p", n = Inf),
  color_scheme = "YlOrBr",
  lab_fill = my.method
)
p

p <- p +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  # add x adn y axes lines and tick maarks
  theme(
    # panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    axis.ticks = element_line(color = "black"),  # Add tick marks
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black", face = "bold")  # Make axis labels visible
  )
p
fig.name <- glue("Figure-4-GeneSetAnalysis-AUCell-NeftelRichads-Neoplastic-Heatmap-alt")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 9, height = 7)
# SeuratExtend::Heatmap $data: long-format plotted z-score matrix (gene sets x Level_4)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))





### ---- waterfall sample_type ----
annot <- "sample_type"
Idents(seurat.subset) <- annot

my.method <- "zscore"
p <- SeuratExtend::Heatmap(
  CalcStats(
    matr, f = Idents(seurat.subset),
    #  "mean", "median", "zscore", "tscore", "p", or "logFC". Default: "zscore".
    method = my.method,
    order = "p", n = Inf),
  color_scheme = "YlOrBr",
  lab_fill = my.method
)
p

p <- SeuratExtend::WaterfallPlot(
  matr, f = Idents(seurat.subset),
  ident.1 = "Recurrent",
  ident.2 = "Primary",
  style = "segment",
  len.threshold = ,
  top.n = Inf
)
p
# add x and y axis lines
p <- p +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  # add x adn y axes lines and tick maarks
  theme(
    # panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    axis.ticks = element_line(color = "black"),  # Add tick marks
    axis.text = element_text(color = "black"),  # Ensure axis text is visible
    axis.title = element_text(color = "black", face = "bold")  # Make axis labels visible
  )



p

fig.name <- glue("Figure-4-GeneSetAnalysis-AUCell-NeftelRichads-Neoplastic-TPvTR-Waterfall")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 6, height = 5)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


# ---- Supp Figures ----
mdata <- subset(seurat.object, subset = Level_1 == "Neoplastic")@meta.data
# mdata <- mdata %>% dplyr::filter(sample_type %in% c("Primary", "Recurrent"))
p <- .barplot_stacked(
  mdata, color.var = "sample_type", group.var = "Level_4", scaled.y = T, my.pal = .color_pal[["sample_type"]]
  )
fig.name <- glue("Figure-S4-StackedBarplot-Neoplastic")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 10, height = 5)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))
