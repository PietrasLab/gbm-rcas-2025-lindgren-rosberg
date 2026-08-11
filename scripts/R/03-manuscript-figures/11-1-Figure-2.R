
# ---- Prepare Env ----
require(conflicted)
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

set.seed(169)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table

plot.dir <- glue("./manuscript-figures/figure-2")
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




# ---- UMAP SNN clustering - All cells ----
table(seurat.object@meta.data$SNN_clusters_all)

reduction.name <- "umap_AllCells"
Idents(seurat.object) <- "SNN_clusters_all"
p <- DimPlot(
  pt.size = 1.2,
  alpha = 1,shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = "SNN_clusters_all",
  cols = .color_pal[["SNN_clusters_all"]],
  label = TRUE,
  raster = FALSE
)
p
fig.name <- "Figure-2-UMAP-AllCells-SNN_clusters_all"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 10, height = 8)
readr::write_csv(
  cbind(Embeddings(seurat.object, reduction.name)[, 1:2],
        seurat.object@meta.data[, "SNN_clusters_all", drop = FALSE]),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)

# ---- Barplot RCAS_call vs SNN clusters_all ----

# p <- .barplot_stacked(
#   seurat.object@meta.data, color.var = "rcas_call",
#   my.pal = .color_pal[["rcas_call"]],
#   group.var = "SNN_clusters_all",
#   )
# p
table(seurat.object@meta.data$Level_1, seurat.object@meta.data$sample_type)
p <- .barplot_stacked(
  seurat.object@meta.data, color.var = "sample_type",
  my.pal = .color_pal[["sample_type"]],
  group.var = "Level_1",
)
p
fig.name <- "Figure-2-Barplot-SampleType-vs-Level_1"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))

p1 <- .barplot_stacked(
  seurat.object@meta.data, color.var = "rcas_call",
  my.pal = .color_pal[["rcas_call"]],
  group.var = "SNN_Clusters_broad",
  )
p1
p2 <- .barplot_stacked(
  seurat.object@meta.data, color.var = "sample_type",
  my.pal = .color_pal[["sample_type"]],
  group.var = "SNN_Clusters_broad",
  )
p2

fig.name <- "Figure-S2-Barplot-RCAS_call-vs-SNN_clusters_broad"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1+p2, width = 11, height = 5)
readr::write_csv(
  dplyr::bind_rows(rcas_call = p1$data, sample_type = p2$data, .id = "panel"),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


p1 <- .barplot_stacked(scaled.y = T,
  seurat.object@meta.data, color.var = "rcas_call",
  my.pal = .color_pal[["rcas_call"]],
  group.var = "SNN_clusters_all",
)
p1
fig.name <- "Figure-S2-Barplot-RCAS_call-vs-SNN_clusters_all"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1, width = 9, height = 4.5)
readr::write_csv(p1$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))

p1 <- .barplot_stacked(
  seurat.object@meta.data, color.var = "sample_type",
  my.pal = .color_pal[["sample_type"]],
  group.var = "SNN_clusters_all",
)
p1
fig.name <- "Figure-S2-Barplot-sample_type-vs-SNN_clusters_all"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1, width = 9, height = 4.5)
readr::write_csv(p1$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


# ---- 2C: UMAP of SNN broad clusters ----
table(seurat.object@meta.data$SNN_Clusters_broad)

p <- DimPlot(pt.size = 1.2,
  alpha = 1,shuffle = T,
  seurat.object,
  reduction =  "umap_AllCells",
  group.by = "SNN_Clusters_broad",
  cols = .color_pal[["SNN_Clusters_broad"]],
  label = TRUE,
  raster = FALSE
)
p
fig.name <- "Figure-2-UMAP-AllCells-SNN_Clusters_broad"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 10, height = 8)
readr::write_csv(
  cbind(Embeddings(seurat.object, "umap_AllCells")[, 1:2],
        seurat.object@meta.data[, "SNN_Clusters_broad", drop = FALSE]),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ---- UMAP SNN_clusters_nonNeoplastic  ----
annot <- "SNN_clusters_nonNeoplastic"
plot.object <- seurat.object %>% subset(subset = Level_1 == "Non-Neoplastic")

Idents(plot.object) <- annot
table(Idents(plot.object))

p <- DimPlot(pt.size = 1.2,
  alpha = 1,shuffle = T,
  plot.object,
  reduction = "umap_NonNeoplastic",
  group.by = "SNN_clusters_nonNeoplastic",
  cols = .color_pal[["SNN_clusters_nonNeoplastic"]],
  label = TRUE,
  raster = FALSE
)
p

fig.name <- "Figure-2-UMAP-SNN_clusters_nonNeoplastic"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 10, height = 8)
readr::write_csv(
  cbind(Embeddings(plot.object, "umap_NonNeoplastic")[, 1:2],
        plot.object@meta.data[, "SNN_clusters_nonNeoplastic", drop = FALSE]),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ---- Alluvial SampleType - RCAS call - SNN_clusters_all - Level 1 ----
# ---- 2D: Alluvial (zoomed) L1 to nonneoplastic L2/L3 ----
### ---- Fix so annotations are not shared over multiple alluvials ----
mdata_alluvial <- seurat.object@meta.data %>%
  dplyr::mutate(all_cells = "all")
mdata_alluvial <- mdata_alluvial %>%
  mutate(
    L1 = Level_1,
    L2 = Level_2,
    L3 = Level_3,
    L3AC = Level_3AC,
    L4 = Level_4,
    L1 = paste0("L1_", L1),
    L2 = paste0("L2_", L2),
    L3 = paste0("L3_", L3),
    L3AC = paste0("L5_", L3AC),
    L4 = paste0("L4_", L4)
  ) # %>%

# Define a color palette for each level
alluvial_pal <- .color_pal
# alluvial_pal[["SNN_clusters_all"]] <- .color_pal[["SNN_clusters_all"]]
# alluvial_pal[["sample_type"]] <- .color_pal[["sample_type"]]

alluvial_pal[["L1"]] <- .color_pal[["Level_1"]]
alluvial_pal[["L2"]] <- .color_pal[["Level_2"]]
alluvial_pal[["L3"]] <- .color_pal[["Level_3"]]
alluvial_pal[["L3AC"]] <- .color_pal[["Level_3AC"]]
alluvial_pal[["L4"]] <- .color_pal[["Level_4"]]
names(alluvial_pal[["L1"]]) <- paste0("L1_", names(alluvial_pal[["L1"]]))
names(alluvial_pal[["L2"]]) <- paste0("L2_", names(alluvial_pal[["L2"]]))
names(alluvial_pal[["L3"]]) <- paste0("L3_", names(alluvial_pal[["L3"]]))
names(alluvial_pal[["L3AC"]]) <- paste0("L5_", names(alluvial_pal[["L3AC"]]))
names(alluvial_pal[["L4"]]) <- paste0("L4_", names(alluvial_pal[["L4"]]))

str(mdata_alluvial)
mdata_alluvial[["L1"]] <- factor(mdata_alluvial[["L1"]], levels = names(alluvial_pal[[ "L1"]]))
mdata_alluvial[["L2"]] <- factor(mdata_alluvial[["L2"]], levels = names(alluvial_pal[[ "L2"]]))
mdata_alluvial[["L3"]] <- factor(mdata_alluvial[["L3"]], levels = names(alluvial_pal[[ "L3"]]))
mdata_alluvial[["L3AC"]] <- factor(mdata_alluvial[["L3AC"]], levels = names(alluvial_pal[[ "L3AC"]]))
mdata_alluvial[["L4"]] <- factor(mdata_alluvial[["L4"]], levels = names(alluvial_pal[[ "L4"]]))




### ---- sampleType to  SNN_allsamples to L1 ----
# Color the SNN clustering split on SNN broad types
str(mdata_alluvial)
colnames(mdata_alluvial)
alluvial.name <- "sampleType to SNN_All to rcas_joint to L1"
annotations <- c("sample_type","rcas_call","SNN_clusters_all","L1")
# annotations <- c("sample_type","SNN_clusters_all")
annot <- "sample_type"
table(mdata_alluvial[[annot]])
table(mdata_alluvial[["SNN_clusters_all"]])
sum(table(mdata_alluvial[["SNN_clusters_all"]]))

g <- .create_alluvial_plot(
  data = mdata_alluvial,
  #filter_value = c("L1_Non-Neoplastic"),
  #filter_column = "Tx08_L1",
  annotations = annotations,
  title = alluvial.name,
  fill_by = annot,
  color_palette = alluvial_pal[[annot]]
)
dev.new()
g


pdf.name <- file.path(plot.dir, glue("Figure-2-alluvial.pdf"))
ggsave(pdf.name, g, width = 12, height = 8)
readr::write_csv(g$data, file.path(plot.dir, "Figure-2-alluvial-sourcedata.csv"))


# ---- 2D: Alluvial (zoomed) L1 to nonneoplastic L2/L3 ----
### ---- Fix so annotations are not shared over multiple alluvials ----
mdata_alluvial <- seurat.object@meta.data %>%
  dplyr::mutate(all_cells = "all")
mdata_alluvial <- mdata_alluvial %>%
  mutate(
    L1 = Level_1,
    L2 = Level_2,
    L3 = Level_3,
    L3AC = Level_3AC,
    L4 = Level_4,
    L1 = paste0("L1_", L1),
    L2 = paste0("L2_", L2),
    L3 = paste0("L3_", L3),
    L3AC = paste0("L5_", L3AC),
    L4 = paste0("L4_", L4)
  ) # %>%

# Define a color palette for each level
alluvial_pal <- .color_pal
# alluvial_pal[["SNN_clusters_all"]] <- .color_pal[["SNN_clusters_all"]]
# alluvial_pal[["sample_type"]] <- .color_pal[["sample_type"]]

alluvial_pal[["L1"]] <- .color_pal[["Level_1"]]
alluvial_pal[["L2"]] <- .color_pal[["Level_2"]]
alluvial_pal[["L3"]] <- .color_pal[["Level_3"]]
alluvial_pal[["L3AC"]] <- .color_pal[["Level_3AC"]]
alluvial_pal[["L4"]] <- .color_pal[["Level_4"]]
names(alluvial_pal[["L1"]]) <- paste0("L1_", names(alluvial_pal[["L1"]]))
names(alluvial_pal[["L2"]]) <- paste0("L2_", names(alluvial_pal[["L2"]]))
names(alluvial_pal[["L3"]]) <- paste0("L3_", names(alluvial_pal[["L3"]]))
names(alluvial_pal[["L3AC"]]) <- paste0("L5_", names(alluvial_pal[["L3AC"]]))
names(alluvial_pal[["L4"]]) <- paste0("L4_", names(alluvial_pal[["L4"]]))

str(mdata_alluvial)
mdata_alluvial[["L1"]] <- factor(mdata_alluvial[["L1"]], levels = names(alluvial_pal[[ "L1"]]))
mdata_alluvial[["L2"]] <- factor(mdata_alluvial[["L2"]], levels = names(alluvial_pal[[ "L2"]]))
mdata_alluvial[["L3"]] <- factor(mdata_alluvial[["L3"]], levels = names(alluvial_pal[[ "L3"]]))
mdata_alluvial[["L3AC"]] <- factor(mdata_alluvial[["L3AC"]], levels = names(alluvial_pal[[ "L3AC"]]))
mdata_alluvial[["L4"]] <- factor(mdata_alluvial[["L4"]], levels = names(alluvial_pal[[ "L4"]]))


## ---- Alluvials non-neoplastic ----
### ---- non neoplastic regular alluvial ----
#### ---- Level 1 to Level 2 ----
str(mdata_alluvial)
alluvial.name <- "NonNeoplastic Level1 to Level2"
annotations <- c("L1","L2" )
annot <- "L2"
table(mdata_alluvial[[annot]])


g <- .create_alluvial_plot(
  data = mdata_alluvial,
  filter_value = c("L1_Non-Neoplastic"),
  filter_column = "L1",
  annotations = annotations,
  title = alluvial.name,
  fill_by = annot,
  color_palette = alluvial_pal[[annot]]
)
g
pdf.name <- file.path(plot.dir, glue("Figure-2-alluvial_L1_to_L2.pdf"))
ggsave(pdf.name, g, width = 12, height = 8)
readr::write_csv(g$data, file.path(plot.dir, "Figure-2-alluvial_L1_to_L2-sourcedata.csv"))


#### ----- Zoom NonNeoplstic Level2 to Level3 ----
str(mdata_alluvial)
alluvial.name <- "Zoom NonNeoplastic Level2 to Level3"
table(mdata_alluvial$Level_3, mdata_alluvial$Level_2)
table(mdata_alluvial$SNN_Clusters_broad, mdata_alluvial$rcas_call)
table(mdata_alluvial$SNN_Clusters_broad, mdata_alluvial$Level_1)


annot <- "Level_3"
filter_column = "Level_2"
table(mdata_alluvial[[annot]])
table(mdata_alluvial[[filter_column]])
g <- .create_zoom_alluvial_plot(
  data = mdata_alluvial,
  annot = annot,
  scale.all.ratio = 2,
  filter_column = filter_column,
  filter_value = c("Neural", "Immune","Vascular/Stromal"),
  title = alluvial.name,
  color_palette = alluvial_pal[[annot]],
  project_downward = T

)
g
pdf.name <- file.path(plot.dir, glue("Figure-2-alluvial_L2_to_L3.pdf"))
ggsave(pdf.name, g, width = 10, height = 10)
readr::write_csv(g$data, file.path(plot.dir, "Figure-2-alluvial_L2_to_L3-sourcedata.csv"))





# ---- Heatmap. TopGenes L3 non-neoplastic  ----
# Load DEGs
deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-NonNeoplastic-Level_3.csv")
deg_df <- .FactorizeMdata(deg_df)
str(deg_df$Level_3)
table(deg_df$Level_3)

topN <- 50
# annot <- "cluster"
topGenes <- deg_df %>%
  arrange(Level_3) %>%
  group_by(Level_3) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=50)
table(topGenes$cluster)
topGenes <- .FactorizeMdata(topGenes)


fig.name <- glue("Figure-2-Heatmap_FindMarkers_NonNeoplastic_Level_3_top50")

p <- Seurat::DoHeatmap(
  group.by = "Level_3",
  group.colors = .color_pal[["Level_3"]],
  draw.lines = TRUE,
  lines.width = 50,
  seurat.object %>% subset(Level_1 == "Non-Neoplastic"),
  features = topGenes$gene,
  label = FALSE
)
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
png.name <- file.path(plot.dir, glue("{fig.name}.png"))

p <- p + ggtitle(glue("{fig.name}")) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")[2:8]))
ggsave(pdf.name, p, width = 16, height = 12)
ggsave(png.name, p, width = 16, height = 12)
# Source data: per-cluster (Level_3) mean of the plotted scaled expression (Feature x cluster).
# The full per-cell DoHeatmap matrix is ~14M rows / ~1 GB — not repo-suitable.
readr::write_csv(
  p$data |> dplyr::filter(!is.na(Identity)) |>
    dplyr::group_by(Feature, Identity) |>
    dplyr::summarise(mean_scaled_expression = mean(Expression), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Identity, values_from = mean_scaled_expression),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))



# ---- DotPlot TopGenes L3 non-neoplastic Dotplot topGenes ----

deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-NonNeoplastic-Level_3.csv")
deg_df <- .FactorizeMdata(deg_df)
str(deg_df$Level_3)
table(deg_df$Level_3)
topN <- 3
topGenes <- deg_df %>%
  arrange(Level_3) %>%
  group_by(Level_3) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=topN)
table(topGenes$cluster)
topGenes <- .FactorizeMdata(topGenes)

p <- SCpubr::do_DotPlot(cluster = F,
  diverging.direction = 1,
  diverging.palette = "Blues",
  zscore.data = T,
  group.by = "Level_3", flip = F,
  seurat.object %>% subset(Level_1 == "Non-Neoplastic"),
  features = topGenes$gene)

p +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(
    expand = expansion(add = 0.5)  # shifts tick marks relative to grid centers
  )
p
fig.name <- file.path("Figure-2-DotPlot_FindMarkers_NonNeoplastic_Level_3_top3")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 14, height = 6.2)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


# ---- DotPlot TopGenes L1: Neoplastic vs NonNeoplastic  ----

deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Level_1.csv")
deg_df <- .FactorizeMdata(deg_df)
str(deg_df$Level_1)
table(deg_df$Level_1)
topN <- 25
topGenes <- deg_df %>%
  dplyr::filter(Level_1 == "Neoplastic") %>%
  #arrange(Level_1) %>%
  group_by(Level_1) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=topN)
table(topGenes$cluster)
topGenes <- .FactorizeMdata(topGenes)

fig.name <- file.path("Figure-2-DotPlot_FindMarkers_Neoplastic_top25")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 7, height = 7)
scCustomize::Clustered_DotPlot(
    colors_use_idents = .color_pal[["Level_3"]],
    seurat.object %>% subset(Level_1 != "Ambiguous"),
    features = topGenes$gene,
    colors_use_exp = c("#FFF5F7", "#FFE5EC", "#FAC5C3", "#FAA0A0", "#F44336", "#B71C1C"),
    group.by = "Level_3",
    flip = F
  )
dev.off()
# ggsave(pdf.name, p, width = 14, height = 6.2)
# Source data: Clustered_DotPlot draws via ComplexHeatmap (no ggplot $data). The plotted
# values are per-Level_3 scaled average expression + percent expressed; reproduce them
# with Seurat::DotPlot on the same subset/genes/grouping.
readr::write_csv(
  Seurat::DotPlot(
    seurat.object %>% subset(Level_1 != "Ambiguous"),
    features = topGenes$gene, group.by = "Level_3"
  )$data,
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)



# ---- Barplot sample_type Level_3  (Tumors) ----
# order based on proportion Recurrent per Level_3
mdata <- seurat.object@meta.data %>%
  dplyr::filter(Level_1 == "Non-Neoplastic") %>%
  .FactorizeMdata(.)
mdata.t <- mdata %>%
  dplyr::filter(disease_state == "Tumor") %>%
  .FactorizeMdata(.)



rec_order <- mdata.t %>%
  dplyr::count(Level_3, sample_type) %>%                 # counts per group
  group_by(Level_3) %>%
  dplyr::mutate(sum = sum(n)) %>%
  dplyr::mutate(prop = n / sum(n)) %>%                   # within-Level_3 proportions
  dplyr::filter(sample_type == "Recurrent") %>%          # keep Recurrent only
  # tidyr::complete(Level_3, fill = list(prop = 0)) %>%  # if some levels lack Recurrent
  # arrange(prop) %>%
  arrange(sum) %>%
  pull(Level_3)

# 2) relevel the factor
# mdata.t <- mdata.t %>%
#   mutate(Level_3 = factor(Level_3, levels = as.character(rec_order)))
mdata <- mdata %>%
  mutate(Level_3 = factor(Level_3, levels = as.character(rec_order)))

p1 <- .barplot_stacked(
  mdata %>% dplyr::filter(disease_state == "Tumor"),
  color.var = "sample_type",
  my.pal = .color_pal[["sample_type"]],
  group.var = "Level_3",

)
p1
p2 <- .barplot_stacked(
  mdata %>% dplyr::filter(disease_state == "Healthy"),
  color.var = "sample_type",
  my.pal = .color_pal[["sample_type"]],
  group.var = "Level_3",

)
p2

fig.name <- "Figure-2-Barplot-SampleType-vs-Level3"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1+p2, width = 8, height = 4)
readr::write_csv(
  dplyr::bind_rows(Tumor = p1$data, Healthy = p2$data, .id = "panel"),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ---- Barplots Level_3  ----
# order based on proportion Recurrent per Level_3
mdata <- seurat.object@meta.data %>%
  dplyr::filter(Level_1 != "Ambiguous") %>%
  .FactorizeMdata(.)

p1 <- .barplot_stacked(
  scaled.y = T,
  mdata,
  color.var = "Level_3",
  my.pal = .color_pal[["Level_3"]],
  group.var = "sample_type"
)
p1
p2 <- .barplot_stacked(
  mdata %>% dplyr::filter(Level_1 != "Neoplastic"),
  color.var = "Level_3",
  my.pal = .color_pal[["Level_3"]],
  group.var = "sample_type",
  scaled.y = T
)
p2


fig.name <- "Figure-2-Barplot-Level3-vs-SnampleType"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1+p2, width = 9, height = 6)
readr::write_csv(
  dplyr::bind_rows(all = p1$data, `non-neoplastic` = p2$data, .id = "panel"),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)



# ---- UMAP Level 3, All cells ----

seurat.object@meta.data$SNN_clusters_all
seurat.object@meta.data$SNN_clusters_nonNeoplastic

annot <- "Level_3"
fig.name <- glue("Figure-2G-UMAP-AllCells-Level_3")
reduction.name <- "umap_AllCells"
Idents(seurat.object) <- "Level_3"


p <- DimPlot(pt.size = 1.2,
  alpha = 1,shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = "Level_3",
  cols = .color_pal[["Level_3"]],
  label = TRUE,
  raster = FALSE
)
p

fig.name <- "Figure-2-UMAP-AllCells-Level_3"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 10, height = 8)
readr::write_csv(
  cbind(Embeddings(seurat.object, "umap_AllCells")[, 1:2],
        seurat.object@meta.data[, "Level_3", drop = FALSE]),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ----- Tabulate Level 3 -----

# Raw counts
tab <- table(seurat.object@meta.data$Level_3, seurat.object@meta.data$sample_type)
# Column-wise percentages
tab_colpct <- prop.table(tab, margin = 2) * 100
# Round and view
round(tab_colpct, 1)

seurat.sub <- subset(seurat.object, Level_1 == "Non-Neoplastic")
tab <- table(seurat.sub@meta.data$Level_3, seurat.sub@meta.data$sample_type)
# Column-wise percentages
tab_colpct <- prop.table(tab, margin = 2) * 100
# Round and view
round(tab_colpct, 1)




# ---- Supp Figures ----
mdata <- subset(seurat.object, subset = Level_1 == "Non-Neoplastic")@meta.data
# mdata <- mdata %>% dplyr::filter(sample_type %in% c("Primary", "Recurrent"))
p <- .barplot_stacked(
  mdata, color.var = "sample_type", group.var = "Level_3", scaled.y = T, my.pal = .color_pal[["sample_type"]]
)
p
fig.name <- glue("Figure-S2-StackedBarplot-NonNeoplastic")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 6, height = 8)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


