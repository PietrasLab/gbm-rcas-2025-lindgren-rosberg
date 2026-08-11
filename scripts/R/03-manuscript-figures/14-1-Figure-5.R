
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
set.seed(169)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table

plot.dir <- glue("./manuscript-figures/figure-5")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

# ---- Prepare objects ----

## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features with added metadata")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)
# add level annotations



# ----- Astro/Astro-neoplastic focused annotation ------
# Seelct only glia and neoplastic cells
annot <- "Level_4AC"
table(seurat.object$Level_2)
table(seurat.object$Level_4ACM)
pathways.object <- seurat.object
pathways.object <- subset(pathways.object, subset = Level_1 != "Ambiguous")
pathways.object <- subset(pathways.object, subset = Level_2 %in% c("Neural","Neoplastic"))
pathways.object <-  pathways.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nrow(pathways.object))  # all genes (v5: Inf trips 1:Inf on large objects)
colnames(pathways.object@meta.data)
pathways.object <- .seuratFactorizeMdata(pathways.object)
Idents(pathways.object) <- annot
table(Idents(pathways.object))

# make a Astro/Astro-neoplastic focused annotation
annot
astro.vec <- c("Astrocyte TE","Astrocyte NT","Astrocyte R")
mdata <- pathways.object@meta.data %>%
  mutate(
    annot = Level_3AC,                               # start with Level_3AC
    # merge TE + NT before other recoding
    annot = if_else(annot %in% c("Astrocyte TE","Astrocyte NT"),
      "Astrocyte TE/NT", annot),
    # recode non-neoplastic cells not in astro.vec → Non-Neoplastic (other)
    annot = case_when(
      Level_1 == "Non-Neoplastic" & !annot %in% c("Astrocyte TE/NT","Astrocyte R")
      ~ "Non-Neoplastic (other)",
      TRUE ~ annot
    )
  ) %>%
  mutate(
    annot = if_else(Level_4 == "Neopl-ACR", "Neoplastic ACR", annot)
  ) %>%
  # rename Neoplastic to Neoplastic (other)
  mutate(
    annot = if_else(annot == "Neoplastic",
      "Neoplastic (other)", annot)
  )
pathways.object$annot <- mdata$annot

# turn into factor with desired ordering
base_levels <- c("Non-Neoplastic (other)", "Astrocyte TE/NT","Astrocyte R","Neoplastic ACR","Neoplastic (other)")
other_levels <- base::setdiff(sort(unique(pathways.object$annot)), base_levels)
pathways.object$annot <- factor(pathways.object$annot, levels = c(base_levels, other_levels))
Idents(pathways.object) <- "annot"
table(Idents(pathways.object))

# Non-Neoplastic (other)        Astrocyte TE/NT            Astrocyte R
# 25235                   4397                    713
# Neoplastic ACR     Neoplastic (other)
# 443                  33284




# ----- UMAP ----
my.pal <- .color_pal[["Level_4ACM"]]
u <- which(!names(my.pal) %in% c("Astrocyte R","Neopl-ACR"))
my.pal[u] <- "#B3B3B3"
levels(seurat.object$"Level_4ACM") %in% names(my.pal)
#  c(
#   "Non-Neoplastic (other)" = "#BBD8B8",
#   "Astrocyte TE/NT" = "#BBD8B8", # "#82AC7C"
#   "Astrocyte R" = "#B8D53D",
#   "Neoplastic ACR" = "#F86814",
#   "Neoplastic (other)" = "#F4B7C7"
# )

p <- DimPlot(
  pt.size = 1, alpha = 1, shuffle = T,
  seurat.object,
  reduction = "umap_AllCells",
  group.by = "Level_4ACM",
  cols = my.pal,
  label = FALSE,
  raster = FALSE
)
p <- p + ggtitle(glue("ACR and Neoplastic ACR"))
p
ggsave(
  plot = p,
  filename = glue("{plot.dir}/Figure-5-UMAP-Level_4ACM-ACR-NeoplACR.pdf"),
  width = 11, height = 7
)
readr::write_csv(
  cbind(Embeddings(seurat.object, "umap_AllCells")[, 1:2],
        seurat.object@meta.data[, "Level_4ACM", drop = FALSE]),
  file.path(plot.dir, "Figure-5-UMAP-Level_4ACM-ACR-NeoplACR-sourcedata.csv")
)


# ----- Figure DotPlot DEGs ------
### ---- Load Astro reactive and Neopl-ACR signatures ----
deg_AC_vs_Other <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-NonNeoplastic-Level_3.csv") %>%
  dplyr::filter(cluster == "Astrocyte") %>%
  dplyr::arrange(p_val_adj) %>%
  slice_head(n=5)%>% pull(gene)
deg_ACR_vs_AstroTENT <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv") %>%
  dplyr::filter(cluster == "Astrocyte R") %>%
  dplyr::arrange(p_val_adj) %>%
  slice_head(n=5)%>% pull(gene)
# acr.sig.genes <-  c("S100a6", "Sulf2", "Cebpd", "S1pr3",
#   "S100a11", "Cd44", "Ifitm3", "C4b", "Gfap", "Vim")
deg_NACR_vs_OtherNeopl <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Neoplastic-Level_4.csv") %>%
  dplyr::filter(cluster == "Neopl-ACR") %>%
  dplyr::arrange(p_val_adj) %>%
  slice_head(n=10)%>% pull(gene)
deg_NACR_vs_Astro <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindMarkers-NeoplACR_vs_Astrocyte.csv") %>%
  dplyr::filter(cluster == "Neopl-ACR") %>%
  dplyr::arrange(p_val_adj) %>%
  slice_head(n=10) %>% pull(gene)
acr_acrl_genes <- unique(c(deg_AC_vs_Other, (deg_NACR_vs_OtherNeopl), rev(deg_ACR_vs_AstroTENT), rev(deg_NACR_vs_Astro)))

Idents(pathways.object) <- factor(
  pathways.object$annot,
  levels = c(
    "Neoplastic (other)",
    "Neoplastic ACR",
    "Astrocyte R",
    "Astrocyte TE/NT",
    "Non-Neoplastic (other)"
  )
)
table(pathways.object$annot, pathways.object$Level_4ACM)
p3 <- SCpubr::do_DotPlot(
  # use_viridis = T,
  # viridis.palette = "C",
  diverging.direction = 1,
  sequential.palette = "YlOrRd",
  plot.grid = T,
  flip=T,
  sample = pathways.object,
  features = rev(acr_acrl_genes)
  )

nx <- nlevels(Idents(pathways.object))
ny <- length(unique(acr_acrl_genes))

p3 <- p3 +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_discrete(position = "right", expand = expansion(add = 0.5)) +
  geom_vline(xintercept = if (nx > 1) seq(1.5, nx - 0.5, by = 1) else numeric(0),
    color = "grey75", linewidth = 0.4) +
  geom_hline(yintercept = if (ny > 1) seq(1.5, ny - 0.5, by = 1) else numeric(0),
    color = "grey85", linewidth = 0.3) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p3

  fig.name <- glue("Figure-5-ACR-NACR-topGenes-DEGs-DotHeatmap")
  pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
  p3 <- p3 + ggtitle(glue("{fig.name}"))
  ggsave(pdf.name, p3, width = 4.5, height = 12)
  readr::write_csv(p3$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


  # Flip it
  grp_lvls <- levels(Idents(pathways.object))

  p3 <- SCpubr::do_DotPlot(
    sequential.palette = "YlOrRd",
    plot.grid = TRUE,
    flip = FALSE,
    sample = pathways.object,
    features = acr_acrl_genes
  ) +
    scale_y_discrete(
      limits   = (grp_lvls),      # <<< reverse groups on y
      position = "right",
      expand   = expansion(add = 0.5)
    ) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    theme(
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor   = element_blank()
    )
  p3

  fig.name <- glue("Figure-5-ACR-NACR-topGenes-DEGs-DotHeatmap-flip")
  pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
  p3 <- p3 + ggtitle(glue("{fig.name}"))
  ggsave(pdf.name, p3, width = 12, height = 4.5)
  readr::write_csv(p3$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))



## ---- Pathway analysis ----

library(ComplexHeatmap)
library(circlize)
conflicts_prefer(dplyr::count)
matr <- readRDS(file.path("./results/07-2-gene-set-analysis/07-2-GeneSetAnalysis-AUCell-Richards-AC-NACR.rds"))
str(matr)
rownames(matr) <- gsub("\\.", "", rownames(matr))
## Canonical Richards display names: the two
## InHouse_BulkRNAseq_2019 lists are the FULL Richards dev/injury signature; the two
## *_Richards lists are the SHORT (top-250) version of the SAME signature. One axis, two depths.
rownames(matr) <- gsub("InHouse_BulkRNAseq_2019_DevelopmentalGSC",  "Richards_Developmental_2021",        rownames(matr))
rownames(matr) <- gsub("InHouse_BulkRNAseq_2019_InjuryResponseGSC", "Richards_InjuryResponse_2021",       rownames(matr))
rownames(matr) <- gsub("^Developmental_Richards$",                  "Richards_Developmental_2021_short",  rownames(matr))
rownames(matr) <- gsub("^Injury_Response_Richards$",                "Richards_InjuryResponse_2021_short", rownames(matr))
## --- assume you already built: matr, pathways.object, base_levels, sel_rows ---
sub.mat <- CalcStats(
  matr, f = Idents(pathways.object),
  order = "p", n = 10)
rownames(sub.mat)
sel_rows <- rownames(sub.mat)

# subset rows in your chosen order and z-score per row
matr_sel    <- matr[sel_rows, , drop = FALSE]
matr_sel_z  <- t(scale(t(matr_sel)))

# groups (boxes)
grp_vec <- Idents(pathways.object)[colnames(matr_sel_z)]
grp_vec <- factor(grp_vec, levels = base_levels)

# bring in per-cell annotations, aligned to matrix columns
meta <- pathways.object@meta.data[colnames(matr_sel_z), , drop = FALSE]

# make sure Level_3 and sample_type are factors with sensible (or existing) levels
if (!is.factor(meta$Level_3))     meta$Level_3     <- factor(meta$Level_3)
if (!is.factor(meta$sample_type)) meta$sample_type <- factor(meta$sample_type)
if (!is.factor(meta$Level_4ACM))  meta$Level_4ACM  <- factor(meta$Level_4ACM)

# (i) KEEP BOXES; (ii) ORDER inside each box by Level_3 then sample_type (no clustering)
col_order_by_group <- lapply(levels(grp_vec), function(g){
  idx <- which(grp_vec == g)
  # order by Level_3 then sample_type within this slice
  ord <- order(meta$Level_4ACM[idx], meta$sample_type[idx], na.last = TRUE)
  idx[ord]
})
col_order <- do.call(c, col_order_by_group)

# (ii) white → bright orange monochrome (adjust anchors if you want more contrast)
col_fun <- colorRamp2(c(-2, 0, 2), c("#FFFFFF", "#FFD28A", "#FF4A00"))

# (iii) top annotation bars
# reuse your palettes if available; otherwise auto-generate
get_pal <- function(x, fallback = NULL){
  ux <- levels(x)
  if (!is.null(fallback)) return(fallback[ux])
  # auto palette (distinct hues)
  cols <- grDevices::hcl.colors(length(ux), "Sunset")
  setNames(cols, ux)
}

col_Level3    <- if (exists(".color_pal") && !is.null(.color_pal[["Level_3"]]))     .color_pal[["Level_3"]]     else get_pal(meta$Level_3)
col_sample    <- if (exists(".color_pal") && !is.null(.color_pal[["sample_type"]])) .color_pal[["sample_type"]] else get_pal(meta$sample_type)
col_Level4ACM <- if (exists(".color_pal") && !is.null(.color_pal[["Level_4ACM"]]))  .color_pal[["Level_4ACM"]]  else get_pal(meta$Level_4ACM)

ha_top <- ComplexHeatmap::HeatmapAnnotation(
  Level_3     = meta$Level_3,
  sample_type = meta$sample_type,
  Level_4ACM  = meta$Level_4ACM,
  col = list(
    Level_3     = col_Level3,
    sample_type = col_sample,
    Level_4ACM  = col_Level4ACM
  ),
  annotation_name_side = "left",
  gp = grid::gpar(col = NA) # no borders on anno tiles
)

# main heatmap
ht <- ComplexHeatmap::Heatmap(
  matr_sel_z,
  name = "z-score",
  col = col_fun,
  column_split     = grp_vec,   # five boxes
  column_order     = col_order, # your custom order inside each box
  cluster_columns  = FALSE,     # no clustering; we set order
  cluster_rows     = TRUE,     # keep your row order (sel_rows)
  show_row_names   = TRUE,
  show_column_names= FALSE,
  top_annotation   = ha_top,
  border = TRUE,
  row_split = 5,
  use_raster = TRUE, raster_device = "png", raster_quality = 2
)

draw(ht, heatmap_legend_side = "right", merge_legends = TRUE)

fig.name <- glue("Figure-5-GeneSetAnalysis-AUCell-NeftelRichads-Heatmap-ACR-focus")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 18, height = 6)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()
# Source data: per-cluster (Level_4ACM) mean of the plotted z-scores (gene sets x cluster).
# The full per-cell z-matrix is ~34 MB — the per-cluster mean is the interpretable content.
local({
  g  <- as.character(meta$Level_4ACM)[col_order]
  mz <- matr_sel_z[, col_order, drop = FALSE]
  gm <- sapply(sort(unique(g)), function(k) rowMeans(mz[, g == k, drop = FALSE]))
  readr::write_csv(tibble::rownames_to_column(as.data.frame(gm), "geneset"),
                   file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))
})
readr::write_csv(
  data.frame(
    cell = colnames(matr_sel_z)[col_order],
    group = as.character(grp_vec[col_order]),
    Level_3 = as.character(meta$Level_3[col_order]),
    sample_type = as.character(meta$sample_type[col_order]),
    Level_4ACM = as.character(meta$Level_4ACM[col_order])
  ),
  file.path(plot.dir, glue("{fig.name}-colannot-sourcedata.csv"))
)




# ---- Figure. Boxplot cell numbers ----
### ---- Astrocytes  Tumor AstroOnly-----
annotations <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R", "Tumor ACR-sig")
# Process Metadata
mdata.plot <- seurat.object@meta.data %>%
  dplyr::filter(sample_type %in% c("Primary", "Recurrent")) %>%
  dplyr::filter(Level_4 %in% c("Astrocyte TE", "Astrocyte NT", "Astrocyte R",  "Neopl-ACR"))
mdata.plot <- .FactorizeMdata(mdata.plot)
#conflicts_prefer(base::sum)
p <- .barplot_stacked(
  my.pal = .color_pal[["Level_4"]],
  plot_df = mdata.plot,
  group.var = "sample_type",
  color.var = "Level_4",
  scaled.y = F,
  scaled.y.wihtin.color = F
)
p +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )+
  scale_y_discrete(expand = expansion(add = 100 ))

fig.name <- glue("Figure-5-Barplot_Astro_nCells_TPvTR")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))

ggsave(pdf.name, p, width = 5, height = 7)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


# One scaled vatiant as well
mdata.plot <- seurat.object@meta.data %>%
  dplyr::filter(Level_4 %in% c("Astrocyte TE", "Astrocyte NT", "Astrocyte R",  "Neopl-ACR"))
mdata.plot <- .FactorizeMdata(mdata.plot)
#conflicts_prefer(base::sum)
p <- .barplot_stacked(
  my.pal = .color_pal[["Level_4"]],
  plot_df = mdata.plot,
  group.var = "sample_type",
  color.var = "Level_4",
  scaled.y = T,
  scaled.y.wihtin.color = F
)
p +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )+
  scale_y_discrete(expand = expansion(add = 100 ))

fig.name <- glue("Figure-5-Barplot_Astro_Cells_TPvTR-scaled")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))

ggsave(pdf.name, p, width = 5, height = 7)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))
