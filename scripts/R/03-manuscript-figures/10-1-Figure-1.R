
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

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
set.seed(169)



## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table
plot.dir <- "./manuscript-figures/figure-1"
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)


# ---- Prepare objects ----

## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features.")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
seurat.object <-  seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData(verbose = T) %>%
  .seuratFactorizeMdata(.)
colnames(seurat.object@meta.data)



# ---- Barplot samples vs type ----
p <- .barplot_stacked(
  plot_df = seurat.object@meta.data,
  group.var = "sample_type",
  color.var = "sample_id",
  my.pal = .color_pal[["sample_id"]]
)
p
fig.name <- "Figure-1-sample-id-vs-type-barplot"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))


# ---- UMAP -AllCells-Sample_id-Split ----
# UMAPs 64k Sample characteristics
# UMAP visualization (All cells in QC filtered data, 04-02)
# Color by Sample type

mdata <- seurat.object@meta.data
str(mdata)
reduction.name <- "umap_AllCells"
annot <- "sample_id"
Idents(seurat.object) <- annot
table(seurat.object@meta.data[[annot]])


# Healthy: keep BRAIN_* names, everything else "other"
seurat.object$Healthy <- ifelse(grepl("^BRAIN_",Idents(seurat.object)),
  as.character(Idents(seurat.object)), "other")
# Primary: keep TP_* names, everything else "other"
seurat.object$Primary <- ifelse(grepl("^TP_", Idents(seurat.object)),
  as.character(Idents(seurat.object)), "other")
# Recurrent: keep TR_* names, everything else "other"
seurat.object$Recurrent <- ifelse(grepl("^TR_", Idents(seurat.object)),
  as.character(Idents(seurat.object)), "other")


p1 <- DimPlot(
  # order = names(peRcebe::percebe_color_pal[[annot]]),
  pt.size = 0.8, alpha = 0.5,
  shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = "Healthy",
  cols = .color_pal[[annot]],
  #cols = my.pal,
  label = FALSE,
  raster = FALSE
)
p2 <- DimPlot(
  #order = names(peRcebe::percebe_color_pal[[annot]]),
  pt.size = 0.8, alpha = 0.5,
  shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = "Primary",
  cols = .color_pal[[annot]],
  label = FALSE,
  raster = FALSE
)
p3 <- DimPlot(
  #order = names(peRcebe::percebe_color_pal[[annot]]),
  pt.size = 0.8, alpha = 0.5,
  shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = "Recurrent",
  cols = .color_pal[[annot]],
  label = FALSE,
  raster = FALSE
)

fig.name <- "Figure-1-umap-AllCells-Sample_id-Split"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1+p2+p3, width = 34, height = 9)
readr::write_csv(
  cbind(
    Embeddings(seurat.object, reduction.name)[, 1:2],
    seurat.object@meta.data[, c("Healthy", "Primary", "Recurrent"), drop = FALSE]
  ),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ----- Violin transgenes sample type ----
table(seurat.object[["disease_state"]], useNA = "ifany")
str(seurat.object[["disease_state"]])
levels(seurat.object@meta.data[,"disease_state"])

Idents(seurat.object) <- "sample_type"
genes <- rev(c("Tg-RFP-1","Tg-RFP-2","Tg-hPDGFB-nHA"))
fig.name <- "Figure-1-Violin-Transgenes"

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
for(gene in genes){
  p <- Seurat::VlnPlot(
    fill.by = "ident",
    flip=T,
    seurat.object, pt.size = 0,
    stack = F,
    cols = .color_pal[["sample_type"]],
    features = gene,
    group.by = "sample_type"
  )
  print(p)
  #ggsave(pdf.name, p, width = 7, height = 7)
  }
dev.off()
readr::write_csv(
  cbind(
    FetchData(seurat.object, vars = genes, slot = "data"),
    seurat.object@meta.data[, "sample_type", drop = FALSE]
  ),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)


# ----- UMAPCAS call ----
## UMAP and Barplot
table(seurat.object@meta.data$rcas_call)
p <- .barplot_stacked(
  plot_df = seurat.object@meta.data,
  group.var = "sample_type",
  color.var = "rcas_call",
  my.pal = .color_pal[["rcas_call"]]
  )
fig.name <- "Figure-1-Barplot-rcas_call-by-sample_type"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))

reduction.name <- "umap_AllCells"
Idents(seurat.object) <- annot

p <- DimPlot(pt.size = 0.25,

  alpha = 1,shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = "rcas_call",
  cols = .color_pal[["rcas_call"]],
  label = TRUE,
  raster = FALSE
)
p
fig.name <- "Figure-1-umap-AllCells-rcas_call"
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 12, height = 9)
readr::write_csv(
  cbind(
    Embeddings(seurat.object, reduction.name)[, 1:2],
    seurat.object@meta.data[, "rcas_call", drop = FALSE]
  ),
  file.path(plot.dir, glue("{fig.name}-sourcedata.csv"))
)



# ---- VISUM HD GEX ----


# ---- Visium HD Signatures ----
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette
# 10x recommends the use of 8um binned data for analysis
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
require(glue)
require(cli)
require(hdf5r)
# pak::pkg_install("arrow")
require(arrow)
?SeuratData
?Seurat::Load10X_Spatial
library(SeuratWrappers)
library(Banksy)
options(future.globals.maxSize = 1000 * 1024^2)  # Set limit to 1000 MiB (1 GiB), adjust as needed


if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(spacexr)


## ---- Load Visium data ----
data.processed.paths = list(
  "visium-hd-healthy" = "./data/processed/seurat/seurat_visium-hd-healthy-16um_v1.0.rds",
  "visium-hd-tumor" = "./data/processed/seurat/seurat_visium-hd-tumor-16um_v1.0.rds"
)

sd.list <- list()
def.assay <- "Spatial.016um"
for (data.name in names(data.processed.paths)) {
  data.file <- data.processed.paths[[data.name]]
  sd.list[[data.name]] <- readRDS(file = data.file)
  DefaultAssay(sd.list[[data.name]]) <- def.assay
}




### ---- Olig2 ----
gene <- "Olig2"
v_healthy <- FetchData(sd.list[["visium-hd-healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["visium-hd-tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "visium-hd-healthy"
fig.name <- glue("Figure-1-{data.name}-{gene}")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()
# Tumor
data.name <- "visium-hd-tumor"
fig.name <- glue("Figure-1-{data.name}-{gene}")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()



### ---- Tg-hPDGFB-nHA ----
# 1) Pick a gene and compute one global max (robust to outliers)
gene <- "Tg-hPDGFB-nHA"
v_healthy <- FetchData(sd.list[["visium-hd-healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["visium-hd-tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "visium-hd-healthy"
fig.name <- glue("Figure-1-{data.name}-{gene}")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()

# Tumor
data.name <- "visium-hd-tumor"
fig.name <- glue("Figure-1-{data.name}-{gene}")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()




### ---- Tg-RFP ----

# 1) Pick a gene and compute one global max (robust to outliers)
gene <- "CUSTOMPROBE-RFP"
v_healthy <- FetchData(sd.list[["visium-hd-healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["visium-hd-tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "visium-hd-healthy"
fig.name <- glue("Figure-S1-{data.name}-{gene}")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()

# Tumor
data.name <- "visium-hd-tumor"
fig.name <- glue("Figure-S1-{data.name}-{gene}")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
pdf(pdf.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()





# ---- Supplemental  ----
## RFP1 vs RFP2
RFP1 <- FetchData(seurat.object, vars = "Tg-RFP-1", slot = "data")[,1]
RFP2 <- FetchData(seurat.object, vars = "Tg-RFP-2", slot = "data")[,1]
cor.test(RFP1, RFP2, method = "pearson")
df.rfp <- data.frame(sample_type = seurat.object@meta.data$sample_type, RFP1 = RFP1, RFP2 = RFP2)
p <- ggplot(data = df.rfp, aes(x = RFP1, y = RFP2)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  xlab("RFP1 Expression") +
  ylab("RFP2 Expression") +
  ggtitle("Correlation of RFP1 and RFP2 Expression") +
  facet_wrap(~ sample_type)

  # geom_smooth(method = "lm", color = "blue", se = FALSE)
p
pdf.name <- file.path(plot.dir, "Figure-S1-RFP1-vs-RFP2-correlation.pdf")
ggsave(pdf.name, p, width = 13, height = 5)
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-RFP1-vs-RFP2-correlation-sourcedata.csv"))

## ---- UMAP hPDGFRb Expression ----

fig.name <- "Figure-S1-UMAP-hPDGFB-expression"
reduction.name <- "umap_AllCells"

grep("Tg-", rownames(seurat.object), value = TRUE)
grep("rcas", colnames(seurat.object@meta.data), value = TRUE, ignore.case = TRUE)
gene <- "Tg-hPDGFB-nHA"
seurat.data <- Seurat::NormalizeData(seurat.object)

p <- .seuratFeaturePlotHexbin(
  seurat.object = seurat.object,
  feature=gene,
  reduction = reduction.name
)

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 12, height = 9)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))
#
# SCpubr::do_ViolinPlot(seurat.object, features ="Kcnn4", group.by = "Level_3AC")
# SCpubr::do_ViolinPlot(seurat.object, features ="Sox10", group.by = "Level_4")
# dev.new()

# RFP 1
fig.name <- "Figure-S1-UMAP-RFP-1-expression"
reduction.name <- "umap_AllCells"

grep("Tg-", rownames(seurat.object), value = TRUE)
grep("rcas", colnames(seurat.object@meta.data), value = TRUE, ignore.case = TRUE)
gene <- "Tg-RFP-1"
seurat.data <- Seurat::NormalizeData(seurat.object)

p <- .seuratFeaturePlotHexbin(
  seurat.object = seurat.object,
  feature=gene,
  reduction = reduction.name
)

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 12, height = 9)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))

# RFP 2
fig.name <- "Figure-S1-UMAP-RFP-2-expression"
reduction.name <- "umap_AllCells"

grep("Tg-", rownames(seurat.object), value = TRUE)
grep("rcas", colnames(seurat.object@meta.data), value = TRUE, ignore.case = TRUE)
gene <- "Tg-RFP-2"
seurat.data <- Seurat::NormalizeData(seurat.object)

p <- .seuratFeaturePlotHexbin(
  seurat.object = seurat.object,
  feature=gene,
  reduction = reduction.name
)

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 12, height = 9)
readr::write_csv(p$data, file.path(plot.dir, glue("{fig.name}-sourcedata.csv")))




## ---- rcas probes: expression counts binned ----
probeCalls <- readr::read_csv("./results//03-qc/03-4-rcas-transGenes-probeCalls_full_table.csv")
head(probeCalls)

probeCalls$sample_type <- seurat.object@meta.data[probeCalls$cell_name,"sample_type"]
probeCalls <- .FactorizeMdata(probeCalls)
p1 <- .barplot_stacked(probeCalls, group.var = "sample_type", color.var = "hPDGFB_bin", my.pal = .color_pal[["hPDGFB_bin"]])
p2 <- .barplot_stacked(probeCalls, group.var = "sample_type", color.var = "RFP_bin", my.pal = .color_pal[["hPDGFB_bin"]])


table(probeCalls$rcas_both_3bin)
p3 <- .barplot_stacked(probeCalls, group.var = "rcas_both_3bin", color.var = "hPDGFB_bin", my.pal = .color_pal[["hPDGFB_bin"]])
p4 <- .barplot_stacked(probeCalls, group.var = "rcas_both_3bin", color.var = "RFP_bin", my.pal = .color_pal[["hPDGFB_bin"]])

p <- cowplot::plot_grid(
  p1, p2,
  p3, p4,
  ncol = 2,
  labels = "AUTO"
)
pdf.name <- file.path(plot.dir, glue("Figure-S1-rcas_probeCalls_barplots.pdf"))
ggsave(pdf.name, p, width = 9, height = 11)
readr::write_csv(
  dplyr::bind_rows(
    A = p1$data, B = p2$data, C = p3$data, D = p4$data,
    .id = "panel"
  ),
  file.path(plot.dir, "Figure-S1-rcas_probeCalls_barplots-sourcedata.csv")
)



## B: how the two RFP probes agree (binned counts)
rfp.df <- table(probeCalls$RFP_1_bin, probeCalls$RFP_2_bin)

rfp.df <- as.data.frame(rfp.df)
names(rfp.df) <- c("RFP_1_bin", "RFP_2_bin", "n")

p <- ggplot(rfp.df,
  aes(x = RFP_1_bin, y = RFP_2_bin, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    trans = "log10"
  ) +
  labs(x = "RFP probe 1 (binned counts)",
    y = "RFP probe 2 (binned counts)",
    fill = "Cells",
    title = "Agreement of RFP probe 1 and 2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p


## split RFP by sample_type onto probeCalls --------------------------

## 1) Define bin levels explicitly so all plots share axes ---------------
rfp_levels <- c("count-00", "count-01", "count-02-05", "count-06+")

## 2) Make table per sample_type, keeping all bin combinations -----------
rfp.df <- probeCalls  %>%
  dplyr::filter(!is.na(sample_type)) %>%   # restrict to the final annotated cells (drop cells not in the 64,804 object; avoids a spurious "NA" facet)
  dplyr::mutate(
    RFP_1_bin   = factor(RFP_1_bin, levels = rfp_levels),
    RFP_2_bin   = factor(RFP_2_bin, levels = rfp_levels),
    sample_type = factor(sample_type, levels = c("Healthy","Primary","Recurrent"))
  ) %>%
  dplyr::count(sample_type, RFP_1_bin, RFP_2_bin, .drop = FALSE)
# .drop = FALSE keeps all combinations of factor levels (even n = 0)

## 3) Plot: faceted by sample_type, continuous white–yellow–red fill ----
p <- ggplot(rfp.df,
  aes(x = RFP_1_bin, y = RFP_2_bin, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1
    #trans = "log10"
  ) +
  labs(
    x = "RFP probe 1 (binned counts)",
    y = "RFP probe 2 (binned counts)",
    fill = "Cells",
    title = "Agreement of RFP probes by sample type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ sample_type)
p <- p +
  scale_fill_gradientn(
    colours = c("white", "#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"),
    limits  = c(0, 5000),   # cap the palette at 10k
    oob = scales::squish       # values >10k become 10k color
  )

pdf.name <- file.path(plot.dir, glue("Figure-S1-rcas-RFP-probe-agreement-by-sample-type.pdf"))
pdf(pdf.name, width = 8, height = 4)
p
dev.off()
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-rcas-RFP-probe-agreement-by-sample-type-sourcedata.csv"))


###  2. RFP_bin × hPDGFB_bin heatmaps
bin_levels <- c("count-00", "count-01", "count-02-05", "count-06+")

df.bin <- probeCalls %>%
  mutate(
    RFP_bin    = factor(RFP_bin,    levels = bin_levels),
    hPDGFB_bin = factor(hPDGFB_bin, levels = bin_levels)
  ) %>%
  dplyr::count(RFP_bin, hPDGFB_bin, .drop = FALSE)

p <- ggplot(df.bin,
  aes(x = RFP_bin, y = hPDGFB_bin, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  labs(
    x = "RFP (binned counts)",
    y = "hPDGFB (binned counts)",
    fill = "Cells",
    title = "Joint distribution of RFP and hPDGFB (all samples)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p +
  scale_fill_gradientn(
    colours = c("white", "#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"),
    limits  = c(0, 10000),   # cap the palette at 10k
    oob = scales::squish       # values >10k become 10k color
  )
p

df.bin <- probeCalls %>%
  dplyr::filter(!is.na(sample_type)) %>%   # restrict to the final annotated cells (drop cells not in the 64,804 object; removes the spurious "NA" facet)
  mutate(
    RFP_bin    = factor(RFP_bin,    levels = bin_levels),
    hPDGFB_bin = factor(hPDGFB_bin, levels = bin_levels),
    sample_type = factor(sample_type, levels = c("Healthy","Primary","Recurrent"))
  ) %>%
  dplyr::count(sample_type, RFP_bin, hPDGFB_bin, .drop = FALSE)

p <- ggplot(df.bin,
  aes(x = RFP_bin, y = hPDGFB_bin, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 3) +
  labs(
    x = "RFP (binned counts)",
    y = "hPDGFB (binned counts)",
    fill = "Cells",
    title = "Joint distribution of RFP and hPDGFB by sample type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ sample_type)
p <- p +
  scale_fill_gradientn(
    colours = c("white", "#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"),
    limits  = c(0, 5000),   # cap the palette at 10k
    oob = scales::squish       # values >10k become 10k color
  )
p
pdf.name <- file.path(plot.dir, glue("Figure-S1-rcas-RFP-hPDGFB-joint-distribution-by-sample-type.pdf"))
pdf(pdf.name, width = 8, height = 4)
p
dev.off()
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-rcas-RFP-hPDGFB-joint-distribution-by-sample-type-sourcedata.csv"))


# 3-way table: RFP_bin × hPDGFB_bin × rcas_both_3bin
rcas_all <- probeCalls %>%
  mutate(
    RFP_bin    = factor(RFP_bin,    levels = bin_levels),
    hPDGFB_bin = factor(hPDGFB_bin, levels = bin_levels)
  ) %>%
  dplyr::count(RFP_bin, hPDGFB_bin, rcas_both_3bin, .drop = FALSE)

rcas_all <- rcas_all %>%
  group_by(RFP_bin, hPDGFB_bin) %>%
  summarise(
    n_total   = sum(n),
    rcas_class = rcas_both_3bin[which.max(n)],
    frac_major = max(n) / n_total,
    .groups = "drop"
  )

p <- ggplot(rcas_all,
  aes(x = RFP_bin, y = hPDGFB_bin, fill = n_total)) +
  geom_tile(color = "white") +
  geom_text(aes(
    label = sprintf("%s\n%0.0f%%", rcas_class, frac_major * 100)
  ), size = 3) +
  labs(
    x = "RFP (binned counts)",
    y = "hPDGFB (binned counts)",
    fill = "Cells",
    title = "Final RCAS class per (RFP, hPDGFB) bin (all samples)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

rcas_major_all <- rcas_all %>%
  mutate(
    rcas_class = factor(rcas_class,
      levels = c("rcas_neg", "rcas_ambiguous", "rcas_pos")),
    label = sprintf("%s\nn=%d\n%0.0f%%",
      rcas_class, n_total, frac_major * 100)
  )


p <- ggplot(rcas_major_all,
  aes(x = RFP_bin, y = hPDGFB_bin)) +
  # base heatmap: cell counts
  geom_tile(aes(fill = n_total), color = "white") +
  # coloured outline: RCAS status
  geom_tile(aes(color = rcas_class), fill = NA, size = 1) +
  # text: class, n, %
  geom_text(aes(label = label), size = 3) +
  scale_color_manual(
    values = c(
      rcas_neg       = "#FAF7D2FF",
      rcas_ambiguous = "#D6BB3B",
      rcas_pos       = "#F00000"
    ),
    na.value = NA,
    name = "RCAS class\n(outline)"
  ) +
  labs(
    x = "RFP (binned counts)",
    y = "hPDGFB (binned counts)",
    fill = "Cells",
    title = "Final RCAS class per (RFP, hPDGFB) bin (all samples)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p +
  scale_fill_gradientn(
    colours = c("white", "#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"),
    limits  = c(0, 5000),   # cap the palette at 10k
    oob = scales::squish       # values >10k become 10k color
  )
p

# wrapped plots per sample_type
rcas.df <- probeCalls %>%
  mutate(
    RFP_bin     = factor(RFP_bin,    levels = bin_levels),
    hPDGFB_bin  = factor(hPDGFB_bin, levels = bin_levels),
    sample_type = factor(sample_type)
  ) %>%
  dplyr::count(sample_type, RFP_bin, hPDGFB_bin, rcas_both_3bin, .drop = FALSE) %>%
  group_by(sample_type, RFP_bin, hPDGFB_bin) %>%
  summarise(
    n_total   = sum(n),
    n_max     = max(n),
    # majority class as factor
    rcas_class = rcas_both_3bin[which.max(n)],
    .groups = "drop"
  ) %>%
  mutate(
    frac_major = if_else(n_total > 0, n_max / n_total, NA_real_),
    # handle empty bins + convert safely
    rcas_class = if_else(
      n_total > 0,
      as.character(rcas_class),
      NA_character_
    ),
    rcas_class = factor(rcas_class,
      levels = c("rcas_neg", "rcas_ambiguous", "rcas_pos")
    ),
    label = if_else(
      is.na(rcas_class),
      "0",
      sprintf("%s\nn=%d\n%0.0f%%", rcas_class, n_total, frac_major * 100)
    )
  )


p <- ggplot(rcas.df,
  aes(x = RFP_bin, y = hPDGFB_bin)) +
  geom_tile(aes(fill = n_total), color = "white") +
  geom_tile(aes(color = rcas_class), fill = NA, size = 1) +
  geom_text(aes(label = label), size = 3) +
  scale_color_manual(
    values = c(
      rcas_neg       = "#FAF7D2FF",
      rcas_ambiguous = "#D6BB3B",
      rcas_pos       = "#F00000"
    ),
    na.value = NA,
    name = "RCAS class\n(outline)"
  ) +
  labs(
    x = "RFP (binned counts)",
    y = "hPDGFB (binned counts)",
    fill = "Cells",
    title = "Final RCAS class per (RFP, hPDGFB) bin by sample type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ sample_type)


p <- p +
  scale_fill_gradientn(
    colours = c("white", "#FFFFB2", "#FECC5C", "#FD8D3C", "#E31A1C"),
    limits  = c(0, 5000),   # cap the palette at 10k
    oob = scales::squish       # values >10k become 10k color
  )
p

pdf.name <- file.path(plot.dir, glue("Figure-S1-rcas-RFP-hPDGFB-and-RCAS-distribution-by-sample-type.pdf"))
pdf(pdf.name, width = 18, height = 7)
p
dev.off()
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-rcas-RFP-hPDGFB-and-RCAS-distribution-by-sample-type-sourcedata.csv"))


## ---- rcas gene module ----
rcas_deg_df <- readr::read_csv(file = "./results/03-qc/03-6-rcas-findallmarkers.csv")
## ---- Caclulate rcas-pos module score ----
module.name <- "rcas_module_score"
gene.vec <- rcas_deg_df$gene[1:25]
seurat.object <- seurat.object %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = nrow(seurat.object)) %>%  # all genes (v5: Inf trips 1:Inf on large objects)
  ScaleData()



### ---- plot rcas module score ----

Idents(seurat.object) <- "sample_type"
fig.name <- "Figure-S1-Violin-rcas-score-sample-type"

#pdf(pdf.name, width = 8, height = 6)
  p <- Seurat::VlnPlot(
    fill.by = "ident",
    flip=T,
    seurat.object, pt.size = 0,
    stack = F,
    cols = .color_pal[["sample_type"]],
    features = "rcas_module_score",
    group.by = "sample_type"
  )
p <- p +
  ggtitle("RCAS module score") +
  geom_hline(yintercept = 0.1, linetype="dashed", color = "red") +
  geom_hline(yintercept = -0.1, linetype="dashed", color = "red")


ggsave(file.path(plot.dir, "Figure-S1-Violin-rcas-score.pdf"), p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-Violin-rcas-score-sourcedata.csv"))


### ---- violin rcas call ----
p <- Seurat::VlnPlot(
  fill.by = "ident",
  flip=T,
  seurat.object, pt.size = 0,
  stack = F,
  cols = .color_pal[["rcas_module_bin"]],
  features = "rcas_module_score",
  group.by = "rcas_module_bin"
  )
p <- p +
  ggtitle("RCAS module score") +
  geom_hline(yintercept = 0.1, linetype="dashed", color = "red") +
  geom_hline(yintercept = -0.1, linetype="dashed", color = "red")
table(seurat.object@meta.data$rcas_module_bin)

ggsave(file.path(plot.dir, "Figure-S1-violin-rcas-score_by-rcasModuleBin.pdf"), p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-violin-rcas-score_by-rcasModuleBin-sourcedata.csv"))


### ---- barplot sample type per bin ----
mdata <- seurat.object@meta.data
# mdata <- mdata %>% dplyr::filter(sample_type %in% c("Primary", "Recurrent"))
p <- .barplot_stacked(
  mdata, color.var = "sample_type", group.var = "rcas_module_bin",
  scaled.y = F, my.pal = .color_pal[["sample_type"]]
)
p
ggsave(file.path(plot.dir, "Figure-S1-barplot-sampleType_by-rcasModuleBin.pdf"), p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-barplot-sampleType_by-rcasModuleBin-sourcedata.csv"))


### ---- barplot sample type per bin ----
mdata <- seurat.object@meta.data
# mdata <- mdata %>% dplyr::filter(sample_type %in% c("Primary", "Recurrent"))
p <- .barplot_stacked(
  mdata, color.var = "sample_type", group.var = "rcas_call",
  scaled.y = F, my.pal = .color_pal[["sample_type"]]
)
p
ggsave(file.path(plot.dir, "Figure-S1-barplot-sampleType_by-rcasCall.pdf"), p, width = 5, height = 8)
readr::write_csv(p$data, file.path(plot.dir, "Figure-S1-barplot-sampleType_by-rcasCall-sourcedata.csv"))



### --- Clustered HEatmap Plot individual genes ----
cell.order <- seurat.object@meta.data %>%
  arrange(rcas_module_score) %>%
  rownames()

# 2) Extract scaled data for the genes x cells you will plot
#    (assumes DefaultAssay(seurat.object) is the assay you want)
mat <- FetchData(
  seurat.object,
  vars = gene.vec,
  cells = cell.order,
  layer = "scale.data"
)

# Now mat is a data.frame: cells × genes (64k × 25)
# We must transpose it so genes = rows
mat <- t(as.matrix(mat))   # now 25 × 64k

# sanity check: should be 25 x 64k-ish
print(dim(mat))
# [1] 25 64804  (for example)

## 3. Cluster genes only (25 rows) --------------------------------------
# Euclidean distance between genes (rows)
d_genes  <- dist(mat, method = "euclidean")   # 25 x 25 distance
hc_genes <- hclust(d_genes, method = "ward")
gene.order <- rownames(mat)[hc_genes$order]

## 4. DoHeatmap with clustered rows, fixed cell order -------------------
seurat.object[["all_samples"]] <- "all"

p <- DoHeatmap(
  object     = seurat.object,
  features   = gene.order,     # <- clustered genes
  group.by   = "all_samples",
  cells      = cell.order,     # <- ordered by rcas_module_score
  draw.lines = TRUE,
  raster     = FALSE
) +
  ggtitle(fig.name) +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")[2:8])
  )
# p
pdf.name <- file.path(plot.dir, paste0(fig.name, ".pdf"))
ggsave(pdf.name, p, width = 18, height = 5)
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(mat[gene.order, cell.order]), "gene"),
  file.path(plot.dir, paste0(fig.name, "-sourcedata.csv"))
)


p <- DoHeatmap(
  object     = seurat.object,
  features   = gene.order,     # <- clustered genes
  group.by   = "all_samples",
  cells      = cell.order,     # <- ordered by rcas_module_score
  draw.lines = TRUE,
  raster     = TRUE
) +
  ggtitle(fig.name) +
  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")[2:8])
  )
# p
pdf.name <- file.path(plot.dir, paste0(fig.name, "_rasterized.pdf"))
ggsave(pdf.name, p, width = 18, height = 5)
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(mat[gene.order, cell.order]), "gene"),
  file.path(plot.dir, paste0(fig.name, "_rasterized-sourcedata.csv"))
)


p3 <- .plot_annotation_grid(
  plot_annotations = c("rcas_call","rcas_module_bin","sample_type"),
  #sort_levels = as.character(x_levels),
  meta_data = seurat.object@meta.data[cell.order,],
  reverse = T,
  color_palettes = .color_pal, keep_input_order = T
)
p3
pdf.name <- file.path(plot.dir, paste0(fig.name, "-AnnotationBar.pdf"))
ggsave(pdf.name, p3, width = 18, height = 5)
readr::write_csv(
  tibble::rownames_to_column(
    seurat.object@meta.data[cell.order, c("rcas_call", "rcas_module_bin", "sample_type"), drop = FALSE],
    "cell_name"
  ),
  file.path(plot.dir, paste0(fig.name, "-AnnotationBar-sourcedata.csv"))
)



# meta ordered exactly as in the heatmap
meta_rcas <- seurat.object@meta.data[cell.order, , drop = FALSE] %>%
  dplyr::mutate(
    cell_index = seq_len(n()),      # x-position
    rcas_module_bin = factor(rcas_module_bin)
  )

# waterfall / lollipop style plot
p_rcas_score <- ggplot(meta_rcas,
  aes(x = cell_index, y = rcas_module_score,
    colour = rcas_module_bin)) +
  # segments
  geom_segment(aes(xend = cell_index, y = 0, yend = rcas_module_score),
    linewidth = 0.2) +
  # new: gray dots per cell
  geom_point(
    aes(x = cell_index, y = rcas_module_score),
    inherit.aes = FALSE,
    color = "gray70",
    alpha = 0.5,
    size = 0.8
  ) +
  scale_color_manual(values = .color_pal[["rcas_module_bin"]]) +
  labs(
    x = NULL,
    y = "RCAS module score",
    colour = "RCAS module bin"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_hline(yintercept = 0.1,  linetype="dashed", color = "red") +
  geom_hline(yintercept = -0.1, linetype="dashed", color = "red")

p_rcas_score
pdf.name <- file.path(plot.dir, paste0(fig.name, "-RCAS-module-score.pdf"))
ggsave(pdf.name, p_rcas_score, width = 18, height = 3)
readr::write_csv(
  meta_rcas[, c("cell_index", "rcas_module_score", "rcas_module_bin")],
  file.path(plot.dir, paste0(fig.name, "-RCAS-module-score-sourcedata.csv"))
)

