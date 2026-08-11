

# ---- Visium HD Signatures ----
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



## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
plot.dir <- "./manuscript-figures/figure-1"
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data (Visium HD healthy and tumor) ----
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


