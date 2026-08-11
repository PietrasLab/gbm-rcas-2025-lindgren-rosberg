
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


# ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")


## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
results.dir <- "./results/06-visium-hd"
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data (Visium HD healthy and tumor) ----
# use 16um bins here
data.obj.paths = list(
  "healthy" = "./data/processed/seurat/seurat_visium-hd-healthy-16um_v1.0.rds",
  "tumor" = "./data/processed/seurat/seurat_visium-hd-tumor-16um_v1.0.rds"

)

sd.list <- list()
def.assay <- "Spatial.016um"
for (data.name in names(data.obj.paths)) {
  data.file <- data.obj.paths[[data.name]]
  sd.list[[data.name]] <- qs::qread(file = data.file)
  DefaultAssay(sd.list[[data.name]]) <- def.assay
}


# ---- Basic QC -----
for( data.name in names(sd.list) ){
  object <- sd.list[[data.name]]
  cli::cli_alert_info("Processing {.var {data.name}}.")

  colnames(object@meta.data)
  table(object@meta.data$orig.ident)

  # note that many spots have very few counts, in-part
  # due to low cellular density in certain tissue regions
  pdf.name  <- file.path(results.dir, glue("06-1-nCount-plots-{data.name}.pdf"))
  pdf(pdf.name, width = 8, height = 4)
  p <- VlnPlot(object, features = "nCount_Spatial.016um", pt.size = 0) + theme(axis.text = element_text(size = 12)) + NoLegend()
  plot(p)
  p <- Seurat::SpatialFeaturePlot(object, features = "nCount_Spatial.016um") + theme(legend.position = "right")
  plot(p)
  dev.off()
}


# ---- Transcripts individual genes ----
# Custom plotting function

genes <- c(
  "Tg-hPDGFB-nHA",
  "CUSTOMPROBE-RFP",
  "Gfap"
)
stopifnot(all(genes %in% rownames(object)))




### ---- Tg-hPDGFB-nHA ----
# 1) Pick a gene and compute one global max (robust to outliers)
gene <- "Tg-hPDGFB-nHA"
v_healthy <- FetchData(sd.list[["healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "healthy"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object  = sd.list[[data.name]],
  gene.i  = gene,
  limits  = global_limits,   # identical for both plots
  slot    = "data"           # or "counts"/"scale.data" as you prefer
)
print(p); dev.off()

# Tumor
data.name <- "tumor"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()

### ---- RFP ----
# 1) Pick a gene and compute one global max (robust to outliers)
gene <- "CUSTOMPROBE-RFP"
v_healthy <- FetchData(sd.list[["healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "healthy"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object  = sd.list[[data.name]],
  gene.i  = gene,
  limits  = global_limits,   # identical for both plots
  slot    = "data"           # or "counts"/"scale.data" as you prefer
)
print(p); dev.off()

# Tumor
data.name <- "tumor"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()


### ---- Olig2 ----
# 1) Pick a gene and compute one global max (robust to outliers)
gene <- "Olig2"
v_healthy <- FetchData(sd.list[["healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "healthy"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()

# Tumor
data.name <- "tumor"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()



### ---- Gfap ----
# 1) Pick a gene and compute one global max (robust to outliers)
gene <- "Gfap"
v_healthy <- FetchData(sd.list[["healthy"]], vars = gene, slot = "data")[,1]
v_tumor   <- FetchData(sd.list[["tumor"]],   vars = gene, slot = "data")[,1]
# Use a robust cap so a few hot spots don’t blow out the legend
global_max <- as.numeric(quantile(c(v_healthy, v_tumor), probs = 0.99, na.rm = TRUE))
global_limits <- c(0, global_max)     # keep 0 on the low end color

# Healthy
data.name <- "healthy"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  #palette=c("#E4E3DB", "#BFE8DA", "#E3F6EF", "#FEE08B", "#F7B36B", "#EA8969", "#E07BA3", "#9E3D8D"),
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()

# Tumor
data.name <- "tumor"
file.name <- file.path(results.dir, glue("06-1-feature-{data.name}-{gene}.pdf"))
pdf(file.name, width = 8, height = 6)
p <- .spatialFeaturePlot(
  object = sd.list[[data.name]],
  gene.i = gene,
  limits = global_limits,
  slot   = "data"
)
print(p); dev.off()






