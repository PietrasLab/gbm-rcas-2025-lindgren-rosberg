# ---- Visium HD Data Processing ----
# Process visium-hd from raw (spaceranger h5) to R object
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette
# 10x recommends the use of _at_least_ 8um binned data for analysis. We settle for 8 and 16um in downstream analyses

# These Correspond to the seurat files on zenodo
# https://zenodo.org/records/21607378

# seurat_visium-hd-healthy-16um_v1.0.rds
# seurat_visium-hd-healthy-8um_v1.0.rds
# seurat_visium-hd-tumor-16um_v1.0.rds
# seurat_visium-hd-tumor-8um_v1.0.rds


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


## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
results.dir <- "./results/05-cell-annotation"


# ---- LOAD RAW DATA and save seurat objects ----
raw.data.paths <- list(
    healthy_8um = c("./data/processed/spaceranger/10x-visium-hd/HB_1089",8),
    healthy_16um = c("./data/processed/spaceranger/10x-visium-hd/HB_1089",16),
    tumor_8um = c("./data/processed/spaceranger/10x-visium-hd/TP_1083",8),
    tumor_16um = c("./data/processed/spaceranger/10x-visium-hd/TP_1083",16)
)

?Seurat::Load10X_Spatial
# h5.file <- list.files(raw.data.path, pattern = "feature_slice.h5", full.names = FALSE)

spatial.list <- list()
for (i in seq_along(raw.data.paths)) {
    analysis.name <- names(raw.data.paths)[i]
    seurat.object <- Seurat::Load10X_Spatial(
        data.dir = raw.data.paths[[i]][1],
        bin.size = as.numeric(raw.data.paths[[i]][2])
        )
    norm.data <- NormalizeData(seurat.object)

    seurat.file <- file.path("./data/processed/seurat", glue("01-visium-hd-{analysis.name}.qs"))
    qs::qsave(norm.data, file = seurat.file)
}


