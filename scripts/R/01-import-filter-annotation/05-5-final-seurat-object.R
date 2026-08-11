# Create final filtered object for downstream analyses in manuscript
# Should match latest version of dataset on zenodo
# e.g. seurat_flex_filtered_v1.0.rds


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

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
results.dir <- "./results/05-cell-annotation"

## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")


# ---- Prepare objects ----

## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features.")
seurat.object <- qs2::qs_read("./data/processed/seurat/02-2-flex-dsUMI-s16-c92737.qs")
seurat.object$cell_barcode <- colnames(seurat.object)



## ---- Metadata: QC Cell Filter  -----
cli::cli_alert_info("Loading cell QC filtering metadata.")
cell.metadata <- readr::read_csv(file = "./results/05-cell-annotation/05-4-cell-metadata-annotations-levels.csv")
seurat.object <- AddMetaData(seurat.object, metadata = cell.metadata %>% column_to_rownames(var = "cell_barcode"))
seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )
seurat.object

# An object of class Seurat
# 19073 features across 64804 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 1 layer present: counts

### ---- ADD UMAPs ----
# add main umap for 64k cells

umap_df <- readr::read_csv("./results/04-reductions-snn/04-1-umap-m1v2-n16-Pds10k-C1-G1.csv")
seurat.object[["umap_AllSamples"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(tibble::column_to_rownames(umap_df, "cell_barcode")),
    # key = "umapAllSamples_",
    assay = DefaultAssay(seurat.object)
)

umap_df <- readr::read_csv("./results/05-cell-annotation/05-1-umap-nonneoplastic-filtered.csv")
seurat.object[["umap_NonNeoplastic"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(tibble::column_to_rownames(umap_df, "cell_barcode")),
    # key = "umapNN_",
    assay = DefaultAssay(seurat.object)
)

# ---- Factorize metadata ----
seurat.object <- .seuratFactorizeMdata(seurat.object)
table(seurat.object$Level_1, useNA = "ifany")
table(seurat.object$Level_2, useNA = "ifany")
table(seurat.object$Level_3, useNA = "ifany")
table(seurat.object$Level_3AC, useNA = "ifany")
table(seurat.object$Level_4, useNA = "ifany")
# mdata <- seurat.object@meta.data



# create "ACM" annotations with AC TE and NT in one cluster
mdata <- seurat.object@meta.data %>%
    mutate(Level_3ACM = if_else(Level_3AC %in% c("Astrocyte TE","Astrocyte NT"), "Astrocyte TE/NT", Level_3AC)) %>%
    mutate(Level_4ACM = if_else(Level_4 %in% c("Astrocyte TE","Astrocyte NT"), "Astrocyte TE/NT", Level_4))
table(mdata$Level_3ACM)
seurat.object@meta.data <- mdata
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)
table(seurat.object[["Level_3ACM"]])
table(seurat.object[["Level_4ACM"]])


# ---- Save Seurat Object ----
qs2::qs_save(seurat.object, "./data/processed/05-5-seurat-filtered-dsUMI-s16-c64804.qs")




