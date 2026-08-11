## Reductions, umaps etc

# ---- Env ----
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

set.seed(169)



# ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

# ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")


# ---- Load data -----
cli::cli_alert_info("Loading UMI-downsampled Seurat object with 92737 cells and 19073 features.")
seurat.object <- qs2::qs_read("./data/processed/seurat/02-2-flex-dsUMI-s16-c92737.qs")
colnames(seurat.object)
seurat.object@meta.data$cell_barcode <- colnames(seurat.object)
# 19073 features across 92737 samples within 1 assay


# ---- Load sample and probe feature metadata ----
pdata <- readRDS("./metadata/00-sample-metadata-flex-s16.Rds")
cli::cli_alert_info("Loading FlexProbes feature metadata.")
fdata <- read.delim("./metadata/00-probeset-flex.csv", header = T, sep=",")
str(fdata)
head(fdata)
gene.list <- readRDS(file = "./references/genesets/gene-list.Rds")



## ---- QC Filter Cell Metadata  -----
### defined in ./data/metadata/03-Cell-Filter-Dataframe.csv
#   percent_mito < 10
#   percent_haemoglobin < 5
#   nFeature_RNA > q05 & > q95 (sample adaptive)
#   nCount_RNA > q05 & > q95 (sample adaptive)
#   scDblFinder - remove cells flagged as non singlets (expected doublet rate 10%)
#   Balance cells from different samples - Downsample cells in the Tumor Recurrent group
filter.df <- readr::read_csv("./metadata/03-cell-filters.csv")
table(filter.df$final_cell_filter)
stopifnot(identical(filter.df$cell_barcode, colnames(seurat.object)))
seurat.object@meta.data$final_cell_filter <- filter.df$final_cell_filter



## --- Add transgene binned counts ----
# binned data from raw counts
# use as proxy for transgene expression
# use to obtina transgene/neoplastic gene module
# from 03-4
rcas.bin.df <- readr::read_csv("./results/03-qc/03-5-rcas-transGenes-probeCalls_rcas_both_3bin.csv")


# ---- QC Filter cells to final object   ----
# filter to th 64k cells
seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )


# ---- Gene filter  ----
# standard gene filter for reductions
# remove, y genes, mitochondria, haemoglobin, y chromosome
# also remove rcas flex probes to lessen the influence of these hugh expressing probes
# "Cst3" # identified as outlier - extreme high experressing in many cells and samples
G1 <- .CreateSeuratFilter(
  min.cells.expressed = 25,
  features.remove = c(
    gene.list$mitochondria,
    gene.list$haemoglobin,
    gene.list$rcas_flex_probes,
    gene.list$y_chr,
    "Cst3")
)

cli::cli_alert_info("Applying gene filter")
seurat.object <- seurat.object %>% .seuratFilterFoo(filter.object = G1)
seurat.object
cli::cli_alert_info("Filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")

nfeatures <- 2000
dims_pca <- 30
run.name <- "rcas_flex_n16"

# gsub multiple underscore with one underscore
cli::cli_alert_info(paste0("run.name: ", cli::col_br_yellow(run.name)))

# print Dims of final object
cli::cli_alert_info("Final object dims")
print(dim(seurat.object))
# --- 15687 genes, 64804 cells

# Normalize and find variable features
cli::cli_alert_info("Normalizing and finding variable features")
seurat.object <- seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nfeatures)

reduction.name <- "umap_04_1_filtered"
cli::cli_alert_info("Running PCA and UMAP with {dims_pca} dimensions and reduction name: {reduction.name}")
seurat.object <- seurat.object %>%
  Seurat::ScaleData(verbose = T) %>%
  Seurat::RunPCA(verbose=TRUE, npcs = dims_pca) %>%
  Seurat::RunUMAP(reduction.name = reduction.name,
    dims = 1:dims_pca, verbose = T)


# ---- Plot UMAP ----
colnames(seurat.object@meta.data)
cli::cli_alert_info("Plotting UMAP")
p <- Seurat::DimPlot(seurat.object, cols = .color_pal[["orig_ident"]],
  reduction = reduction.name, group.by = "orig.ident", label = TRUE, label.size = 3) +
  ggplot2::ggtitle("UMAP of Seurat object with 64804 cells") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p

.seura






