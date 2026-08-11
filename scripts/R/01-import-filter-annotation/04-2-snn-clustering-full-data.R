

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
results.dir <- "./results/04-reductions-snn"
dir.create(results.dir, showWarnings = F, recursive = T)

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


# ---- Add METADATA ----
cli::cli_alert_info("Loading cell QC filtering metadata.")
# Add fitering metadata to object

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
table(seurat.object@meta.data$final_cell_filter)

# ---- Filter data ----
# filter to th 64k cells

# An object of class Seurat
# 19073 features across 64804 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features
cli::cli_alert_info("Filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")
## ---- filter cell QC ----
seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )
seurat.object



# 64804.qs
# An object of class Seurat
# 19073 features across 64804 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 1 layer present: counts
# 1 dimensional reduction calculated: umap_04_1_filtered
cli::cli_alert_info("Pre-filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")


# --- Add rcas transgene binned counts and module score  ----
rcas.module.df <- readr::read_csv(file = "./results/03-qc/03-6-rcas-module-score.csv" )
stopifnot(identical(rcas.module.df$cell_barcode, colnames(seurat.object)))
seurat.object@meta.data <- seurat.object@meta.data %>%
  dplyr::left_join(rcas.module.df, by = c("cell_barcode")) %>%
  column_to_rownames(var = "cell_barcode")

colnames(seurat.object@meta.data)
seurat.object@meta.data$cell_barcode <- colnames(seurat.object)
colnames(seurat.object)



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
G1@features.remove %in% rownames(seurat.object)
cli::cli_alert_info("Applying gene filter")
seurat.object <- seurat.object %>% .seuratFilterFoo(filter.object = G1)
seurat.object
cli::cli_alert_info("Filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")
# ℹ Filtered Seurat object contains 64804 cells and 15687 features.
run.name <- "4_dsUMI_s16_c64804_SNN"
nfeatures <- 2000
dims_pca <- 30
# reduction.prefix <- "UMAP_04_1_filtered"


# print Dims of final object
cli::cli_alert_info("Final object dims")
print(dim(seurat.object))

# Normalize and find variable features
cli::cli_alert_info("Normalizing and finding variable features")
seurat.object <- seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nfeatures)

seurat.object <- seurat.object %>%
  Seurat::ScaleData(verbose = T) %>%
  Seurat::RunPCA(verbose=TRUE, npcs = dims_pca) # %>%

# Run SNN graph
seurat.object <- FindNeighbors(
  seurat.object,
  dims = 1:dims_pca,
  k.param = 60,
  prune.SNN = 1/15,
  graph.name = c("kNN_k60","SNN_p15")
)
graph.name <- "SNN_p15"
names(seurat.object@graphs)
res <- 0.5
seurat.object <- FindClusters(
  seurat.object,
  graph.name = graph.name,
  resolution = res,
  algorithm = 1
)
graph.cluster <- paste0(graph.name,"_res.",res)
seurat.object@meta.data$SNN_clusters_all <- sapply(
  as.numeric(
    seurat.object@meta.data[,graph.cluster]), function(x){
      if (x < 10) {
        return(paste0("0", x))
      } else {
        return(as.character(x))
      }
    })
table(seurat.object@meta.data$SNN_clusters_all)
# 01    02    03    04    05    06    07    08    09    10    11    12    13    14    15    16    17    18
# 13798  6412  5660  3762  3687  3341  3015  2898  2881  2700  2693  2360  2276  1566  1136  1049   985   972
# 19    20    21    22    23    24    25    26
# 897   689   525   423   381   328   229   141

# ---- Add cluster info to metadata ----
snn.df <- seurat.object@meta.data[,c("cell_barcode","SNN_clusters_all")]
# save csv with cluster assignments
write.csv(snn.df, file = file.path(results.dir, "04-2-SNN-clusters-all.csv"), row.names = FALSE, quote = FALSE)


