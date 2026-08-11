## ---- Step 4: run Allen MapMyCells on non-Neoplastic subset ----
# Use the Allen Brain Atlas to annotate the non-neoplastic clusters
# RRID:SCR_024672
# https://portal.brain-map.org/atlases-and-data/bkp/mapmycells
# Map the entire dataset

### Save h5ad object to perform for MapMyCells
# You should now have a matrix called "count_matrix" a sample vector called "obs" and a gene vector called "var2"
# seurat get matrix raw counts


## UMAP of the 64804 cells  included after filtering.
## The UMAP is performed on gene-filtered data, 15688 genes
## nfeatures <- 2000
# dims_pca <- 30

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

# ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
results.dir <- "./results/05-cell-annotation"
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



## ---- Metadata: QC Cell Filter  -----
cli::cli_alert_info("Loading cell QC filtering metadata.")
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
seurat.object <- AddMetaData(seurat.object, metadata = filter.df %>% dplyr::select(cell_barcode, final_cell_filter) %>% column_to_rownames(var = "cell_barcode"))
table(seurat.object@meta.data$final_cell_filter)
cli::cli_alert_info("Pre-filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")
seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )
cli::cli_alert_info("Filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")
seurat.object

# n object of class Seurat
# 19073 features across 64804 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 1 layer present: counts

### --- Metadata: RCAS transgene   ----
rcas.module.df <- readr::read_csv(file = "./results/03-qc/03-6-rcas-module-score.csv" )
stopifnot(identical(rcas.module.df$cell_barcode, colnames(seurat.object)))
# seurat.object@meta.data <- seurat.object@meta.data %>%
#   dplyr::left_join(rcas.module.df, by = c("cell_barcode")) %>%
#   column_to_rownames(var = "cell_barcode")
seurat.object <- AddMetaData(
  seurat.object, metadata = rcas.module.df
    %>% column_to_rownames(var = "cell_barcode"))
colnames(seurat.object@meta.data)

### ---- Metadata: SNN clustering (full filtered dataset) ----
# SNN clustreing of all 64804 cells.
# These clusters, together with RCAS transgene expression annotations form the basis for cell type annotation.
snn.df <- readr::read_csv(file = "./results/04-reductions-snn/04-2-SNN-clusters-all.csv")
snn.df
stopifnot(identical(snn.df$cell_barcode, colnames(seurat.object)))

seurat.object <- AddMetaData(
  seurat.object, metadata = snn.df  %>%
    column_to_rownames(var = "cell_barcode"))
table(seurat.object@meta.data$SNN_clusters_all)

snn.broad <- readr::read_csv(file = "./results/05-cell-annotation/05-1-snn-clusters-broad.csv")
snn.broad
seurat.object <- AddMetaData(seurat.object, metadata = snn.broad  %>% column_to_rownames(var = "cell_barcode"))
table(seurat.object[["SNN_Clusters_broad"]])


### ---- Metadata: SNN clustering (non-neoplastic subset of dataset) ----
# ADD SNN clusters and Cluster Names for non-neoplastic cells (From Script 04-3)
# save csv with cluster assignments
snn.df <- readr::read_csv(file = "./results/05-cell-annotation/05-1-snn-clusters-non-neoplastic.csv")
stopifnot(all(snn.df$cell_barcode %in% colnames(seurat.object)))
seurat.object <- AddMetaData(seurat.object, metadata = snn.df  %>% column_to_rownames(var = "cell_barcode"))
colnames(seurat.object@meta.data)



### ---- ADD UMAPs ----
# add main umap for 64k cells
umap_df <- readr::read_csv("./results/04-reductions-snn/04-1-umap-m1v2-n16-Pds10k-C1-G1.csv")
seurat.object[["umap_AllSamples"]] <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(tibble::column_to_rownames(umap_df, "cell_barcode")),
  # key = "umapAllSamples_",
  assay = DefaultAssay(seurat.object)
)
# emb <- as.matrix(tibble::column_to_rownames(umap_df, "cell_barcode"))
# Use the canonical reduction name/key that both packages accept
# colnames(emb) <- c("umap_1","umap_2")
# seurat.object[["umap"]] <- Seurat::CreateDimReducObject(
#   embeddings = emb,
#   key       = "umap_",              # valid key
#   assay     = DefaultAssay(seurat.object)
# )
umap_df <- readr::read_csv("./results/05-cell-annotation/05-1-umap-nonneoplastic-filtered.csv")
seurat.object[["umap_NonNeoplastic"]] <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(tibble::column_to_rownames(umap_df, "cell_barcode")),
  # key = "umapNN_",
  assay = DefaultAssay(seurat.object)
)


# ---- Genreate hd5ad for Allen MapMyCells -----
## Use downsampled data, all genes, Non Neoplastic cells only,
sd.nonneoplastic <- subset(seurat.object, subset = rcas_module_bin == "rcas_module_neg")
# 19073 features across 32379 samples within 1 assay
sd.nonneoplastic <- subset(sd.nonneoplastic, subset = rcas_counts_bin != "rcas_pos")
# 19073 features across 32093 samples within 1 assay
table(sd.nonneoplastic@meta.data$SNN_Clusters_broad)
sd.nonneoplastic <- subset(sd.nonneoplastic, subset = SNN_Clusters_broad != "Neoplastic Clusters")
# 19073 features across 30345 samples within 1 assay

## ---- Create counts matrix  ----
count_matrix <- GetAssayData(sd.nonneoplastic, layer = "counts")
count_matrix <- t(count_matrix)
str(count_matrix)
head(count_matrix[1:10, 1:5])
colnames(sd.nonneoplastic@meta.data)
obs <- sd.nonneoplastic@meta.data[, c("orig.ident","sample_type")]
head(obs)
obs$cell_barcode <- rownames(obs)
obs <- obs %>% dplyr::select(cell_barcode, orig.ident, sample_type)
# Get probe feature data with ENSEMBL ids
var2 <- data.frame(symbol = colnames(count_matrix), stringsAsFactors = F)
str(var2) # 19073
str(fdata)
df <- dplyr::left_join(var2, fdata, by=c("symbol" = "external_gene_name"))
str(df)
df <- df %>% dplyr::select(ensembl_gene_id, symbol)
str(df)
identical(df$symbol, var2$symbol)
var2 <- df
colnames(count_matrix) <- var2$ensembl_gene_id
table(duplicated(var2$ensembl_gene_id)) # 24
# Some few probes have duplicated ensembl ids - remove
u <- duplicated(var2$ensembl_gene_id)
var2 <- var2[!u,]
count_matrix <- count_matrix[,!u]
# remove 1 gene w NA ensembl_gene_id
u <- which(is.na(var2$ensembl_gene_id))
var2[u,]
count_matrix <- count_matrix[, -u]
var2 <- var2[-u,]
rownames(var2) <- var2$ensembl_gene_id
dim(var2)
dim(count_matrix)
# [1] 30345 19059

## ---- Export via anndata object ----
# Use reticulate and conda to generate anndata object
# This to output your variable in a compressed h5ad file,
# check the size, and then split the output file into multiple files for upload to MapMyCells if the size exceeds 2GB.
reticulate::use_python("/replace/with/path/mambaforge/envs/gbm-rcas-pyminer/bin/python", required = TRUE)
reticulate::py_config()

# Import libraries
library(anndata)
class(count_matrix)
# Convert count matrix to a CSR (row-based access) sparse matrix and save in anndata format.
#count_matrix = scipy.sparse.csr_matrix(count_matrix)
str(count_matrix)
ad = AnnData(
  X=count_matrix,
  obs=obs,
  var=var2
)
ad
# Write to your location of choice using compression
output_path = './data/processed/anndata/05-2-flex-anndata-non-neoplastic-dsUMI10k-c30345.h5ad'
write_h5ad(ad, output_path, compression='gzip')


# ---- Import Allen MapMyCells results ----
# https://portal.brain-map.org/atlases-and-data/bkp/mapmycells
allen.celltypes <- readr::read_csv("./results/05-cell-annotation/05-2-allen-non-neoplastic-allen-HierarchicalMapping.csv", skip=4)
allen.celltypes # 30,345
table(allen.celltypes$class_name)
table(allen.celltypes$subclass_name)
# Create a allen_celltpe column
# stratify neurons only into Glut, GABA and Dopa
# for remaining use the subclass_name
library(stringr)
allen.celltypes <- allen.celltypes %>%
  mutate(
    allen_celltype = case_when(
      str_detect(class_name, "Glut") ~ "Glut",
      str_detect(class_name, "GABA") ~ "GABA",
      str_detect(class_name, "Dopa") ~ "Dopa",
      TRUE                           ~ subclass_name    # all other classes: use subclass
    ),
    # safety fallback if subclass_name is missing/blank
    allen_celltype = if_else(is.na(allen_celltype) | allen_celltype == "",
      class_name, allen_celltype)
  )
table(allen.celltypes$allen_celltype)

## Save allen_celltype metadata
allen.meta <- allen.celltypes %>%
  dplyr::select(cell_id, allen_celltype) %>%
  dplyr::rename(cell_barcode = cell_id)
head(allen.meta)
readr::write_csv(allen.meta, "./metadata/05-2-non-neoplastic-allen-mapping-proc.csv")


## add Allen celltypes to Seurat object ----
stopifnot(identical(allen.meta$cell_barcode, colnames(sd.nonneoplastic)))
sd.nonneoplastic <- AddMetaData(sd.nonneoplastic, metadata = allen.meta %>% column_to_rownames(var = "cell_barcode"))
table(sd.nonneoplastic$allen_celltype)
sd.nonneoplastic <- AddMetaData(sd.nonneoplastic, metadata = allen.celltypes %>% column_to_rownames(var = "cell_id"))

## Use Allen Celltypes to Set names to the SNN clusters (SNN_clusters_nonNeoplastic)
my.tab <- table(sd.nonneoplastic$allen_celltype, sd.nonneoplastic$SNN_clusters_nonNeoplastic) # 0 NAs
write.csv(as.data.frame.matrix(my.tab) , row.names=T, file = "./results/05-cell-annotation/05-2-snn-cluster-non-neoplastic_vs_AllenMapMyCells.csv")
# write one tab with neurons kapt in original allen classes, allen cluster_name
my.tab <- table(sd.nonneoplastic$cluster_name, sd.nonneoplastic$SNN_clusters_nonNeoplastic) # 0 NAs
write.csv(as.data.frame.matrix(my.tab) , row.names=T, file = "./results/05-cell-annotation/05-2-snn-clusters-non-neoplastic_vs_AllenMapMyCellsCluster")
# tabulate disease state,



