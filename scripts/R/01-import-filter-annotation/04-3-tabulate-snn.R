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



## ---- 05 Step 5: FGSEA and MCAA of Non-Neoplastic subclusters ----
mdata <- seurat.object@meta.data

# Loop all non-neoplastic SNN clusters
colnames(seurat.object@meta.data)
## Use downsampled data, all genes, Non Neoplastic cells only,
sd.nonneoplastic <- subset(seurat.object, subset = rcas_module_bin == "rcas_module_neg")
# 19073 features across 32379 samples within 1 assay
sd.nonneoplastic <- subset(sd.nonneoplastic, subset = rcas_counts_bin != "rcas_pos")
# 19073 features across 32093 samples within 1 assay
table(sd.nonneoplastic@meta.data$SNN_Clusters_broad)
sd.nonneoplastic <- subset(sd.nonneoplastic, subset = SNN_Clusters_broad != "Neoplastic Clusters")
# 19073 features across 30345 samples within 1 assay


table(sd.nonneoplastic$SNN_clusters_nonNeoplastic, useNA = "ifany")
sum(table(sd.nonneoplastic$SNN_clusters_nonNeoplastic, useNA = "ifany")) # 30345
# [1] 30345
active.ident = "SNN_clusters_nonNeoplastic"
Idents(sd.nonneoplastic) <- active.ident
table(Idents(sd.nonneoplastic))


# --- Table of non-neoplastic SNN clusters vs Disease State ----

my.tab <- table(sd.nonneoplastic$disease_state, sd.nonneoplastic$SNN_clusters_nonNeoplastic) # 0 NAs
write.csv(as.data.frame.matrix(my.tab) , row.names=T, file = "./results/04-reductions-snn/04-3-table-snn-non-neopl_disease-state.csv")
# perform chi-sq test and odds ratios to identofy clusters associated with disease state
chisq.test(my.tab) # X-squared = 22994, df = 26, p-value < 2.2e-16
# my.tab = table(disease_state, cluster)
pvals <- sapply(colnames(my.tab), function(cl){
    mat <- matrix(c(my.tab["Healthy", cl],
        sum(my.tab["Healthy", ]) - my.tab["Healthy", cl],
        my.tab["Tumor", cl],
        sum(my.tab["Tumor", ]) - my.tab["Tumor", cl]),
        nrow = 2, byrow = TRUE)
    fisher.test(mat)$p.value
})
pvals_adj <- p.adjust(pvals, method = "BH") # multiple testing correction
results <- data.frame(cluster = names(pvals),
    pval = pvals,
    pval_adj = pvals_adj,
    Healthy = my.tab["Healthy", ],
    Tumor = my.tab["Tumor", ])
results %>% arrange(pval_adj)
# A number of clusters are highly enriched with respect to disease state (Healthy vs tumor)
# Of these clusters 04,06,12,16, 09,13,15, 17 and 23 are particularity enriched for non neoplastic cells from tumor tissue
#  Of these, all but three are Lymphoid or Myeloid populations
# The other thre is an OPC, an Astrocyte, and an Endothelial cluster
# Add this to the cluster names (TA = Tumor Associated)

# Tabulate vs sample_type, identify groups that come from specific samples
my.tab <- table(sd.nonneoplastic$sample_type, sd.nonneoplastic$SNN_clusters_nonNeoplastic) # 0 NAs
write.csv(as.data.frame.matrix(my.tab) , row.names=T, file = "./results/04-reductions-snn/04-3-table-snn-non-neopl_sample-type.csv")
