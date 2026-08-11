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



## ---- Normalize data ----
sd.nonneoplastic <- sd.nonneoplastic %>%
  NormalizeData()

table(sd.nonneoplastic$SNN_clusters_nonNeoplastic)
sd.nonneoplastic$SNN_clusters_nonNeoplastic <- paste0("cl", as.character(sd.nonneoplastic$SNN_clusters_nonNeoplastic))


# ---- FGSEA PER clusterset ----
# See if any of the recursive subclusters identify as a known cell type
# subclusters are compared to other cells (As a group) using FindAllMarkers
cli::cli_h2("Starting FGSEA")

# Load cell type list used for fgsea
# Cell signatures curated from multiple databases

# DATABASES
# CM2db = CellMarker 2.0 database.  Congxue Hu, Tengyue Li, Yingqi Xu, Xinxin Zhang, Feng Li, Jing Bai, Jing Chen, Wenqi Jiang, Kaiyue Yang, Qi Ou, Xia Li, Peng Wang, Yunpeng Zhang, CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D870–D876, https://doi.org/10.1093/nar/gkac947
#   http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_annotation.jsp.
# PGdb = PanglaoDB = PanglaoDB database. Oscar Franzén, Li-Ming Gan, Johan L M Björkegren, PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data, Database, Volume 2019, 2019, baz046, doi:10.1093/database/baz046
#   https://panglaodb.se/
# Oconnor  = Cell cycle signatures = O'Connor SA, Feldman HM, Arora S, Hoellerbauer P, Toledo CM, Corrin P, Carter L, Kufeld M, Bolouri H, Basom R, Delrow J, McFaline-Figueroa JL, Trapnell C, Pollard SM, Patel A, Paddison PJ, Plaisier CL. Neural G0: a quiescent-like state found in neuroepithelial-derived cells and glioma. Mol Syst Biol. 2021 Jun;17(6):e9522. doi: 10.15252/msb.20209522. PMID: 34101353; PMCID: PMC8186478.
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC8186478/#sec52
# Azim = HuBMAP Azimuth Motor Cortex clusters . The data used in this reference was generated as part of the effort by Yao, Liu, Xie, Fischer, et al, bioRxiv 2020 to generate an integrated transcriptomic and epigenomic atlas of the primary motor cortex in mouse. To generate the reference used in Azimuth, we applied a reference-based integration strategy to harmonize the 10x v3 snRNA-seq data across 12 individual mice. Annotations are provided at the level of class (e.g. GABAergic), subclass (e.g. L2/3 IT), cluster (e.g. L2/3 IT_3), and cross-species cluster (e.g. Sst_6)
#   https://azimuth.hubmapconsortium.org/references/#Human%20-%20Motor%20Cortex
# GBmap = Meta-analysis. Human Glioblastoma and brain cell types. Ruiz-Moreno C, Salas SM, Samuelsson E, Minaeva M, Ibarra I, Grillo M, Brandner S, Roy A, Forsberg-Nilsson K, Kranendonk MEG, Theis FJ, Nilsson M, Stunnenberg HG. Charting the single-cell and spatial landscape of IDH-wild-type glioblastoma with GBmap. Neuro Oncol. 2025 Oct 14;27(9):2281-2295. doi: 10.1093/neuonc/noaf113. PMID: 40312969; PMCID: PMC12526130.
#     https://pmc.ncbi.nlm.nih.gov/articles/instance/12526130/bin/noaf113_suppl_supplementary_tables_s1-s5.xlsx
celltype_list <-readRDS(file = "./references/genesets/fgsea-celltypes.Rds")

names(celltype_list)

fgsea.dir <- file.path(results.dir, "05-3-fgsea-results")
dir.create(fgsea.dir, showWarnings = FALSE, recursive = TRUE)

.wrapper_fgsea(
  seurat.object = sd.nonneoplastic, # pre-filtered data
  # cluster.label =  clusterset, # vector
  cluster.label =  "SNN_clusters_nonNeoplastic", # vector
  gene.sets.list = celltype_list,
  results.path = fgsea.dir,
  analysis.name = glue("celltype-signatures")
)

fgsea.plot.file <- file.path(results.dir,  "05-3-dotplots-fgsea.pdf")
pdf(fgsea.plot.file, width = 20, height = 20)
.wrapper_fgsea_plot(
  fgsea.dir = fgsea.dir,
  annotation = "SNN_clusters_nonNeoplastic",
  seurat.object = sd.nonneoplastic,
  analysis.name = glue("celltype-signatures")
  )
dev.off()



# ---- MCA PER clusterset ----
cli::cli_h2("Starting MCA")
mca.dir <- file.path(results.dir, "05-3-mca-results")
dir.create(mca.dir, showWarnings = FALSE, recursive = TRUE)

.wrapper_mca(results.dir = mca.dir,
  clustersets = "SNN_clusters_nonNeoplastic",
  sd = sd.nonneoplastic
)
.wrapper_mca_plot(
  results.dir = mca.dir,
  sd.ref = sd.nonneoplastic,
  clustersets =  c("SNN_clusters_nonNeoplastic")
)
