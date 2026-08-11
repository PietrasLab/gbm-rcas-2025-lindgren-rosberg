
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


## ---- filter cell QC ----
cli::cli_alert_info("Pre-filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")
seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )
seurat.object
# filter to th 64k cells
# An object of class Seurat
# 19073 features across 64804 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features
cli::cli_alert_info("Filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")


## --- Metadata: RCAS transgene   ----
rcas.module.df <- readr::read_csv(file = "./results/03-qc/03-6-rcas-module-score.csv" )
stopifnot(identical(rcas.module.df$cell_barcode, colnames(seurat.object)))
# Create an rcas-call column based on rcas module score and rcas counts
rcas.module.df$rcas_call <- dplyr::case_when(
  rcas.module.df$rcas_module_bin == "rcas_module_neg" &
    rcas.module.df$rcas_counts_bin == "rcas_pos" ~ "rcas_pos",

  rcas.module.df$rcas_module_bin == "rcas_module_neg" ~ "rcas_neg",
  rcas.module.df$rcas_module_bin == "rcas_module_pos" ~ "rcas_pos",
  rcas.module.df$rcas_module_bin == "rcas_module_amb" ~ "rcas_ambiguous"
)
table(rcas.module.df$rcas_call, useNA = "ifany")

seurat.object@meta.data <- seurat.object@meta.data %>%
  dplyr::left_join(rcas.module.df, by = c("cell_barcode")) %>%
  column_to_rownames(var = "cell_barcode")

colnames(seurat.object@meta.data)
seurat.object@meta.data$cell_barcode <- colnames(seurat.object)
colnames(seurat.object)



### ---- Metadata: SNN clustering ----
# SNN clustreing of all 64804 cells.
# These clusters, together with RCAS transgene expression annotations form the basis for cell type annotation.
snn.df <- readr::read_csv(file = "./results/04-reductions-snn/04-2-SNN-clusters-all.csv")
snn.df
str(snn.df)
stopifnot(identical(snn.df$cell_barcode, colnames(seurat.object)))
seurat.object@meta.data <- seurat.object@meta.data %>%
  dplyr::left_join(snn.df, by = c("cell_barcode")) %>%
  column_to_rownames(var = "cell_barcode")
colnames(seurat.object@meta.data)
seurat.object@meta.data$cell_barcode <- colnames(seurat.object)
colnames(seurat.object)
table(seurat.object@meta.data$SNN_clusters_all)


### ---- ADD UMAP ----
# add main umap for 64k cells
umap_df <- readr::read_csv("./results/04-reductions-snn/04-1-umap-m1v2-n16-Pds10k-C1-G1.csv")
umap_df
umap_df <- umap_df %>% tibble::column_to_rownames(var = "cell_barcode")
seurat.object@reductions[["umap_all"]] <- Seurat::CreateDimReducObject(as.matrix(umap_df))


## ---- define gene filter for SNN and UMAPS ----
G1 <- .CreateSeuratFilter(
  min.cells.expressed = 25,
  features.remove = c(
    gene.list$mitochondria,
    gene.list$haemoglobin,
    gene.list$rcas_flex_probes,
    gene.list$y_chr,
    "Cst3")
)

# ---- Cell Annotation ----
# Multiple steps of cell annotations
mdata <- seurat.object@meta.data
str(mdata)


## ---- Step 1: UMAP SNN clusters (full data) - all 64k cells----
#  script 04-2-snn-clustering-full-data.R
#  metadata annotation seurat.object@meta.data$SNN_clusters_all
#  26 SNN clusters
table(seurat.object@meta.data$SNN_clusters_all)

# plot SNN clusters on UMAP
p <- Seurat::DimPlot(
  seurat.object, reduction = "umap_all",
  group.by = "SNN_clusters_all", label = TRUE, label.size = 5) +
  ggtitle("SNN clusters all cells") +
  theme_minimal() +
  theme(legend.position = "none")
p
ggsave(filename = file.path(results.dir, "05-1-umap-snn-clusters-all.png"), plot = p, width = 6, height = 5, dpi = 300)




## ---- Step 2: Define SNN_Clusters_broad: Broad SuperTypes - snn clusters into Tumor, Non-Tumor and Mixed----
# Based on RCAS expression proxy - let this define neoplastic cells
# Some clusters include a mixture of Tumor and Non-Tumor cells
library(dplyr)
library(ggplot2)
table(mdata$rcas_counts_bin)
table(mdata$rcas_module_bin)


p1 <- .barplot_stacked(
  my.pal = .color_pal[["rcas_counts_bin"]],
  scaled.y = T,
  mdata, group.var = "SNN_clusters_all",
  color.var =  "rcas_counts_bin") +
  ggtitle("RCAS counts by SNN clusters (all cells)") +
  xlab("SNN clusters all cells") +
  ylab("Percentage of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
p2 <- .barplot_stacked(
  my.pal = .color_pal[["rcas_module_bin"]],
  scaled.y = T,
  mdata, group.var = "SNN_clusters_all",
  color.var =  "rcas_module_bin") +
  ggtitle("RCAS counts by SNN clusters (all cells)") +
  xlab("SNN clusters all cells") +
  ylab("Percentage of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(patchwork)
# ggsave(p1 + p2, filename = file.path(results.dir, "05-1-1RCAS-by-SNN-clusters-all-cells.png"), width = 10, height = 5, dpi = 300)

#
p3 <- .barplot_stacked(
  my.pal = .color_pal[["sample_type"]],
  scaled.y = T,
  mdata, group.var = "SNN_clusters_all",
  color.var =  "sample_type") +
  ggtitle("Sample_Type by SNN clusters (all cells)") +
  xlab("SNN clusters all cells") +
  ylab("Percentage of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(patchwork)
ggsave(p1 + p2 + p3, filename = file.path(results.dir, "05-1-barplots-snn-clusters-all.pdf"), width = 18, height = 5, dpi = 300)


# use the rcas_counts_bin and. rcas_module_bin to type the SNN clusters roughly into
# Neoplastic, Mixed and Non-neoplastic
mdata <- mdata %>%
  mutate(
    SNN_Clusters_broad = case_when(
      SNN_clusters_all %in% c("01","03","05","10","11","12") ~ "Neoplastic Clusters",
      SNN_clusters_all %in% c("06","15","16","21") ~ "Mixed Clusters",
      TRUE ~ "Non-neoplastic Clusters"
    )
  )
# seurat.object$SNN_Clusters_broad <- mdata$SNN_Clusters_broad
seurat.object <- AddMetaData(seurat.object, metadata = tibble(mdata) %>% dplyr::select(cell_barcode,SNN_Clusters_broad) %>% column_to_rownames(var = "cell_barcode"))
table(seurat.object$SNN_Clusters_broad)
# save SNN clusters Broad assignment
write.csv(mdata %>% dplyr::select(cell_barcode,SNN_Clusters_broad), file = file.path(results.dir, "05-1-snn-clusters-broad.csv"), row.names = FALSE, quote = FALSE)


# check how rcas_module_score distributes in the SNN clusters
p1 <- peRcebe_violin(mdata, score.var = "rcas_module_score", group.var = "SNN_clusters_all")
# check how rcas_module_score distributes in the SNN_Clusters_broad
p2 <- peRcebe_violin(mdata, score.var = "rcas_module_score", group.var = "SNN_Clusters_broad",  my.pal = .color_pal$SNN_Clusters_broad)
ggsave(p1 + p2, filename = file.path(results.dir, "05-1-1RCAS-module-score-by-SNN-clusters-broad.png"), width = 10, height = 5, dpi = 300)

# plot the broad SNN clusters on UMAP
p1 <- Seurat::DimPlot(cols = .color_pal$rcas_counts_bin,
  seurat.object, reduction = "umap_all",
  group.by = "rcas_counts_bin", label = TRUE, label.size = 5) +
  ggtitle("RCAS counts bin") +
  theme_minimal() +
  theme(legend.position = "none")
p1
p2 <- Seurat::DimPlot(cols = .color_pal$rcas_module_bin,
  seurat.object, reduction = "umap_all",
  group.by = "rcas_module_bin", label = TRUE, label.size = 5) +
  ggtitle("RCAS module bin") +
  theme_minimal() +
  theme(legend.position = "none")
p3<- Seurat::DimPlot(
  seurat.object, reduction = "umap_all",
  group.by = "SNN_clusters_all", label = TRUE, label.size = 5) +
  ggtitle("SNN clusters all cells") +
  theme_minimal() +
  theme(legend.position = "none")
p4 <- Seurat::DimPlot(cols = .color_pal$SNN_Clusters_broad,
  seurat.object, reduction = "umap_all",
  group.by = "SNN_Clusters_broad", label = TRUE, label.size = 5) +
  ggtitle("SNN clusters all cells") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(p1 + p2 + p3 + p4, filename = file.path(results.dir, "05-1-1UMAPS-SNN-clusters-all.png"), width = 10, height = 10, dpi = 300)



## ---- Step 3: SNN Non-Neoplastic cells  ----

# non-neoplastic cells are rcas negative/ambient cells from the non-neoplastic clusters and rcas negative from the
# Selet non-neoplastic (non RCAS positive cells) from the non-neoplastic and mixed clusters
# perform a separate SNN clustering of these cells, then use publicly available databases to annotate these clusters
table(seurat.object$rcas_module_bin, seurat.object$rcas_counts_bin)

# check that nothing is left unassigned
table(seurat.object$rcas_call, useNA = "always")
sd.nonneoplastic <- subset(seurat.object,
  subset = SNN_Clusters_broad != "Neoplastic Clusters" & rcas_call == "rcas_neg")
sd.nonneoplastic

# # 19073 features across 30345 samples within 1 assay

## XX NonTumor UMAP (v4). umap_m1v2_n16_SubSet_Tx05_nonTumor_v4_C1_G1
## XX ## RCAS neg v4: Tx05_nonTumor_v4
## XX 06-SNN-Clustering-n16-C1-m1v2.R
## XX This i) to reduce the level of ambiguous, perhaps tumor OPC cells but ii) mostly to make the selection of nonTumor not be dependent on Tx05_s1c subclusterings that was
## XX used in nonTumor v2. These seemed not logical and problematic when selecting the Tumo cell branch of clusters (leading to double classifications.)
## XX SNN_p15_res05_Tx05_nonTumor_v4_C1G1
## XX I: Select only cells <0 of the RCAS top25 signature `rcas_sig_top25_3bin == "rcas_sig_top25_neg"` generated from `rcas_3bin_pos_vs_neg_top25pos_andTg < -0.1`)
## XX II: for stringency, remove all cells that are `rcas_pos` in `sd.subset@meta.data[["rcas_both_3bin"]]`, these are few but there are some 286.
## XX II: Do not include cells in strict Tumor cluster, only keep clusters that are inlcude non tumor cells. Most of these have a significant proportion of healthy cells, or a significant proportion of rcas neg cells (these are the opc-6, ocp-16 and neuron-21 clusters, as well as the Astro-16 cluster  )
run.name <- "SNN_clusters_nonNeoplastic"
dims_pca <- 30
nfeatures <- 2000
sd.nonneoplastic <- sd.nonneoplastic %>%
  .seuratFilterFoo(G1) %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nfeatures) %>%
  Seurat::ScaleData(verbose = T) %>%
  Seurat::RunPCA(verbose=TRUE, npcs = dims_pca)

# Run SNN graph
sd.nonneoplastic <- FindNeighbors(
  sd.nonneoplastic,
  dims = 1:dims_pca,
  k.param = 30,
  # prune.SNN = 1/15,
  graph.name = c("kNN_k30","SNN_p15")
)
graph.name <- "SNN_p15"
names(sd.nonneoplastic@graphs)
res <- 0.5
sd.nonneoplastic <- FindClusters(
  sd.nonneoplastic,
  graph.name = graph.name,
  resolution = res,
  algorithm = 1
)
graph.cluster <- paste0(graph.name,"_res.",res)
sd.nonneoplastic@meta.data[,run.name] <- sapply(
  as.numeric(
    sd.nonneoplastic@meta.data[,graph.cluster]), function(x){
      if (x < 10) {
        return(paste0("0", x))
      } else {
        return(as.character(x))
      }
    })
table(sd.nonneoplastic@meta.data[,run.name])
# 01   02   03   04   05   06   07   08   09   10   11   12   13   14   15   16   17   18   19   20
# 4987 2833 2793 2408 1989 1807 1565 1409 1384  977  973  940  879  810  746  670  532  435  379  366
# 21   22   23   24   25   26   27
# 326  286  272  229  136  122   92

### ---- Save nonNeoplastic SNN cluster info to file ----
snn.df <- sd.nonneoplastic@meta.data[,c("cell_barcode",run.name)]
str(snn.df)
# save csv with cluster assignments
write.csv(snn.df, file = file.path(results.dir, "05-1-snn-clusters-non-neoplastic.csv"), row.names = FALSE, quote = FALSE)

### ---- UMAP NonNeoplastic ----
reduction.name <- "umap_05_1_NonNeoplastic_filtered"
cli::cli_alert_info("Running PCA and UMAP with {dims_pca} dimensions and reduction name: {reduction.name}")
sd.nonneoplastic <- sd.nonneoplastic %>%
  Seurat::RunUMAP(reduction.name = reduction.name,
    dims = 1:dims_pca, verbose = T)

# save UMAP as csv
umap.df <- as.data.frame(sd.nonneoplastic@reductions[[reduction.name]]@cell.embeddings)
colnames(umap.df) <- c("UMAP_1", "UMAP_2")
umap.df$cell_barcode <- rownames(umap.df)
umap.df <- umap.df[,c("cell_barcode", "UMAP_1", "UMAP_2")]
# umap.df$cell_barcode <- rownames(umap.df) %>%
readr::write_csv(umap.df, "./results/05-cell-annotation/05-1-umap-nonneoplastic-filtered.csv")



### ---- Plot Barplots and UMAP for  NonNeoplastic ----
library(dplyr)
library(ggplot2)
table(sd.nonneoplastic$SNN_clusters_nonNeoplastic)
table(sd.nonneoplastic$rcas_call)
mdata <- sd.nonneoplastic@meta.data


p1 <- .barplot_stacked(
  my.pal = .color_pal[["rcas_counts_bin"]],
  scaled.y = T,
  mdata,
  group.var = "SNN_clusters_nonNeoplastic",
  color.var =  "rcas_counts_bin") +
  ggtitle("RCAS counts by SNN clusters (nonNeoplastic only)") +
  xlab("SNN_clusters_nonNeoplastic") +
  ylab("Percentage of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
p2 <- .barplot_stacked(
  my.pal = .color_pal[["rcas_module_bin"]],
  scaled.y = T,
  mdata, group.var = "SNN_clusters_nonNeoplastic",
  color.var =  "rcas_module_bin") +
  ggtitle("RCAS counts by SNN clusters (nonNeoplastic only)") +
  xlab("SNN_clusters_nonNeoplastic") +
  ylab("Percentage of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(patchwork)
# ggsave(p1 + p2, filename = file.path(results.dir, "05-1-1RCAS-by-SNN-clusters-all-cells.png"), width = 10, height = 5, dpi = 300)

#
p3 <- .barplot_stacked(
  my.pal = .color_pal[["sample_type"]],
  scaled.y = T,
  mdata, group.var = "SNN_clusters_nonNeoplastic",
  color.var =  "sample_type") +
  ggtitle("Sample_Type by SNN clusters (nonNeoplastic only)") +
  xlab("SNN_clusters_nonNeoplastic") +
  ylab("Percentage of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3
ggsave(p1 + p2 + p3, filename = file.path(results.dir, "05-1-barolots-snn-clusters-non-neoplastic.pdf"), width = 18, height = 5, dpi = 300)

# plot the broad SNN clusters on UMAP
p1 <- Seurat::DimPlot(cols = .color_pal$rcas_counts_bin,
  sd.nonneoplastic, reduction = "umap_05_1_NonNeoplastic_filtered",
  group.by = "rcas_counts_bin", label = TRUE, label.size = 5) +
  ggtitle("RCAS counts bin") +
  theme_minimal() +
  theme(legend.position = "none")
p1
p2 <- Seurat::DimPlot(cols = .color_pal$rcas_module_bin,
  sd.nonneoplastic, reduction = "umap_05_1_NonNeoplastic_filtered",
  group.by = "rcas_module_bin", label = TRUE, label.size = 5) +
  ggtitle("RCAS module bin") +
  theme_minimal() +
  theme(legend.position = "none")
p3<- Seurat::DimPlot(
  sd.nonneoplastic, reduction = "umap_05_1_NonNeoplastic_filtered",
  group.by = "SNN_clusters_nonNeoplastic", label = TRUE, label.size = 5) +
  ggtitle("SNN clusters non-neoplastic cells") +
  theme_minimal() +
  theme(legend.position = "none")
p4 <- Seurat::DimPlot(
  cols = .color_pal$SNN_Clusters_broad,
  sd.nonneoplastic, reduction = "umap_05_1_NonNeoplastic_filtered",
  group.by = "SNN_Clusters_broad", label = TRUE, label.size = 5) +
  ggtitle("SNN clusters non-neoplastic") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(p1 + p2 + p3 + p4, filename = file.path(results.dir, "05-1-umap-snn-clusters-non-neoplastic.pdf"), width = 10, height = 10, dpi = 300)

