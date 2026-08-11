# compile Filtered Seurat Object with umap and annotations

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

# seurat.object <- AddMetaData(
#   seurat.object, metadata = rcas.module.df
#   %>% column_to_rownames(var = "cell_barcode"))
colnames(seurat.object@meta.data)
table(rcas.module.df$rcas_module_bin)
# for call, override when module says negative but call says positve
rcas.module.df$rcas_call <- dplyr::case_when(
  rcas.module.df$rcas_module_bin == "rcas_module_neg" & rcas.module.df$rcas_counts_bin == "rcas_pos" ~ "rcas_pos",
  rcas.module.df$rcas_module_bin == "rcas_module_neg" ~ "rcas_neg", rcas.module.df$rcas_module_bin == "rcas_module_pos" ~ "rcas_pos",
  rcas.module.df$rcas_module_bin == "rcas_module_amb" ~ "rcas_ambiguous"
)
seurat.object <- AddMetaData(seurat.object, metadata = rcas.module.df  %>% column_to_rownames(var = "cell_barcode"))
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

umap_df <- readr::read_csv("./results/05-cell-annotation/05-1-umap-nonneoplastic-filtered.csv")
seurat.object[["umap_NonNeoplastic"]] <- Seurat::CreateDimReducObject(
  embeddings = as.matrix(tibble::column_to_rownames(umap_df, "cell_barcode")),
  # key = "umapNN_",
  assay = DefaultAssay(seurat.object)
)



# ---- Cell type Annotations ----

# Add as new metadata column
# Data is fist split in broad cluster types (04-2-SNN-clusters-all), tumor cells, mixed and non-tumor cells.
# Then proceed by recyrsive clustering of non-neoplastic cells in non-tumor and mixed clusters (05-3-SNN-clusters-nonNeoplastic)
# Use these clusters to annotate non-neoplastic cells, using Mappings from Allen MapMyCells  (05-2), fgsea and MCA (05-3)
# Tumor/Neoplastic cells are subgrouped using the cluster assignments from 04-2-SNN-clusters-all
# Ambiguous cells are rcas positive or ambiguous but that cluster within non-neoplastic clusters,
# or cells from healthy tissue (and rcas negative) that cluster within neoplastic clusters.

mdata <- seurat.object@meta.data
table(mdata$rcas_call)
table(mdata$SNN_Clusters_broad)
table(mdata$SNN_clusters_nonNeoplastic)
table(mdata$SNN_clusters_all)

## ---- Level 1  -----

# Level 1 is:
# 1. Non.Neoplastic: Cells that were included in 05-3-SNN-clusters-nonNeoplastic.csv
#     These are Non-neoplastic and Mixed Clusters (from SNN of full data set)
#     30345 Cells
#     These cells are classificed as RCAS-negative (not positive or negative) AND are in NonNeoplastic or 4 Mixed SNN clusters (broad)
#     Non-neoplastic cells are re-clustered to get a more fine tuned SNN cluster map
#     This SNN lays ground for Level2 and Level3 mapping of non-neoplastic cell types
# 2. Neoplastic;
#     33727 cells. A pretty non-conservative approach
#     These cells are all cells within the 10 Neoplastic SNN clusters (broad)
#     Also, ambiguous and positive cells that were in the 4 Mixed SNN clusters (broad) are added
#     Cells from healthy tissue are not included and set to Ambiguous (RARE!, n=72)
# 3. Ambiguous:
#     732 cells
#     Cells that were classified as Ambiguous or Positive in the RCAS-call AND were in NonNeoplastic clusters (n=660)
#     OR Cells from healthy tissue are not included and set to Ambiguous (RARE!, n=72)



# create Level_1 (set all to neoplastic)
mdata$Level_1 <- "Neoplastic"
# Set all Non-neoplastic (has value in mdata$SNN_clusters_nonNeoplastic) to "Non-Neoplastic"
mdata$Level_1[!is.na(mdata$SNN_clusters_nonNeoplastic)] <- "Non-Neoplastic"
table(mdata$Level_1)

# Now set Ambiguous cells. Cells that are rcas positive or ambiguous but that cluster within non-neoplastic clusters
# A conservative approach - these may be technical artefacts or rare true positives
mdata$Level_1[
  mdata$SNN_Clusters_broad == "Non-neoplastic Clusters" &
    mdata$rcas_call != "rcas_neg"
] <- "Ambiguous" # 660 cells that are rcas positive but that cluster within apparent non-tnumor clusters. Set as Ambiguous - may be either or
table(mdata$Level_1)

mdata$Level_1[
  mdata$Level_1 == "Neoplastic" &
    mdata$sample_type == "Healthy"
] <- "Ambiguous" # 72 cells from healhy tissue but that cluster with neoplastic cells
table(mdata$Level_1)
# Ambiguous     Neoplastic Non-Neoplastic
# 732          33727          30345



## ---- Level 2 and Level 3  ----
# split the non-neoplastic cells to Glial, Immune/Microglia,and Vascular
# The non-neoplastic cells are clustered separately
# 05-3-SNN-clusters-nonNeoplastic.
# Cluster annotations are performed by MCA, fgsea amd Allen MapMyCells tool
# Also, look at proportions of tumor vs healthy brain cells in clusters (e.g. 5-2-RCAS-by-SNN-clusters-all-cells.pdf)

# Cluster cell namings map for 05-3-SNN-clusters-nonNeoplastic (27 clusters)
cluster_map <- c(
  # Glial lineage
  "02" = "02_Astro_TE",         # Healthy-dominated. Astrocytes, telencephalon (CTX/HIP etc.).
  "07" = "07_Astro_NT",         # Healthy-dominated. Astrocytes, non-telencephalic.
  "15" = "15_Astro_TA",         # Tumor-associated tissue–enriched astrocytes (TA). Cells clusters together with tumor cells.

  "24" = "24_Ependymal",        # Healthy-dominated. Ependymal cells.

  # Neuronal
  "11" = "11_Inh_interneuron",  # Healthy-exclusive inhibitory interneurons (GABA-dominant, ~5% glut); OB/STR/CTX origin.
  "18" = "18_Neuron_mixed",     # Mixed neurons (Glut + GABA + Dopa) with a few neoplastic cells intermixed.

  # OPC / Oligo
  "03" = "03_Oligodendrocyte",  # Healthy-dominated. Mature myelinating oligodendrocytes.
  "20" = "20_COP_TA",           # COP-like (PDGFRA–), OPC/Oligo; TA-enriched; co-clusters with neoplastic cells.
  "14" = "14_OPC",              # Healthy-dominated. OPCs (mostly healthy tissue).
  "12" = "12_OPC_TA",           # Tumor-associated tissue–exclusive OPC-like; co-clusters with neoplastic cells.

  # Immune – lymphoid & pDCs
  "13" = "13_NKT",              # Tumor-associated tissue–dominated T and NK cells.
  "25" = "25_B_pDC",            # Mix of B cells and pDCs (split at Level 3: pDC → Dendritic; B → NKTB).

  # Immune – myeloid & microglia
  "01" = "01_Microglia",        # Healthy-exclusive resident microglia.
  "09" = "09_Microglia",        # Healthy-dominated resident microglia.
  "08" = "08_Microglia_TAM",    # Tumor-dominated resident-like TAM (microglia phenotype).
  "06" = "06_Microglia_TAM",    # Tumor-exclusive bone marrow–derived TAM (macrophage phenotype).
  "16" = "16_Macrophage_BMD_TAM", # BMD-derived TAM.
  "04" = "04_Macrophage_BMD_TAM", # Tumor-exclusive BMD-derived TAM.
  "17" = "17_Dendritic",        # Dendritic cells (cDC-like; non-pDC).
  "19" = "19_BAM",              # Healthy-dominated border-associated macrophages (BAM).
  "26" = "26_Neutrophil",      # Tumor-associated tissue–dominated neutrophils.

  # Vascular & stromal
  "10" = "10_Choroid",          # Healthy-exclusive choroid plexus epithelium.
  "05" = "05_Endothelial",      # Dominated by Healthy Brain endothelial (BBB) cells capillary/venous ECs, Also, a subset of cells with mural/pericyte/smc markers and Myofibroblasts .
  "23" = "23_Endothelial",      # Mostly Tumor-associated endothelial cells, mixed cluster, Name Mural at this level
  "21" = "21_Mural",            # Mural cells: Pericyte and SMCs, pericytes / smooth muscle cells (healthy-dominated).
  "22" = "22_Fibroblast",       # Fibroblasts (present in healthy and tumor tissue).
  "27" = "27_Fibroblast_meningeal" # Healthy-dominated, meningeal fibroblasts.
)

mdata$SNN_clusters_nonNeoplastic_celltypes_l3 <- plyr::mapvalues(
  x   = seurat.object$SNN_clusters_nonNeoplastic,
  from= names(cluster_map),
  to  = base::unname(cluster_map)
)
sum(table(mdata$SNN_clusters_nonNeoplastic_celltypes))  # 30345
table(mdata$SNN_clusters_nonNeoplastic_celltypes_l3, useNA = "ifany")
seurat.object <- AddMetaData(seurat.object, metadata = tibble(mdata)
  %>% dplyr::select(cell_barcode, SNN_clusters_nonNeoplastic_celltypes_l3)  %>% column_to_rownames(var = "cell_barcode"))

annot <- "SNN_clusters_nonNeoplastic_celltypes_l3"
reduction.name <- "umap_NonNeoplastic"
Idents(seurat.object) <- annot

sd.sub <- seurat.object %>% subset(., SNN_clusters_nonNeoplastic != "NA")

table(seurat.object[["SNN_clusters_nonNeoplastic_celltypes_l3"]])
p1 <- Seurat::DimPlot(
  sd.sub,
  # cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot,
  label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p1
annot <- "SNN_clusters_nonNeoplastic"
p2 <- Seurat::DimPlot(
  sd.sub,
  #cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p2
reduction.name <- "umap_AllSamples"
annot <- "SNN_clusters_nonNeoplastic"
p3 <- Seurat::DimPlot(
  sd.sub,
  #cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p3
reduction.name <- "umap_AllSamples"
annot <- "SNN_clusters_nonNeoplastic_celltypes_l3"
p4 <- Seurat::DimPlot(
  sd.sub,
  #cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p4
# ggsave(p1 + p2 + p3 + p4,
#   filename = "./results/05-cell-annotation/05-4-umap-snn-nonneoplastic-clusters.png",
#   width = 20, height = 20, units = "in", dpi = 300)



### ---- Level 2 buckets. broad ----
# Buckets: Glial, Neuronal, Immune, Vascular/Stromal
cluster_map <- c(
  # Neural cells: Glial cells and cells of the nervous system (Astro/OLG/Neuronal/Ependymal/OECs etc)
  "02" = "Neural", "07" = "Neural", "15" = "Neural", # Astrocyte subsets, glial lineage
  "24" = "Neural",               # Ependymal, glial  .
  "03" = "Neural", "20" = "Neural", "14" = "Neural", "12" = "Neural", # opc, olig, glial
  "11" = "Neural", "18" = "Neural", # Neuronal and possibly uncharacterized glial cells, mixed neuorons

  # Immune (lymphoid + myeloid + microglia)
  "13" = "Immune", "25" = "Immune", # Lymphoid/pDC
  "01" = "Immune", "09" = "Immune", "08" = "Immune", "06" = "Immune", # microglia
  "16" = "Immune", "04" = "Immune", "17" = "Immune", "19" = "Immune",  #Macrophage/dendritic
  "26" = "Immune", # neutrophils

  # Vascular / Stromal
  "10" = "Vascular/Stromal",    # Choroid plexus epithelium.
  "05" = "Vascular/Stromal", "21" = "Vascular/Stromal",
  "22" = "Vascular/Stromal", "23" = "Vascular/Stromal",
  "27" = "Vascular/Stromal"
)

mdata$Level_2 <- plyr::mapvalues(
  x   = seurat.object$SNN_clusters_nonNeoplastic,
  from= names(cluster_map),
  to  = base::unname(cluster_map)
)
sum(table(mdata$Level_2))  # 30345
table(mdata$Level_2, useNA = "ifany")

## Add neoplastic/ambiguous
mdata$Level_2 <- ifelse(
  test = mdata$Level_1 %in% c("Neoplastic", "Ambiguous"),
  yes  = mdata$Level_1,
  no   = mdata$Level_2
)
table(mdata$Level_2, useNA = "ifany")
# Ambiguous   Glial/Neuronal           Immune       Neoplastic Vascular/Stromal
# 732            11690            14713            33727             3942


### ---- Level 3 buckets - semi/detailed  ----
cluster_map <- c(
  # Glial lineage
  "02" = "Astrocyte",
  "07" = "Astrocyte",
  "15" = "Astrocyte",

  "24" = "Ependymal",
  "03" = "OPC/COP/OLG",         #  Mature myelinating oligodendrocytes.
  "20" = "OPC/COP/OLG",         # COP-like, PDGFRA–, TA-enriched.
  "14" = "OPC/COP/OLG",         # healthy OPCs (mostly healthy tissue).
  "12" = "OPC/COP/OLG",         # Tumor-associated tissue–exclusive OPC-like; co-clusters with neoplastic OPC like cells

  # Neuronal
  "11" = "Neural",              # Interneurons (GABA-dominant).
  "18" = "Neural",              # Mixed neurons (Glut/GABA/Dopa).

  # Immune – lymphoid & pDCs
  "13" = "NKTB",
  "25" = "Dendritic",           # pDC grouped here; B cells remain under NKTB at L3 per your note.

  # Immune – myeloid & microglia
  "01" = "Microglia",
  "09" = "Microglia",
  "08" = "Microglia",           # Resident-like TAM counted with microglia at this level.
  "06" = "Microglia",
  "16" = "Macrophage",          # BMD-derived TAM.
  "04" = "Macrophage",
  "19" = "Macrophage",          # BAM.
  "17" = "Dendritic",
  "26" = "Neutrophil",

  # Vascular & stromal
  "10" = "Choroid",              # Choroid plexus epithelium.
  "05" = "Endothelial", # Endothelial cells, a few mural cells are remocved from this cluster further
  "23" = "Endothelial",
  "21" = "Mural",        # Pericyte/SMC cells.
  "22" = "Fibroblast", # Fibroblasts
  "27" = "Fibroblast" # Fibroblasts
)

mdata$Level_3 <- plyr::mapvalues(
  x   = seurat.object$SNN_clusters_nonNeoplastic,
  from= names(cluster_map),
  to  = base::unname(cluster_map)
)
sum(table(mdata$Level_3))  # 30345
## Add neoplastic/ambiguous
mdata$Level_3 <- ifelse(
  test = mdata$Level_1 %in% c("Neoplastic", "Ambiguous"),
  yes  = mdata$Level_1,
  no   = mdata$Level_3
)
table(mdata$Level_3, useNA = "ifany")


### ---- Level 3. Outlier cluster splits ----
# A few SNN_clusters_nonNeoplastic include small distinct outlier clusters
# these are manually annotated here and added to the Level 3 annotations.

# Note! cluster 25, pDCs/B-cells, is split into pDCs and B-Cells. These are added to Level3 as Dendritic and NKTB cells, respectively
cells.df <- readr::read_csv(file = "./results/05-cell-annotation/05-4-cluster-25-B-pDC-split.csv")
cells.df
u <- match(cells.df$cell_barcode, mdata$cell_barcode)
stopifnot(all(!is.na(u)))
mdata$Level_3[u] <- cells.df$pDC_or_Bcell
# For Level 3, set Bcells in NKTB bucket and pDCs in Dendritic bucket
u <- which(mdata$Level_3 == "B cell")
mdata$Level_3[u] <- "NKTB"
u <- which(mdata$Level_3 == "pDC")
mdata$Level_3[u] <- "Dendritic"
table(mdata$Level_3, useNA = "ifany")

# NOTE! Mural (unspecified) cell subcluster of cluster 5. (splits from endothelial)
# Set as Mural
cells.df <- readr::read_csv(file = "./results/05-cell-annotation/05-4-cluster-5-mural-split.csv")
cells.df
u <- match(cells.df$cell_barcode, mdata$cell_barcode)
stopifnot(all(!is.na(u)))
mdata$Level_3[u] <- cells.df$Mural_cell
table(mdata$SNN_clusters_nonNeoplastic_celltypes, useNA = "ifany")
table(mdata$Level_3, useNA = "ifany")


## ---- Level 3Ac - Astrocyte specific Subclusters ----
# In this study we make difference between the different astrocyte subclusters
u <- which(mdata$Level_3 == "Astrocyte")
str(u)
table(mdata$Level_3[u])
table(mdata$SNN_clusters_nonNeoplastic_celltypes_l3[u])
mdata$Level_3AC <- mdata$Level_3
u <- which(mdata$SNN_clusters_nonNeoplastic_celltypes_l3 == "02_Astro_TE")
mdata$Level_3AC[u] <- "Astrocyte TE"
u <- which(mdata$SNN_clusters_nonNeoplastic_celltypes_l3 == "07_Astro_NT")
mdata$Level_3AC[u] <- "Astrocyte NT"
u <- which(mdata$SNN_clusters_nonNeoplastic_celltypes_l3 == "15_Astro_TA") # tumor associated reactive astrocytes
mdata$Level_3AC[u] <- "Astrocyte R"
table(mdata$Level_3AC, useNA = "ifany")
table(mdata$Level_3, useNA = "ifany")

# adress the outlier cells in Astrocyte clusters: Olfactory Ensheathing Cells (OECs)
# Very small-unverified group, 34 cells, set as Neural cells in Level3, i.e.
# group for now with the mixed and less characterized Nuronal (and ependymal) cells.
# Set as Neural cells to mark that these are not strictly neuronal cells.
cells.df <- readr::read_csv("./results/05-cell-annotation/05-4-astrocyte-oec-split.csv")
cells.df
u <- match(cells.df$cell_barcode, mdata$cell_barcode)
stopifnot(all(!is.na(u)))
mdata$Level_3[u] <- "Neural" # Cluster with neuronal and OECs
table(mdata$Level_3, useNA = "ifany")
mdata$Level_3AC[u] <- "Neural" # Cluster with neuronal and OECs
table(mdata$Level_3AC, useNA = "ifany")
table(mdata$Level_3, useNA = "ifany")
# astro 5110

## --- Level 4: Include Neoplastic cells ----
# Level_4: add named neoplastic subclusters on top of the Level_3 labels, from the SNN clustering.
mdata$Level_4 <- ifelse(
  test = mdata$Level_1 %in% c("Neoplastic"),
  yes  = mdata$SNN_clusters_all,
  no   = mdata$Level_3AC
)
table(mdata$Level_4, useNA = "ifany")

cluster_map <- c(
  "01" = "Neopl-Bulk",
  "03" = "Neopl-CC-III", # previously CC-B. G2M-phase
  "05" = "Neopl-CC-II", # previously CC-A. S/G2M phase
  "06" = "Neopl-OPC",
  "10" = "Neopl-ECM",
  "11" = "Neopl-CC-I", # previously "Neopl-Bulk-B". S phase
  "12" = "Neopl-RNA-low",
  "15" = "Neopl-COP",
  "16" = "Neopl-ACR",
  "21" = "Neopl-NC"
)

mdata$Level_4 <- plyr::mapvalues(
  x   = mdata$Level_4,
  from= names(cluster_map),
  to  = base::unname(cluster_map)
)
table(mdata$Level_4, useNA = "ifany")
table(mdata$Level_1, useNA = "ifany")

table(mdata$Level_1, mdata$Level_4, useNA = "ifany")

# ---- Save table with set Level 1-3 annotations ----
colnames(mdata)
df <- mdata[, c(
  "cell_barcode",
  "Level_1","Level_2", "Level_3", "Level_3AC", "Level_4",
   "final_cell_filter","rcas_counts_bin","rcas_module_score","rcas_module_bin","rcas_call",
   "SNN_clusters_all","SNN_Clusters_broad","SNN_clusters_nonNeoplastic"
  )]
head(df)
colnames(df)
readr::write_csv(x = df, file = "./results/05-cell-annotation/05-4-cell-metadata-annotations-levels.csv")

## ---- Add as metadata ----
df.levels <- readr::read_csv(file = "./results/05-cell-annotation/05-4-cell-metadata-annotations-levels.csv")
seurat.object <- AddMetaData(seurat.object, metadata = tibble(df.levels)  %>% column_to_rownames(var = "cell_barcode"))
seurat.object <- .seuratFactorizeMdata(seurat.object)
table(seurat.object$Level_1, useNA = "ifany")
table(seurat.object$Level_2, useNA = "ifany")
table(seurat.object$Level_3, useNA = "ifany")
table(seurat.object$Level_3AC, useNA = "ifany")
table(seurat.object$Level_4, useNA = "ifany")


#
# table(seurat.object$Level_1, useNA = "ifany")
#
# > table(seurat.object$Level_1, useNA = "ifany")
#
# Non-Neoplastic     Neoplastic      Ambiguous
# 30345          33727            732
# > table(seurat.object$Level_2, useNA = "ifany")
#
# Neural           Immune Vascular/Stromal       Neoplastic        Ambiguous
# 11690            14713             3942            33727              732
# > table(seurat.object$Level_3, useNA = "ifany")
#
# Astrocyte   OPC/Oligo      Neural   Ependymal   Microglia  Macrophage   Dendritic  Neutrophil        NKTB
# 5110        4909        1442         229        9587        3457         617         122         930
# Choroid Endothelial       Mural  Fibroblast  Neoplastic   Ambiguous
# 977        1947         640         378       33727         732
# > table(seurat.object$Level_3AC, useNA = "ifany")
#
# Astrocyte NT       Astrocyte TE Astrocyte reactive          OPC/Oligo             Neural
# 1565               2833                746               4909               1408
# Ependymal          Microglia         Macrophage          Dendritic         Neutrophil
# 229               9587               3457                617                122
# NKTB            Choroid        Endothelial              Mural         Fibroblast
# 930                977               1947                640                378
# Neoplastic          Ambiguous
# 33727                732

# ---- Plot UMAPs ----
cli::cli_alert_info("Plotting UMAP")
reduction.name <- "umap_AllCells"
annot <- "sample_type"
p0 <- Seurat::DimPlot(
  seurat.object,
  cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p0
annot <- "Level_1"
p1 <- Seurat::DimPlot(
  seurat.object,
  cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p1


annot <- "Level_2"
p2 <- Seurat::DimPlot(
  seurat.object,
  cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p2


annot <- "Level_3"
p3 <- Seurat::DimPlot(
  seurat.object,
  cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p3

annot <- "Level_3AC"
p4 <- Seurat::DimPlot(
  seurat.object,
  cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p4
annot <- "Level_4"
p5 <- Seurat::DimPlot(
  seurat.object,
  cols = .color_pal[[annot]],
  reduction = reduction.name,
  group.by = annot, label = TRUE, label.size = 3) +
  ggplot2::ggtitle(glue("umap_05_4_{reduction.name}")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")
p5

# ggsave(p0 + p1+ p2 + p3+ p4 + p5,
#   filename = "./results/05-cell-annotation/05-4-umap-Levels-1-to-4.png",
#   width = 28, height = 20, units = "in", dpi = 300)




