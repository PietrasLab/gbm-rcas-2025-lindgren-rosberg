## Obtain a gene signature for Neoplastic cell - expression of transgenes and transgene gene module
## What genes are co-expresed by cells positive for the transgene
## Use the binned proxy defined from count values for rcas transgene cell calling

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


## --- Add transgene binned counts ----
# binned data from raw counts
# use as proxy for transgene expression
# use to obtina transgene/neoplastic gene module
# from 03-4
rcas.bin.df <- readr::read_csv("./results/03-qc/03-5-rcas-transGenes-probeCalls_rcas_both_3bin.csv")
stopifnot(identical(rcas.bin.df$cell_name, colnames(seurat.object)))
seurat.object@meta.data$rcas_bin <- rcas.bin.df$rcas_both_3bin
table(seurat.object@meta.data$rcas_bin)



# ---- Add Metadata----
# filter to th 64k cells
# seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )
cli::cli_alert_info("Filtered Seurat object contains {ncol(seurat.object)} cells and {nrow(seurat.object)} features.")
# An object of class Seurat
# 19073 features across 64804 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features

# Compare rcas_bin - posive vs  negative
cli::cli_alert_info("Subsetting to rcas_bin positive and negative cells.")
Idents(seurat.object) <- "rcas_bin"
table(Idents(seurat.object))
# rcas_neg rcas_ambiguous       rcas_pos
# 28477          12665          23662

# ---- subset to rcas_neg and rcas_pos and find rcas module ----
cli::cli_alert_info("Keeping rcas_neg and rcas_pos cells only.")
keep.idents = c("rcas_neg","rcas_pos")
sd.tmp <- subset(seurat.object, idents = keep.idents )
# 19073 features across 75037 samples within 1 assay
table(Idents(seurat.object))
# rcas_neg rcas_pos
# 37964    37073


# ---- FindAllMarkers - RCAS Module ----
# set ident of seurat.object and sd_full to the input_subset

cli::cli_alert("Normalize input subset data")
sd.tmp <- sd.tmp %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

logfc.threshold = 1
cli::cli_h2("FindAllMarkers")
genelists_findallmarkers <- FindAllMarkers(
  min.diff.pct = 0.1,
  min.pct = 0.1,
  sd.tmp,
  logfc.threshold = logfc.threshold,
  pct.2.cut = 0.1,
  only.pos = TRUE,
  verbose = T
)

# sort on avg_log2FC
deg_df <- genelists_findallmarkers%>%
  dplyr::filter(avg_log2FC > logfc.threshold) %>%
  dplyr::arrange( desc(cluster), desc(avg_log2FC) )
deg_df$gene[1:25]
# save  deg_df as csv
write.csv(deg_df, file = "./results/03-qc/03-6-rcas-findallmarkers.csv", row.names = FALSE)
write.csv(deg_df, file = "./manuscript-tables/03-6-rcas-findallmarkers.csv", row.names = FALSE)
rm(sd.tmp)


# ---- Caclulate rcas-pos module score ----
module.name <- "rcas_module_score"
gene.vec <- deg_df$gene[1:25]

# [1] "Tg-hPDGFB-nHA" "Tg-RFP-2"      "Tg-RFP-1"
# [4] "Hcrtr2"        "Slc38a4"       "Moxd1"
# [7] "Tecta"         "Cdkn2a"        "Mab21l1"
# [10] "Col19a1"       "En1"           "Tafa4"
# [13] "Igf2bp2"       "Raet1d"        "Col11a1"
# [16] "Dll3"          "Fbn2"          "Enpp3"
# [19] "H60b"          "Hmga2"         "Vit"
# [22] "Dlk1"          "Chrm3"         "Pnlip"
# [25] "Depdc1a"


# ---- Filter data ----
# filter to th 64k cells
seurat.object <- subset(seurat.object, subset = final_cell_filter == "keep" )
seurat.object
seurat.object <- Seurat::NormalizeData(seurat.object)

cli::cli_alert_info("Calculating module score for {module.name} with {length(gene.vec)} genes.")
seurat.object <- Seurat::AddModuleScore(
  object = seurat.object,
  features = list(gene.vec),
  name = module.name,
  verbose=T
)
colnames(seurat.object@meta.data)

module.df <- seurat.object@meta.data %>%
  dplyr::rename(!!module.name := !!sym(glue("{module.name}1"))) %>%
  dplyr::select(cell_barcode, rcas_bin, !!module.name)

# Create a binned value based on positve/negative (or undetermined) from the module score
module.df <- module.df %>%
  dplyr::mutate(rcas_module_bin = dplyr::case_when(
    !!sym(module.name) <= -0.1 ~ "rcas_module_neg",
    !!sym(module.name) >=  0.1 ~ "rcas_module_pos",
    TRUE                       ~ "rcas_module_amb"
  )) %>%
  # rename rcas_bin to rcas_counts_bin
  dplyr::rename(rcas_counts_bin = rcas_bin)

table(module.df$rcas_module_bin)
# rcas_module_neg rcas_sig_top25_amb rcas_sig_top25_pos
# 32379               6311              26114
str(module.df)

# save module.df
write_csv(module.df, file = "./results/03-qc/03-6-rcas-module-score.csv" )









