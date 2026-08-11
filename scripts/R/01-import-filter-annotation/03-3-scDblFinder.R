# ----  scDblFinder ----
# Expected proportion of doublets
#
# The expected proportion of doublets has little impact on the score, but a very strong impact on where the threshold will be placed (the thresholding procedure simultaneously minimizes classification error and departure from the expected doublet rate). It is specified through the dbr parameter and the dbr.sd parameter (the latter specifies the standard deviation of dbr, i.e. the uncertainty in the expected doublet rate). For 10x data, the more cells you capture the higher the chance of creating a doublet, and Chromium documentation indicates a doublet rate of roughly 1% per 1000 cells captures (so with 5000 cells, (0.01*5)*5000 = 250 doublets), and the default expected doublet rate will be set to this value (with a default standard deviation of 0.015). Note however that different protocols may create considerably more doublets, and that this should be updated accordingly. If you are unsure about the doublet rate, set dbr.sd=1 and the thresholding will be entirely based on the misclassification rates.
#
# Note on Multiplet rates. If 6k recovered cells, expect a 4.8% multipleet rate
# https://cdn.10xgenomics.com/image/upload/v1680118519/support-documents/CG000527_Chromium_FixedRNAProfiling_MultiplexedSamples_UserGuide_Rev_D.pdf
#
# Looking at the observed rate of doublets, set a much higher rate - 10%?

## Basic usage
# Given an object sce of class SingleCellExperiment (which does not contain any empty drops, but hasn’t been further filtered), you can launch the doublet detection with: i.u.

# Here input raw full data (but specifying sample id), tunning algorithm on each separate sample
# specify this using the `samples` argument

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


# ---- Load data -----
# 02-1-flex-raw-s16-c92737.qs is a Seurat object containing raw data from 16 samples
# same as seurat_flex_raw_v1.0.rds on Zenodo
# This can be downloaded from Zenodo:
# https://zenodo.org/records/21607378
# Zenodo file: seurat_flex_raw_v1.0.rds (check for latest version)
# sd.merged <- readRDS("./data/processed/seurat/seurat_flex_raw_v1.0.rds")
sd.merged <- qs2::qs_read("./data/processed/seurat/02-1-flex-raw-s16-c92737.qs")
sd.merged
colnames(sd.merged@meta.data)
# An object of class Seurat
# 19073 features across 92737 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 16 layers present: counts.BRAIN_01, counts.BRAIN_04a, counts.BRAIN_04b, counts.BRAIN_05, counts.BRAIN_06, counts.TP_01b, counts.TP_04, counts.TP_05, counts.TP_06, counts.TP_08, counts.TR_01, counts.TR_02, counts.TR_03, counts.TR_04, counts.TR_06, counts.TR_07

source("./scripts/R/00-0-source-functions.R")
pdata <- readRDS("./metadata/00-sample-metadata-flex-s16.Rds")
fdata <- read.delim("./metadata/00-probeset-flex.csv", header = T, sep=",")
str(fdata)
head(fdata)
stopifnot(identical(fdata$external_gene_name, rownames(sd.merged)))
gene.list <- readRDS(file = "./references/genesets/gene-list.Rds")


# ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")


# ---- Add cell barcodes to metadata ----
cell_barcodes <- paste0(sd.merged@meta.data$sample_id, "_", colnames(sd.merged))
# removing trailing _digit AT the END
cell_barcodes <- gsub("_\\d+$", "", cell_barcodes)
sd.merged@meta.data$cell_barcode <- cell_barcodes
# join layers
sd.merged <- SeuratObject::JoinLayers(sd.merged)
colnames(sd.merged) <- sd.merged@meta.data$cell_barcode
# 19073 features across 92737 samples within 1 assay

# Add fitering metadata to object
# The F03 matches cell filter with adaptive cell q95 cutoffs
# cell.filters = list(
#   percent_mito = "percent_mito < 10",
#   percent_haemo = "percent_haemoglobin < 5",
#   nFeature_low = paste0("nFeature_RNA > ", metadata.stats$nGenes_q05[u]),
#   nFeature_high = paste0("nFeature_RNA < ", metadata.stats$nGenes_q95[u]),
#   nCount_low = paste0("nCount_RNA > ",metadata.stats$nUMI_q05[u]),
#   nCount_high = paste0("nCount_RNA < ",metadata.stats$nUMI_q95[u])
# ))
filter.df <- readr::read_csv("./results/03-qc/03-2-cell-filters.csv")
table(filter.df$F03)
# keep remove
# 80639  12098
stopifnot(identical(filter.df$cell_barcode_sample, unique(colnames(sd.merged))))
sd.merged@meta.data$F03 <- filter.df$F03
# An object of class Seurat


# ---- Run scDblFinder ----
## scDblFinder
set.seed(169)
dbr <- c(0.1)
# sd.names <- c("sd.raw", "sd.ds10k","sd.ds25k")
# sd.names <- c("sd.ds10k","sd.ds25k")
sd.name <- c("sd.merged")
filter.name <- c("")
nfeatures <- 2000


cli::cli_alert_info("processing: {.var {sd.name}}")
cli::cli_alert_info(paste0(
      "  Running scDblFinder on doublet rate: ", cli::col_br_yellow(dbr)))
seurat.object <- get(sd.name)
cli::cli_alert_info(paste0(
  "  Dim Seurat Layered: ", cli::col_br_yellow(dim(seurat.object))
))


run.name <- "dbr10"
cli_alert("Converting Joined Layers to single cell experiment")
sce_obj <- suppressWarnings(seurat.object %>% as.SingleCellExperiment())
cli::cli_alert_info(paste0(
  "  Dim SingleCellExperiment: ", cli::col_br_yellow(dim(sce_obj))
))
gc()

# remove any previous dblFinder results from metadata
u <- grepl("scDblFinder",colnames(colData(sce_obj)))
colData(sce_obj) <- colData(sce_obj)[,!u]

table(SingleCellExperiment::colData(sce_obj)$orig_ident)
start_time <- Sys.time()
cli_alert("Running scDblFinder::scDblFinder: {start_time}")
#BiocParallel::register(BiocParallel::MulticoreParam())
sce_obj <- scDblFinder::scDblFinder(
  # dbr.sd = 1,
  sce_obj,
  samples = SingleCellExperiment::colData(sce_obj)$orig.ident,
  dbr = dbr,
  # BPPARAM=MulticoreParam(4),
  nfeatures = nfeatures,
  verbose = TRUE
)
end_time <- Sys.time()
cli::cli_alert_success(paste0(" Done scDblFinder: ",
  cli::col_blue("(in ", prettyunits::pretty_dt(end_time - start_time), ")"
  )
))
identical(colnames(sce_obj), colnames(seurat.object))
table(sce_obj$scDblFinder.class)
# singlet doublet (note that this differs from the singlet used in manuscript)
# 81285   11452
#BiocParallel::register(BiocParallel::MulticoreParam())
cli_alert("scDblFinder output: {run.name}")
  ## Get the scDbl Finder output columns, rename to match analysis name
sc.cols.i <- grepl("scDblFinder",colnames(colData(sce_obj)))
sc.cols <- colnames(colData(sce_obj))[sc.cols.i]

sc.df <- data.frame(colData(sce_obj)) %>%
  mutate(scDblFinder.class = as.character(scDblFinder.class)) %>%
  dplyr::select(matches("^scDblFinder\\.[^_].*$")) %>%
  # dplyr::select(!scDblFinder.sample) %>%
  rename_all(~ stringr::str_c(.,".",run.name)) %>%
  rename_all(~ stringr::str_replace_all(.x, "[.-]", "_")) %>%
  rename_all(~ stringr::str_replace_all(.x, "_scDblFinder_", "_"))
str(sc.df)
colnames(sce_obj)
sc.df$cell_barcode <- colnames(sce_obj)
table(sc.df$scDblFinder_class)

## ---- save results ----
# write table as csv
readr::write_csv(sc.df, glue("./results/03-3-scDblFinder-results-{run.name}.csv"))

