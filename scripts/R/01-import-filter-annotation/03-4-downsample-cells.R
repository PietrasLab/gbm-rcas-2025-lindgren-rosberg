# From ~/Projects/gbm-rcas-submission/scripts/xx05-Downsample-Cells-Pynorm-n16-m1v2.Rmd
# Balance the number of cells per sample type
# This by Downsample Recurrent sample n to that of Primary (to 22655)


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


# ---- Load data -----
sd.ds10k <- qs2::qs_read("./data/processed/seurat/02-2-flex-dsUMI-s16-c92737.qs")
colnames(sd.ds10k)
sd.ds10k@meta.data$cell_barcode <- colnames(sd.ds10k)
sd.ds10k$cell_barcode


# ---- Source  ----
source("./scripts/R/00-0-source-functions.R")
pdata <- readRDS("./metadata/00-sample-metadata-flex-s16.Rds")
gene.list <- readRDS(file = "./references/genesets/gene-list.Rds")


# 19073 features across 92737 samples within 1 assay

# ---- Filter data -----
## ---- Mito/Haemo/SampleAdaptive qc ----
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
filter.df
table(filter.df$F03)
# keep remove
# 80639  12098
str(filter.df$cell_barcode_sample[is.na(match(filter.df$cell_barcode_sample, colnames(sd.ds10k)))])

stopifnot(identical(filter.df$cell_barcode_sample, unique(colnames(sd.ds10k))))
sd.ds10k@meta.data$F03 <- filter.df$F03

## ---- Filter scDblFinder cells ----
table(filter.df$scDblFinder_class_dbr10)
# doublet singlet
# 11594   81143
stopifnot(identical(filter.df$cell_barcode, unique(colnames(sd.ds10k))))
dbl.filter <- readr::read_csv("./results/03-qc/03-2-cell-filters.csv")
sd.ds10k@meta.data$dbr10 <- filter.df$scDblFinder_class_dbr10
table(sd.ds10k@meta.data$dbr10, sd.ds10k@meta.data$F03)
# keep remove
# doublet  8494   3100
# singlet 72145   8998
sd.ds10k
# An object of class Seurat
# 19073 features across 92737 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 1 layer present: counts

## ---- Subset the data -----
sd.ds10k <- subset(sd.ds10k, subset = F03 == "keep" & dbr10 == "singlet")
sd.ds10k
# An object of class Seurat
# 19073 features across 72145 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 1 layer present: counts


# ---- cells per sample tyoe distribution ----
mdata <- sd.ds10k@meta.data

table(mdata$sample_type)
# Healthy   Primary Recurrent
#   19494     22655     29996
# 5 primary
# 6 recurrent
table(sd.ds10k@meta.data$sample_type)
# Healthy   Primary Recurrent
# 19494     22655     29996


# ---- Downsample Cells ----
set.seed(169)
# separate into subsets of primary/recurrent
sd.merged.p <- subset(sd.ds10k, subset = sample_type == "Primary")
# 22655 samples
sd.merged.r <- subset(sd.ds10k, subset = sample_type == "Recurrent")
# 29996 samples
sd.merged.r
# set active ident to orig,ident and downsample toi the lowest number
sd.merged.r <- SetIdent(sd.merged.r, value = "orig.ident")
table(Idents(sd.merged.r))
# TR_01 TR_02 TR_03 TR_04 TR_06 TR_07
# 4642  5420  4800  5619  4840  4675
sd.merged.r <-  subset(x = sd.merged.r, downsample = 4642)
# 29230 samples
table(sd.merged.r@meta.data$orig.ident)
# set active ident to orig,ident and downsample
# balance the Recurrent total n to Primary total n 22655
sd.merged.r <- SetIdent(sd.merged.r, value = "sample_type")
table(Idents(sd.merged.r))
sd.merged.r <-  subset(x = sd.merged.r, downsample = 22655)
table(sd.merged.r@meta.data$orig.ident)

# TR_01 TR_02 TR_03 TR_04 TR_06 TR_07
#  3966  3932  3962  3946  3958  3932
table(sd.merged.r$sample_type)
# 22655 samples

mdata <- mdata %>%
  dplyr::mutate(
    C1 = ifelse(
      sample_type %in% c("Healthy", "Primary") | cell_barcode %in% colnames(sd.merged.r),
      "keep",
      "remove"
    )
  )
str(mdata)
table(mdata$C1, useNA="always")
table(mdata$orig.ident, mdata$C1)


# ---- Create Data frame with all filters, including Downsampling ----
# ---- Save C1 - Downsampled cell dataframe ----
df <- filter.df %>%
  dplyr::select(-cell_barcode) %>%
  dplyr::rename("cell_barcode" = cell_barcode_sample) %>%
  # scDblFinder_class_dbr10 = "dbr10", set dbr10 keep && F03 keep F03_dbr10
  dplyr::left_join(dbl.filter, by = "cell_barcode") %>%
  dplyr::mutate(
    F03dbr10 = if_else(
      grepl("F03", cell_barcode) & grepl("dbr10", scDblFinder_class_dbr10),
      "keep","remove")) %>%
  dplyr::left_join(mdata %>% dplyr::select(cell_barcode, C1), by = c("cell_barcode")) %>%
  # set NA ro remove for C1
  dplyr::mutate(
    final_cell_filter = ifelse(is.na(C1), "remove", C1)
  )
  # remove the F01
colnames(df)
df <- df[, c("cell_barcode", "F03","scDblFinder_class_dbr10","final_cell_filter")]
tail(as.data.frame(df))
table((df$final_cell_filter))
# keep remove
# 64804  27933
# ---- save the filter data frame
readr::write_csv(df, file = "./results/03-qc/03-4-cell-filters-downsample.csv")
