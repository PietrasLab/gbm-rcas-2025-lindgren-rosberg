# Use the `SeruatDisk::Convert` function package to convert the nUMI downsampled data (from pythin scripts)
# from Scanpy to Seurat object. This is done to ensure that the data is in the correct format for downstream.
# First h5ad is then loaded using `SeuratDisk::LoadH5Seurat` and used to
# creatae a Seurat object using `Seurat::CreateSeuratObject` function.

# ---- Env ----
require(conflicted)
require(tidyverse)
require(qs2)
require(cli)
require(glue)
require(Seurat)


# ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")



# ---- Params ----
data.path = "../data/processed/anndata"
pdata <- readRDS("./metadata/00-sample-metadata-flex-s16.Rds")
# h5.file <- list.files(data.path, pattern = "01-scanpy-n16-dsUMI_.*\\.h5ad", full.names = TRUE)
python.h5ad <- "./data/processed/anndata/01-scanpy-n16-dsUMI_10000.h5"


# ---- Convert downsampled h5 from python script to seurat ----
cli::cli_h1(paste0("Converting ", python.h5ad, " to Seurat object"))
SeuratDisk::Convert(python.h5ad, dest = "h5seurat", overwrite = TRUE)
sd <- SeuratDisk::LoadH5Seurat(gsub(".h5ad", ".h5seurat", python.h5ad))
stopifnot(identical(sd[["RNA"]]@counts, sd[["RNA"]]@data))

# create a seurat object
sd.new <- Seurat::CreateSeuratObject(counts = sd[["RNA"]]@counts, project = "gbm-rcas", min.cells = 0, min.features = 0)
sd.new

# add the pdata object
head(pdata)


# An object of class Seurat
# 19073 features across 92737 samples within 1 assay
# Active assay: RNA (19073 features, 0 variable features)
# 1 layer present: counts


# ---- Add cell qc from python metadata csv file ----
# read the metadata from csv
cli::cli_h1(paste0("Reading metadata from ", gsub(".h5ad", "-metadata.csv", python.h5ad)))
mdata <- read.csv(gsub(".h5ad", "-metadata.csv", python.h5ad))
str(mdata)
colnames(mdata)[1] <- "barcode"
head(mdata)
colnames(sd.new) <- mdata$barcode # add cell ids as column colnames
sd.new
sd.new$cell_barcode <- colnames(sd.new)

# get columns containing "total_counts_" and "total_counts_" followed by nemericals, e.g. "total_counts_1000"
# head(mdata)
# pattern <- "(total_counts_|n_genes_by_counts_)\\d+$"
# u <- grepl(pattern, colnames(mdata))
# mdata <- mdata[,!u]
# head(sd.new@meta.data)
# stopifnot(identical(colnames(sd.new), mdata$cell_barcode))
# sd.new@meta.data$barcode <- colnames(sd.new)

# sd.new@meta.data <- sd.new@meta.data %>%
#   left_join(mdata, by = "barcode") %>%
#   tibble::column_to_rownames("barcode")
# head(sd.new@meta.data)

# Fix sample_id column and set same as orig.ident
sample_ids <- sub("^([^_]+_[^_]+)_.*", "\\1", colnames(sd.new))
table(sample_ids)
sd.new$orig.ident <- factor(sample_ids, levels = unique(sample_ids))
sd.new$sample_id <- factor(sample_ids, levels = unique(sample_ids))

# add the sample pdata to meta.data (match with sample_id)
sd.new@meta.data <- sd.new@meta.data %>%
  left_join(pdata, by = "sample_id") %>%
  tibble::column_to_rownames("cell_barcode")
table(sd.new$sample_id)



# ---- qsave seurat object ----

qs.file <- "./data/processed/seurat/02-2-flex-dsUMI-s16-c92737.qs"
cli::cli_h1(glue("Saving Seurat object to {qs.file}"))
system.time(
  qs2::qs_save(sd.new, qs.file)
)
