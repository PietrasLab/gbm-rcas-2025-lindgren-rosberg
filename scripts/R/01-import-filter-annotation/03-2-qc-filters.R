
# ---- Summary ----
# Evaluate impact on different filtering steps on the Flex data.
# The step used in manusript is F03


# ---- Env ----
require(conflicted)
require(tidyverse)
require(qs2)
require(cli)
require(glue)
require(readr)
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
# library(scCustomize)
# require(peRcebe)
require(CellMetaVerse)
require(kableExtra)
require(scCustomize)


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



# ---- QC PLOTS dir ----
plot.dir <- "./results/03-qc"
if (!dir.exists(plot.dir)) {
  dir.create(plot.dir, recursive = TRUE)
}


# ---- Filter F, Mito and Haemo Filter ----
## ---- define filter dataframe ----
filt.df <- data.frame(
  cell_barcode = colnames(sd.merged),
  cell_barcode_sample = sd.merged@meta.data$cell_barcode,
  F01 = 'remove', F02 = 'remove', F03 = 'remove')
head(filt.df)



## ---- F01: Fixed cutoffs ----
# Remove cells with high mitochondrial and haemoglobin content
# Cell filter F01. Only cells
# percent_mitochondria = "percent_mitochondria < 10",
# percent_haemoglobin = "percent_haemoglobin < 5",
# F01: Mito and Haemo Filter
# Cell filter F01. Only cells
F01 <- .CreateSeuratFilter(
  filter.name = "cell-filter-fixed-cutoffs",
  cell.filters = list(
    percent_mitochondria = "percent_mitochondria < 10",
    percent_haemoglobin = "percent_haemoglobin < 5"
  )
)
colnames(sd.merged@meta.data)
# str(sd.f01@meta.data)
sd.f01 <- .seuratFilterFoo(
  sd.merged,
  filter.object = F01
)
#   19073 genes and 90816 cells
# percent_mitochondria = "percent_mitochondria < 10",
# percent_haemoglobin = "percent_haemoglobin < 5",

## Plot F01. nCount and nFeatures
## nCount: number of UMIs per sample
df <- sd.f01@meta.data
df <- df %>%
  dplyr::select(orig.ident, nCount_RNA) %>%
  dplyr::group_by(orig.ident) %>%
  summarize(sum_nCount = sum(nCount_RNA))

my.pal <- .palette_discrete_bp(df$orig.ident)

ggplot(
  df,
  aes(x = orig.ident,
    y = sum_nCount,
    fill = orig.ident)) +
  scale_fill_manual("sum_nCount", values = my.pal)  +
  geom_bar(stat = "identity", color="black") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))


## nCount: Ridgeplot nCount. number of UMIs
RidgePlot(sd.f01,
  features = "nCount_RNA",
  group.by = "orig.ident") +
  guides(fill = "none") + scale_fill_viridis_d()

## nCount: Violin. number of UMIs per cell
p <- SCpubr::do_ViolinPlot(
  sample = sd.f01, y_cut = 500,
  features = "nCount_RNA", plot_boxplot = F
)
p
# plot logged y axis
p + scale_y_continuous(trans='log10') +
  labs(y = "nUMI (log)")


## Binned cutoff - nCount
i.name <- "nCount_RNA"
b.name <- paste0(i.name, "_bin")
my.breaks <-  c(0,100,1000,10000,50000, 100000,Inf)
my.labels <- paste(
  as.character(my.breaks)[-length(my.breaks)],
  as.character(my.breaks)[-1],
  sep="-")

metadata <- sd.f01@meta.data
df <- metadata %>%
  dplyr::select(orig.ident, !!rlang::sym(i.name)) %>%
  mutate(b.name = cut(!!rlang::sym(i.name),
    breaks = my.breaks,
    labels = my.labels,
    include.lowest = T, right = F))
u <- which(colnames(df) == "b.name")
colnames(df)[u] <- b.name

.barplot_stacked(
  df,
  group.var = "orig.ident",
  color.var = b.name,
  scaled.y = F)

## nFeatures: the number of genes detected in each cell
# plot nFeatures viloin
p <- SCpubr::do_ViolinPlot(
  sample = sd.f01, y_cut = 500,
  features = "nFeature_RNA", plot_boxplot = F
)
p
# plot logged y axis
p + scale_y_continuous(trans='log10') +
  labs(y = "nGene (log)")

RidgePlot(sd.f01,
  features = "nFeature_RNA",
  group.by = "orig.ident") +
  guides(fill = "none") + scale_fill_viridis_d()

## Genes UMIs and Complexity plots
scCustomize::QC_Plots_Genes(
  sd.f01, pt.size = 0, low_cutoff = 1000,
  high_cutoff = 10000)
scCustomize::QC_Plots_UMIs(
  seurat_object = sd.f01,
  low_cutoff = 1200, high_cutoff = 45000,
  pt.size = 0, y_axis_log = T)
scCustomize::QC_Plots_Complexity(
  seurat_object = sd.f01,
  high_cutoff = 0.8, pt.size = 0)


scCustomize::QC_Plot_GenevsFeature(
  seurat_object = sd.f01,
  feature1 = "percent_mitochondria",
  low_cutoff_gene = 1000,
  high_cutoff_gene = 10000, high_cutoff_feature = 10)

# Plot individual QC_Plot_GenevsFeature vs Mitochondrial load...
fname <- paste0("03-2-F01-FixedMitoHaemoFilter.pdf")
file.name <- file.path(plot.dir, fname)
cli_alert_info("plotting to file {.file {file.name}}")

pdf(file = file.name, width=8, height = 7)
for(i in sort(unique(Idents(sd.f01)))){
  p <- QC_Plot_GenevsFeature(
    seurat_object = subset(sd.f01, idents = i),
    feature1 = "percent_mitochondria",
    low_cutoff_gene = 1000, high_cutoff_gene = 10000,
    raster = TRUE,
    high_cutoff_feature = 10)
  plot(p + ggtitle(label = paste("Sample",i)))

}
dev.off()

u <- match(colnames(sd.f01), filt.df$cell_barcode)
str(u)
filt.df[u,]$F01 <- 'keep'
table(filt.df$F01)

## perform UMAP
nfeatures <- 2000
dims_pca <- 30
## AUTO SET PARAMS
regress.name <- "none"
process.name <- "n16_F01"
reduction.name = paste0("umap_", process.name)

# normalize <- TRUE # Always TRUE
# scale.data <- TRUE # always TRUE
my.description <- glue::glue(
  "Filter: raw.
          Process: {process.name}.
          Reduction.name: {reduction.name}
          Norm: TRUE. VarFeatures: {nfeatures}.
          Scale: TRUE
          Regress: {regress.name}
          Dims: {dims_pca}
                "
)
cli_alert(my.description)

cli_alert_info("Run Seurat - Normalize - FindVar - Scale - PCA - UMAP")
seurat.object <- sd.f01 %>%
  SeuratObject::JoinLayers() %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nfeatures) %>%
  Seurat::ScaleData(verbose = T) %>%
  Seurat::RunPCA(verbose=TRUE, npcs = dims_pca) %>%
  Seurat::RunUMAP(reduction.name = reduction.name,
    dims = 1:dims_pca, verbose = T)


colnames(seurat.object@meta.data)
annot <- "orig.ident"
p <- SCpubr::do_DimPlot(
  seurat.object,
  # idents.keep = cell.barcodes,
  # split.by = "input_subset",
  group.by = annot,
  colors.use = .color_pal[["orig_ident"]],
  # legend.position = "none",
  # cells.highlight = cell.barcodes,
  label = FALSE,
  plot_cell_borders = FALSE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)
p
# ggsave()

annot <- "nFeature_RNA"
SCpubr::do_FeaturePlot(
  features = annot,
  seurat.object,
  pt.size = 2,
  raster = TRUE,
  verbose = FALSE,
  plot_cell_borders = TRUE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)


## ---- F02: Mito/Haemo and Fixed nFeatureCount cutoff ----
F02 <- peRcebe_CreateSeuratFilter(
  filter.name = "cell-filter-fixed-cutoffs",
  cell.filters = list(
    percent_mitochondria = "percent_mitochondria < 10",
    percent_haemoglobin = "percent_haemoglobin < 5",
    nFeature_low = "nFeature_RNA > 500",
    nFeature_high = "nFeature_RNA < 10000",
    nCount_low = "nCount_RNA > 1000",
    nCount_high = "nCount_RNA < 100000"
    # scDblFinder = "scDblFinder.class != 'doublet'",
    # decontX_call = "decontX_call_025 == 'dcontx_cell'"
  )
)

# str(sd.f01@meta.data)
sd.f02 <- .seuratFilterFoo(
  sd.merged,
  filter.object = F02
)

# ---❯  19073 genes and 86158 cells

## Add the fixed threshold filter F01 to database
u <- match(colnames(sd.f02), filt.df$cell_barcode)
str(u)
filt.df[u,]$F02 <- 'keep'
table(filt.df$F02)

## add UMAP F02
nfeatures <- 2000
dims_pca <- 30

## AUTO SET PARAMS
regress.name <- "none"
process.name <- "n16_F02"
reduction.name = paste0("umap_", process.name)

cli_alert_info("Run Seurat - Normalize - FindVar - Scale - PCA - UMAP")
seurat.object <- sd.f02 %>%
  SeuratObject::JoinLayers() %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nfeatures) %>%
  Seurat::ScaleData(verbose = T) %>%
  Seurat::RunPCA(verbose=TRUE, npcs = dims_pca) %>%
  Seurat::RunUMAP(reduction.name = reduction.name,
    dims = 1:dims_pca, verbose = T)

colnames(seurat.object@meta.data)
annot <- "orig.ident"
p <- SCpubr::do_DimPlot(
  seurat.object,
  # idents.keep = cell.barcodes,
  # split.by = "input_subset",
  group.by = annot,
  colors.use = peRcebe_percebe_color_pal[["orig_ident"]],
  # legend.position = "none",
  # cells.highlight = cell.barcodes,
  label = FALSE,
  plot_cell_borders = FALSE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)
p
ggsave(p, filename = "./results/03-qc/03-2-F02-umap-orig-ident.png",
  width = 10, height = 8, units = "in", dpi = 300
)
annot <- "nFeature_RNA"
SCpubr::do_FeaturePlot(
  features = annot,
  seurat.object,
  pt.size = 2,
  raster = TRUE,
  verbose = FALSE,
  plot_cell_borders = TRUE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)


## ---- F03: individual Quantiles Filter ----
# Quantiles from stats produced in QC scripts.
# Save the list of cells that are kept
metadata.stats <- readr::read_delim(file = "./results/03-qc/03-1-qc-samples-seurat-raw-s16-c92737.csv")
metadata.stats
str(metadata.stats)

## Loop and filter all samples
# nFeature - nGene q05/q95
# nCount - nUMI q05/q95
colnames(sd.merged@meta.data)
sd.list <- SplitObject( sd.merged )
metadata.stats$nGenes_q05
metadata.stats$nGenes_q95
metadata.stats$nUMI_q05
metadata.stats$nUMI_q95

for (i in 1:length(sd.list) ) {
  u <- match(names(sd.list)[i], metadata.stats$orig.ident)
  filter.i <- .CreateSeuratFilter(
    cell.filters = list(
      percent_mitochondria = "percent_mitochondria < 10",
      percent_haemoglobin = "percent_haemoglobin < 5",
      nFeature_low = paste0("nFeature_RNA > ", metadata.stats$nGenes_q05[u]),
      nFeature_high = paste0("nFeature_RNA < ", metadata.stats$nGenes_q95[u]),
      nCount_low = paste0("nCount_RNA > ",metadata.stats$nUMI_q05[u]),
      nCount_high = paste0("nCount_RNA < ",metadata.stats$nUMI_q95[u])
    ))

  sd.list[[i]] <- .seuratFilterFoo(sd.list[[i]], filter.i)
  # sd.layered.filt[[i]] <- .seuratFilterFoo(sd.layered.filt[[i]], F02)
}

## Join all layers
# generate a string to merge data - for if >2 samples

## get all cells that are kept
cells.f03 <- unlist(lapply(sd.list, colnames))
table(duplicated(cells.f03))
sd.f03 <- subset(sd.merged, cells=cells.f03)
sd.f03
u <- match(colnames(sd.f03), filt.df$cell_barcode)
str(u)
filt.df[u,]$F03 <- 'keep'
table(filt.df$F03)

# keep remove
# 80639  12098
head(filt.df)

nfeatures <- 2000
dims_pca <- 30

## AUTO SET PARAMS
regress.name <- "none"
process.name <- "n16_F03"
reduction.name = paste0("umap_", process.name)

cli_alert_info("Run Seurat - Normalize - FindVar - Scale - PCA - UMAP")
seurat.object <- sd.f03 %>%
  SeuratObject::JoinLayers() %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = nfeatures) %>%
  Seurat::ScaleData(verbose = T) %>%
  Seurat::RunPCA(verbose=TRUE, npcs = dims_pca) %>%
  Seurat::RunUMAP(reduction.name = reduction.name,
    dims = 1:dims_pca, verbose = T)


colnames(seurat.object@meta.data)
annot <- "orig.ident"
p <- SCpubr::do_DimPlot(
  seurat.object,
  # idents.keep = cell.barcodes,
  # split.by = "input_subset",
  group.by = annot,
  colors.use = .color_pal[["orig_ident"]],
  # legend.position = "none",
  # cells.highlight = cell.barcodes,
  label = FALSE,
  plot_cell_borders = FALSE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)
ggsave(p, filename = "./results/03-qc/03-2-F03-umap-orig-ident.png",
  width = 10, height = 8, units = "in", dpi = 300
)

annot <- "orig.ident"
SCpubr::do_DimPlot(
  seurat.object,
  # idents.keep = cell.barcodes,
  # split.by = "input_subset",
  group.by = annot,
  colors.use = .color_pal[["orig_ident"]],
  # legend.position = "none",
  # cells.highlight = cell.barcodes,
  label = FALSE,
  plot_cell_borders = FALSE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)
annot <- "sample_type"
SCpubr::do_DimPlot(
  seurat.object,
  # idents.keep = cell.barcodes,
  # split.by = "input_subset",
  group.by = annot,
  colors.use = .color_pal[[annot]],
  # legend.position = "none",
  # cells.highlight = cell.barcodes,
  label = FALSE,
  plot_cell_borders = FALSE
  #plot_marginal_distributions = TRUE
) + ggtitle(annot)


# --- Save Filter Table ----
str(filt.df)
readr::write_csv(filt.df, "./results/03-qc/03-2-cell-filters.csv")



