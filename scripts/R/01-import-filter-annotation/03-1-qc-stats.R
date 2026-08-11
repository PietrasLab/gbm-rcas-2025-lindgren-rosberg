# Add Some QC basic stats
## Code are taken from example tutorials:
# * [scCustomize](https://samuel-marsh.github.io/scCustomize/index.html).
# * [SCPubRr](https://enblacar.github.io/SCpubr-book-v1/).
# * [NBIS workshop](https://nbisweden.github.io/workshop-scRNAseq/exercises.html).
# * [Harvard Chan Bioinformatics Core](https://hbctraining.github.io/main/) and page [SC  qc](https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html).
# scCustomize contains easy shortcut function to add a measure of cell complexity/novelty that can sometimes be useful to filter low quality cells. The metric is calculated by calculating the result of log10(nFeature) / log10(nCount).


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

# ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")


# ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
results.dir <- "./results/03-qc"
dir.create(results.dir, showWarnings = F, recursive = T)

# ---- Load data -----
# 02-1-flex-raw-s16-c92737.qs is a Seurat object containing raw data from 16 samples
# same as seurat_flex_raw_v1.0.rds on Zenodo
# This can be downloaded from Zenodo:
# https://zenodo.org/records/21607378
# Zenodo file: seurat_flex_raw_v1.0.rds (check for latest version)

# seurat.object <- readRDS("./data/processed/seurat/seurat_flex_raw_v1.0.rds")
seurat.object <- qs2::qs_read("./data/processed/seurat/02-1-flex-raw-s16-c92737.qs")
# An object of class Seurat
# 19073 features across 92737 samples within 1 assay

# ---- Load probe feature data ----
fdata <- read.delim("./metadata/00-probeset-flex.csv", header = T, sep=",")
str(fdata)
head(fdata)
stopifnot(identical(fdata$external_gene_name, rownames(seurat.object)))

## ---- Load a gene list ----
# gene list with transgenes, mitochondrial genes, haemoglobin genes, platelet genes, rcas flex probes etc
gene.list <- readRDS(file = "./references/genesets/gene-list.Rds")


# ---- Collect QC stats -----

## Gene Modules: Percentage expression of probes (scCustomize)
# Load the gene signatures from file

## ---- scCustomize ----
cli::cli_alert("Adding scCustomize::Add_Cell_Complexity")
seurat.object <-  scCustomize::Add_Cell_Complexity(
  object = seurat.object, overwrite = TRUE )
head(seurat.object@meta.data$log10GenesPerUMI)


## ---- PercentageFeatureSet ----
## add percents for custom gene lists
feature_vec <- c(
  "mitochondria","haemoglobin","platlets","rcas_flex_probes",
  "rcas_flex_probes_pdgfb","rcas_flex_probes_rfp")
## loop all feature vecs (gene lists). Calculate percentage of probes
for ( g in  feature_vec){
  g.name <- paste0("percent_",g)
  if (
    length(gene.list[[g]] > 1 )
    # && !(g.name %in% colnames(seurat.object@meta.data))
  ){
    cli_alert("calculate PercentageFeatureSet {.var {g.name}}")
    seurat.object <- Seurat::PercentageFeatureSet(
      seurat.object, features = gene.list[[g]],
      col.name = g.name)
  } else {
    cli_alert("skipping {.var {g.name}}")
    next ( g )
  }}
# str(tail(seurat.object@meta.data))

# ---- Save QC stats to tables ----
metadata <- seurat.object@meta.data
str(metadata)

# Save metadata to csv file
write.csv(metadata, file = file.path(results.dir, "03-1-seurat-raw-s16-c92737-metadata.csv", row.names = TRUE))
# Update seurat object on disk
qs2::qs_save(seurat.object, file = "./data/processed/seurat/02-1-flex-raw-s16-c92737.qs")


# ---- Basic QC per sample ----
# sample-wise QC metrics, used for setting sample-adaptive filtering thresholds
# Hquantile nUMI data that may be used for subsequent filtering of individual data sets.
## Basic QC per sample (mean/median as pdata)
#   From: https://link.springer.com/article/10.1186/s13059-021-02584-9 :
# we filtered out cells that fell outside of the 5 to 95% UMI quantiles in each dataset
metadata <- seurat.object@meta.data %>%
  group_by(orig.ident)
# str(metadata)
# colnames(metadata)

df <- metadata %>%
  dplyr::filter(orig.ident=="TR_07")
quantile(df$nCount_RNA, probs = c(0.05,0.95))
quantile(df$nCount_RNA, probs = c(0.01,0.99))


metadata.stats <- metadata %>%
  summarise(
    nCells = n(),
    nUMI_sum = sum(nCount_RNA, na.rm=T),
    nUMI_mean = mean(nCount_RNA, na.rm=T),
    nUMI_median = median(nCount_RNA, na.rm=T),
    nUMI_q01 = quantile(nCount_RNA, na.rm=T, probs=c(0.01)),
    nUMI_q05 = quantile(nCount_RNA, na.rm=T, probs=c(0.05)),
    nUMI_q95 = quantile(nCount_RNA, na.rm=T, probs=c(0.95)),
    nUMI_q99 = quantile(nCount_RNA, na.rm=T, probs=c(0.99)),
    nGenes_sum = sum(nFeature_RNA, na.rm=T),
    nGenes_mean = mean(nFeature_RNA, na.rm=T),
    nGenes_median = median(nFeature_RNA, na.rm=T),
    nGenes_q01 = quantile(nFeature_RNA, na.rm=T, probs=c(0.01)),
    nGenes_q05 = quantile(nFeature_RNA, na.rm=T, probs=c(0.05)),
    nGenes_q95 = quantile(nFeature_RNA, na.rm=T, probs=c(0.95)),
    nGenes_q99 = quantile(nFeature_RNA, na.rm=T, probs=c(0.99)),
    log10GenesPerUMI_mean = mean(log10GenesPerUMI),
    log10GenesPerUMI_median = median(log10GenesPerUMI),
    percent_hb_mean = mean(percent_haemoglobin),
    percent_hb_median = median(percent_haemoglobin),
    percent_mito_mean = mean(percent_mitochondria),
    percent_mito_median = median(percent_mitochondria)
  )

metadata.stats[,7:13]
# Save metadata stats to csv file
write.csv(metadata.stats, file = "./results/03-qc/03-1-qc-samples-seurat-raw-s16-c92737.csv", row.names = TRUE)


# ---- QC PLOTS ----
pdata <- readRDS("./metadata/00-sample-metadata-flex-s16.Rds")
plot.dir <- "./results/03-qc"
if (!dir.exists(plot.dir)) {
  dir.create(plot.dir, recursive = TRUE)
}

## ----  Total & median cell number ----
ncells <- metadata.stats %>%
  summarise(nCells_per_id = sum(nCells))
glue::glue("Total nCells all {nrow(pdata)} samples: {ncells}")
ncell_median <- metadata.stats %>%
  summarise(nCells_per_id = median(nCells))
glue::glue("Median nCells all {nrow(pdata)} samples: {ncell_median}")


## ---- Plot basic cell numbers per experiment ----

metadata %>%
  ggplot(aes(x=orig.ident, fill=orig.ident)) +
  geom_bar() +
  theme_classic() +
  lims(y = c(0, 15000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") +
  theme(legend.position = "none") +
  geom_hline(yintercept = as.numeric(ncell_median)) +
  geom_hline(yintercept = 5000, linetype='dashed')

p <- Plot_Cells_per_Sample(
  seurat_object = seurat.object,
  group_by = "sample_type")
p +
  lims(y = c(0, 15000)) +
  geom_hline(yintercept = as.numeric(ncell_median)) +
  geom_hline(yintercept = 5000, linetype='dashed')

ggsave(p, filename = file.path(plot.dir, "03-1-nCells-sample.png"),
       width = 8, height = 6, units = "in", dpi = 300)


## ---- Plot nCount and nFeatures ----

### ---- nCount: number of UMIs per sample  ----
df <- seurat.object@meta.data
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


### ----  nCount: Ridgeplot nCount. number of UMIs ----
RidgePlot(seurat.object,
  features = "nCount_RNA",
  group.by = "orig.ident") +
  guides(fill = "none") + scale_fill_viridis_d()

### ----   nCount: Violin. number of UMIs per cell ----
p <- SCpubr::do_ViolinPlot(
  sample = seurat.object, y_cut = 500,
  features = "nCount_RNA", plot_boxplot = F
)
p
# plot logged y axis
p + scale_y_continuous(trans='log10') +
  labs(y = "nUMI (log)")
ggsave(p, filename = file.path(plot.dir, "03-1-nCount-sample-violin.png"),
       width = 8, height = 6, units = "in", dpi = 300)



### ---- nFeatures: the number of genes detected in each cell ----
# plot nFeatures viloin
p <- SCpubr::do_ViolinPlot(
  sample = seurat.object, y_cut = 500,
  features = "nFeature_RNA", plot_boxplot = F
)
p
# plot logged y axis
p <- p + scale_y_continuous(trans='log10') +
  labs(y = "nGene (log)")
p
ggsave(p, filename = file.path(plot.dir, "03-1-nFeatures-per-cell-log.png"),
       width = 8, height = 6, units = "in", dpi = 300)

p <- RidgePlot(seurat.object,
  features = "nFeature_RNA",
  group.by = "orig.ident") +
  guides(fill = "none") + scale_fill_viridis_d()
p
ggsave(p, filename = file.path(plot.dir, "03-1-nFeatures-per-sample-ridge.png"),
       width = 8, height = 6, units = "in", dpi = 300)





### ---- Binned cutoff - nFeatures ----
i.name <- "nFeature_RNA"
b.name <- paste0(i.name, "_bin")
my.breaks <-  c(0,250,500,1000,5000,7500,10000,Inf)
my.labels <- paste(
  as.character(my.breaks)[-length(my.breaks)],
  as.character(my.breaks)[-1],
  sep="-")

metadata <- seurat.object@meta.data
df <- metadata %>%
  dplyr::select(orig.ident, !!rlang::sym(i.name)) %>%
  mutate(b.name = cut(!!rlang::sym(i.name),
    breaks = my.breaks,
    labels = my.labels,
    include.lowest = T, right = F))
u <- which(colnames(df) == "b.name")
colnames(df)[u] <- b.name

p <- .barplot_stacked(
  df,
  group.var = "orig.ident",
  color.var = b.name,
  scaled.y = F)
p
ggsave(p, filename = file.path(plot.dir, "03-1-nFeatures-sample-binned.png"),
       width = 8, height = 6, units = "in", dpi = 300)

p <- .barplot_stacked(
  df,
  group.var = "orig.ident",
  color.var = b.name,
  scaled.y = T)
p
ggsave(p, filename = file.path(plot.dir, "03-1-nFeatures-sample-binned-scaled.png"),
       width = 8, height = 6, units = "in", dpi = 300)



### ---- Features per cell ----
str(pdata)
# str(seurat.object@assays$RNA@counts)
my.mat <- sapply(pdata$sample_id, function(x){
  # Matrix::rowSums(subset(seurat.object, idents = x)[["RNA"]]@layers >= 1)
  Matrix::rowSums(GetAssayData(subset(seurat.object, idents=x), assay = "RNA") >= 1)
})
str(my.mat)
my.melt <- reshape2::melt(my.mat)
str(my.melt)
my.melt[,c("Var2","value")]

p <- .violin(
  plot_df = my.melt,
  group.var = "Var2",
  score.var = "value") +
  ggtitle(label = paste(
    "nCells per Gene and Sample
    with non-zero values"))
p
p <- .violin(
  plot_df = my.melt,
  group.var = "Var2",
  score.var = "value",
  log_y = T) +
  ggtitle(label = paste(
    "nCells per Gene and Sample
    with non-zero values"))

str(my.melt)

df <- my.melt %>%
  mutate(nCells_bin = cut(value,
    breaks = c(0,1,5,10,25,100,1000,Inf),
    labels = c("0","1-5","5-10","10-25","25-100","100-1000",">1000"),
    include.lowest = T, right = T))


p <- .barplot_stacked(
  df,
  group.var = "Var2",
  color.var = "nCells_bin",
  scaled.y = F) +
  ggtitle(label = paste("nCells with non-zero values per Gene")) +
  labs(x = "orig.ident", y = "nGenes")
p
ggsave(p, filename = file.path(plot.dir, "03-1-nCells-gene-binned.png"),
       width = 8, height = 6, units = "in", dpi = 300)


# ----- Check adaptive QC quantile cutoffs ----
# Sample QC cutoffs:
# adaptive vs fixed thresholds

.qcCutoffFoo <- function(my.sample, metadata, metadata.stats){
  table(metadata$orig_ident)
  metadata.stats <- metadata.stats %>%
    dplyr::mutate(orig_ident = orig.ident)
  n.orig <- metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    nrow()

  # nUMI filter (0.01)
  n.filter.99 <- n.orig - metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    dplyr::filter(! nCount_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q01"])) %>%
    dplyr::filter(! nCount_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q99"])) %>%
    nrow()

  # Filter both nGenes and nUMI
  n.all.99 <- n.orig - metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    dplyr::filter(! nCount_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q01"])) %>%
    dplyr::filter(! nCount_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q99"])) %>%
    dplyr::filter(! nFeature_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nGenes_q01"])) %>%
    dplyr::filter(! nFeature_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nGenes_q99"])) %>%
    nrow()

  # nUMI filter (10.05
  n.filter.95 <- n.orig - metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    dplyr::filter(! nCount_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q05"])) %>%
    dplyr::filter(! nCount_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q95"])) %>%
    nrow()

  # Filter both nGenes and nUMI
  n.all.95 <- n.orig - metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    dplyr::filter(! nCount_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q05"])) %>%
    dplyr::filter(! nCount_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q95"])) %>%
    dplyr::filter(! nFeature_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nGenes_q05"])) %>%
    dplyr::filter(! nFeature_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nGenes_q95"])) %>%
    nrow()

  # Filter both nGenes and nUMI
  n.mix <- n.orig - metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    dplyr::filter(! nCount_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q05"])) %>%
    dplyr::filter(! nCount_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nUMI_q95"])) %>%
    dplyr::filter(! nFeature_RNA < as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nGenes_q01"])) %>%
    dplyr::filter(! nFeature_RNA > as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, "nGenes_q99"])) %>%
    nrow()

  ## Standard fixed filter cutoff

  n.standard <- n.orig - metadata %>%
    dplyr::filter(orig_ident==my.sample) %>%
    dplyr::filter(! nCount_RNA < 1000) %>%
    dplyr::filter(! nCount_RNA > 100000) %>%
    dplyr::filter(! nFeature_RNA < 500) %>%
    dplyr::filter(! nFeature_RNA > 10000) %>%
    nrow()

  require(glue)
  cli_alert(" ")
  cli_alert("==============")
  cli_alert_info("Sample {my.sample}")
  cli_alert_info("nUMI_q05: {floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nUMI_q05']))}-{floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nUMI_q95']))}, nGenes_q05 {floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nGenes_q05']))}-{floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nGenes_q95']))}")
  cli_alert_info("nUMI_q01: {floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nUMI_q01']))}-{floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nUMI_q99']))}, nGenes_q01 {floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nGenes_q01']))}-{floor(as.numeric(metadata.stats[metadata.stats$orig_ident==my.sample, 'nGenes_q99']))}")
  cli_alert_info("Total nCells raw: {n.orig}")
  # cli_alert_info("Total nCells filter nCount 0.01 and 0.99 perc: {n.filter.99}")
  cli_alert_info("Total nCells filter nCount AND nFeatures 0.01/0.99 perc: {n.all.99}. ({round(n.all.99/n.orig,2)*100}%)")
  cli_alert_info("Total nCells filter nCount AMD nFeatures 0.05/0.95 perc: {n.all.95} ({round(n.all.95/n.orig,2)*100}%)")
  cli_alert_info("Total nCells filter nCount 0.05/0.95 and nFeatures 0.01/0.99  perc: {n.mix} ({round(n.mix/n.orig,2)*100}%)")
  cli_alert_info("----")
  cli_alert_info("if using standard filter w fixed cutoffs: {n.standard} ({round(n.standard/n.orig,2)*100}%)")
  cli_alert_info("==============")
}

## ----- run qcCutoffFoo -----
# check quantiles - RAW DATA
# chech stats for individual samples before filtering
metadata = seurat.object@meta.data
metadata <- metadata %>% dplyr::mutate(orig_ident = orig.ident)
str(metadata)
for(my.sample in as.character(pdata$sample_id)){
  cli::cli_alert("processing {my.sample}")
  .qcCutoffFoo(my.sample, metadata, metadata.stats)
}



### ---- Genes UMIs and Complexity ----

suppressWarnings( scCustomize::QC_Plots_Genes(
  seurat.object, pt.size = 0, low_cutoff = 1000,
  high_cutoff = 10000))
suppressWarnings(scCustomize::QC_Plots_UMIs(
  seurat_object = seurat.object,
  low_cutoff = 1200, high_cutoff = 45000,
  pt.size = 0, y_axis_log = T))
suppressWarnings(scCustomize::QC_Plots_Complexity(
  seurat_object = seurat.object,
  high_cutoff = 0.8, pt.size = 0))


### ---- Percent of total expression ----
#  Percent Mitochondrial transcripts. A common cutoff is 15-20%
table(seurat.object@meta.data$percent_mitochondria > 20)
table(seurat.object@meta.data$percent_mitochondria > 5)
head(seurat.object@meta.data)
Seurat::RidgePlot(seurat.object,
  y.max = 50,
  features = "percent_mitochondria",
  group.by = "orig.ident") +
  guides(fill = "none") + scale_fill_viridis_d()


i.name <- "percent_mitochondria"
b.name <- paste0(i.name, "_bin")
my.breaks <-  c(0,1,2,5,10,Inf)
my.labels <- paste(
  as.character(my.breaks)[-length(my.breaks)],
  as.character(my.breaks)[-1],
  sep="-")

metadata <- seurat.object@meta.data
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

p <- .barplot_stacked(
  df,
  group.var = "orig.ident",
  color.var = b.name,
  scaled.y = T)
ggsave(p,
  filename = file.path(plot.dir, glue("03-1-percent-mito-bin.png")),
  width = 6, height = 6)






# High Hemoglobin expression per cell (>10-20%) may be
# indication red blood cell contamination
i.name <- "percent_haemoglobin"
b.name <- paste0(i.name, "_bin")
my.breaks <- c(0,0.1,0.25,0.5,1,Inf)
my.labels <- paste(
  as.character(my.breaks)[-length(my.breaks)],
  as.character(my.breaks)[-1],
  sep="-")

metadata <- seurat.object@meta.data
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
.barplot_stacked(
  df,
  group.var = "orig.ident",
  color.var = b.name,
  scaled.y = T)


## ---- Stats when filtering Mito/Haemo ----
## Example Sample QC cutoffs:
# requires:
#   metadata.stats table from above (precalculated quantile cuts)
#   metadata - table with cell ids, sample id, cell nFeatures n Genes
.qcCutoffFoo <- function(my.sample, mito.cut=10, haemo.cut=5){
  n.orig <- metadata %>%
    dplyr::filter(orig.ident==my.sample) %>%
    nrow()

  # nUMI filter (0.01)
  # mito.cut <- 10
  n.filter.mito <- n.orig - metadata %>%
    dplyr::filter(orig.ident==my.sample) %>%
    dplyr::filter(! percent_mitochondria > mito.cut) %>%
    nrow()

  # nUMI filter (0.01)
  # haemo.cut <- 5
  n.filter.haemo <- n.orig - metadata %>%
    dplyr::filter(orig.ident==my.sample) %>%
    dplyr::filter(! percent_haemoglobin > haemo.cut) %>%
    nrow()

  n.filter.both <- n.orig - metadata %>%
    dplyr::filter(orig.ident==my.sample) %>%
    dplyr::filter(! percent_mitochondria > mito.cut) %>%
    dplyr::filter(! percent_haemoglobin > haemo.cut) %>%
    nrow()

  require(glue)
  cli_alert(" ")
  cli_alert("==============")
  cli_alert_info("Sample {my.sample}")

  # cli_alert_info("Total nCells filter nCount 0.01 and 0.99 perc: {n.filter.99}")
  cli_alert_info("Total nCells filter perc_mito >{mito.cut}: {n.filter.mito}. ({round(n.filter.mito/n.orig,2)*100}%)")
  cli_alert_info("Total nCells filter perc_haemo >{haemo.cut}: {n.filter.haemo}. ({round(n.filter.haemo/n.orig,2)*100}%)")
  cli_alert_info("Total Filter Both: {n.filter.both}. ({round(n.filter.both/n.orig,2)*100}%)")

  cli_alert_info("==============")
}


for(my.sample in pdata$sample_id){
  .qcCutoffFoo(my.sample)
}


##  ---- Feature Scatters ----

# All functions contain
scCustomize::QC_Plot_UMIvsGene(
  seurat_object = seurat.object,
  meta_gradient_name = "percent_mitochondria",
  meta_gradient_low_cutoff = 10,
  low_cutoff_gene = 500, high_cutoff_gene = 10000,
  low_cutoff_UMI = 1000, high_cutoff_UMI = 50000)

scCustomize::QC_Plot_UMIvsGene(
  seurat_object = seurat.object,
  meta_gradient_name = "percent_haemoglobin",
  meta_gradient_low_cutoff = 5,
  low_cutoff_gene = 500, high_cutoff_gene = 10000,
  low_cutoff_UMI = 1000, high_cutoff_UMI = 50000)


# nUMI vs nGene (logged)
cli_alert_info("plotting to file {.file {file.name}}")
metadata %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, colour=percent_mitochondria)) +
  geom_point(size=0.1) +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 2000) +
  geom_hline(yintercept = 1000) +
  facet_wrap(~sample_id)



for(i in sort(unique(Idents(seurat.object)))){
  p <- QC_Plot_UMIvsGene(
    meta_gradient_name = "percent_mitochondria",
    seurat_object = subset(seurat.object, idents = i),
    low_cutoff_gene = 800, high_cutoff_gene = 5500,
    low_cutoff_UMI = 1000, high_cutoff_UMI = 100000,
    group.by = "orig.ident", raster = T)
  p + ggtitle(label = paste("Sample",i))
}


p <- scCustomize::QC_Plot_GenevsFeature(
  seurat_object = seurat.object,
  feature1 = "percent_mitochondria",
  low_cutoff_gene = 1000,
  high_cutoff_gene = 10000, high_cutoff_feature = 10)
p <- p + ggtitle(label = "All Samples Mitochondrial load vs nGenes") +
ggsave(p, filename = file.path(plot.dir, "03-1-mito-nGenes-all.png"), width=8, height=7)

# Plot individual QC_Plot_GenevsFeature vs Mitochondrial load...
#fname <- file.path(plot.dir, "03-Mito_vs_nGenes_allSamples.pdf")
#file.name <- file.path(params$results.path, fname)
#cli_alert_info("plotting to file {.file {file.name}}")
#pdf(file = file.name, width=8, height = 7)

for(i in sort(unique(Idents(seurat.object)))){
  p <- scCustomize::QC_Plot_GenevsFeature(
    seurat_object = subset(seurat.object, idents = i),
    feature1 = "percent_mitochondria",
    low_cutoff_gene = 1000, high_cutoff_gene = 10000,
    raster = TRUE,
    high_cutoff_feature = 10)
  plot(p + ggtitle(label = paste("Sample",i)))

}
#dev.off()


# ---- Median QCs per Sample group ----
colnames(seurat.object@meta.data)
median_stats <- scCustomize::Median_Stats(
  seurat_object = seurat.object,
  group_by_var = "orig.ident",
  median_var = c("percent_mitochondria","percent_platlets", "percent_haemoglobin",
    "percent_rcas_flex_probes",
    "percent_rcas_flex_probes_pdgfb", "percent_rcas_flex_probes_rfp"))

median_stats %>% knitr::kable("markdown")

.barplot_y(
  median_stats[-15,], group.var = "orig.ident", y.var = "Median_nCount_RNA")

.barplot_y(
  median_stats[-15,], group.var = "orig.ident", y.var = "Median_nFeature_RNA")


my_groups <- c("sample_type", "sex", "date")
seurat.object@meta.data$date <- paste0("date_",seurat.object@meta.data$date)
table(seurat.object@meta.data$date)
seurat.object@meta.data$percent_mito <- seurat.object@meta.data$percent_mitochondria
for(g in my_groups){
  p <- Plot_Median_Genes(seurat_object = seurat.object, group_by = g)
  plot(p)
  p <- Plot_Median_UMIs(seurat_object = seurat.object, group_by = g)
  p
  p <- Plot_Median_Mito(seurat_object = seurat.object, group_by = g)
  p
  p <- Plot_Median_Other(seurat_object = seurat.object, median_var = "percent_haemoglobin", group_by = g)
  p
}



# Sample Sex
#  In primer, no genes on chrY (all female?)
str(fdata)
table(fdata$chromosome_name)
u <- which(fdata$chromosome_name=="Y")
which(fdata$external_gene_name == "Sry")
# chrY_gene = fdata_table$external_gene_name[fdata_table$chromosome_name == "Y"]
#data.filt$pct_chrY = colSums(data.filt@assays$RNA@counts[chrY_gene, ])/colSums(data.filt@assays$RNA@counts)
#. FeatureScatter(data.filt, feature1 = "Xist", feature2 = "pct_chrY")
#
Seurat::VlnPlot(
  seurat.object,  features = c("Ddx3y"), pt.size = 0, ncol = 3, log=T) + NoLegend()




