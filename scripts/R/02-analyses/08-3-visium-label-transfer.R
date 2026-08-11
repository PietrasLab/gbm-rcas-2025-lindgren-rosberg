# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(glue)

  # Core analysis
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)

  library(SeuratWrappers)
  # library(SeuratData)   # <- enable only if you actually use example datasets
  library(Banksy)
  library(spacexr)

  library(BiocParallel)
  library(gridExtra)
  library(kableExtra)

  # I/O
  library(hdf5r)
  library(arrow)
  library(reticulate)
  library(msigdbr)
  require(purrr)
  require(forcats)

  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(purrr)
  library(homologene)

})

# Resolve common function name conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::first)
options(future.globals.maxSize = 1000 * 1024^2)  # 1 GiB




## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

# ---- Create results dir ----
results.dir <- glue("./results/08-3-visium-label-transfer")
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

## ---- Load scRNA flex seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features.")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)



#---- Load Gene Signatures ----
# richards signatures collection (Nefter and Richards)
# signatures for all Level 4 groups
richards_sigs.hs <- readRDS(file.path("./references/genesets/Richards_NatCancer_2021_GeneSets.rds"))
u <- grep("Neftel|_Richards|InHouse_", names(richards_sigs.hs))
lapply(richards_sigs.hs[u], length)
richards_sigs.hs <- richards_sigs.hs[u]
names(richards_sigs.hs)
richards_sigs.mm <- readRDS(file.path("./references/genesets/Richards_NatCancer_2021_GeneSets_Mouse.rds"))
u <- grep("Neftel|_Richards|InHouse_", names(richards_sigs.mm))
lapply(richards_sigs.mm[u], length)
richards_sigs.mm <- richards_sigs.mm[u]
names(richards_sigs.mm)
hall50 <- list(
  human = msigdbr(species = "Homo sapiens", category = "H"),
  mouse = msigdbr(species = "Mus musculus", category = "H")
)

# rcas gbm sig
annot <- "Level_4"
file.name <- "./results/07-1-deg-findmarkers/07-1-FindAllMarkers-Level_4.csv"
rcas_genes_df <-  readr::read_delim( file = file.name)
sig_rcas_genes <- rcas_genes_df %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

sig_rcas_genes <- sig_rcas_genes %>%
  group_by(cluster) %>%
  summarise(genes = list(gene)) %>%
  tibble::deframe()

# convert to human orthologs
cli::cli_h1("Converting Mouse Genes to Human Orthologs")
str(sig_rcas_genes)
# Convert Entrez IDs to Mouse Gene Symbols
# helper: map a character vector of mouse symbols -> unique human symbols
mouse_vec_to_human_syms <- function(mouse_syms) {
  # returns data.frame with columns mouseGene, humanGene
  hm <- homologene::mouse2human(mouse_syms)
  hm %>%
    transmute(human = humanGene) %>%
    dplyr::filter(!is.na(human), human != "") %>%
    distinct(human) %>%
    pull(human) %>%
    sort()
}


# convert list, preserving names
sig_rcas_genes.hs <- map(sig_rcas_genes, mouse_vec_to_human_syms)
str(lapply(sig_rcas_genes.hs, head))

# Hypoixa hallmark
# hallmark msig pathways. 200 genes
msig.hallmarks.hs <- msigdbr::msigdbr("Homo sapiens", "H") %>%
  SCPA::format_pathways()
names(msig.hallmarks.hs) <- as.character(lapply(msig.hallmarks.hs, function(x)unique( x$Pathway)))
msig.hallmarks.mm <- msigdbr::msigdbr("Mus musculus", "H") %>%
  SCPA::format_pathways()
names(msig.hallmarks.mm) <- as.character(lapply(msig.hallmarks.mm, function(x)unique( x$Pathway)))

## Combine signature
str(richards_sigs.hs)
gene_sigs.hs <- c(richards_sigs.hs, sig_rcas_genes.hs, list("HALLMARK_HYPOXIA"=msig.hallmarks.hs[["HALLMARK_HYPOXIA"]]$Genes))
names(gene_sigs.hs)
gene_sigs.mm <- c(richards_sigs.mm, sig_rcas_genes, list("HALLMARK_HYPOXIA"=msig.hallmarks.mm[["HALLMARK_HYPOXIA"]]$Genes))
names(gene_sigs.hs)






# ---- Label Transfer: TransferAnchors from RNAseq to Visium HD ----
# Aim for higher resolution here, use larger object, 8um bins

## ----  Healthy sample  ----
query.object <- readRDS("./data/processed/seurat/seurat_visium-hd-healthy-8um_v1.0.rds")
analysis.name <- glue::glue("TransferAnchors-visium.hd.h")

### --- Prep Ref and query data objects ----
### Quickfix so that RFP probe is named similaly on both platforms
rownames(query.object)[rownames(query.object) == "CUSTOMPROBE-RFP"] <- "Tg-RFP-1"
# Process query object (visium data set)
query.object <- query.object %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

### ---- Integration with RNA seq ref.object ----
seurat.object # full rnaseq object
### Filter so that it matches sample setup in spatial
ref.object <- seurat.object # rnaseq
ref.object <- subset(ref.object, subset = sample_type %in% c("Healthy"))
ref.object <- subset(ref.object, subset = Level_1 != c("Ambiguous"))
ref.object # 19073 features across 19494 samples within 1 assay
table(ref.object[["Level_3AC"]])

ref.object@meta.data[["TransferLabels"]] <- ref.object@meta.data[["Level_3AC"]] %>%
  forcats::fct_collapse(
    "Endothelial/Mural" = c("Endothelial","Mural"),
    .other_level = NULL  # keep others unchanged
  )
table(ref.object@meta.data[["TransferLabels"]] )
ref.object <- .seuratFactorizeMdata(ref.object)



base::intersect(rownames(ref.object), rownames(query.object))
# subsetlect ref cells to predict in spatial. Use only cells from Healthy brain (and exclude ambient if any)
Idents(ref.object) <- "TransferLabels"
table(Idents(ref.object))
# idents.keep.ref <- 25 # discard if an identity has less than this number of cells (none has here)
# cell_counts <- table(Idents(ref.object))
# idents_to_keep <- names(cell_counts[cell_counts >= idents.keep.ref])
# ref <- subset(ref.object, idents = idents_to_keep)
# table(Idents(ref.object))
ref.object <- ref.object %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
  NormalizeData() %>%
  ScaleData() %>%
  .seuratFactorizeMdata()

anchors <- FindTransferAnchors(
  reference = ref.object,
  query = query.object,
  normalization.method = "LogNormalize",
  npcs = 50
)
table(ref.object[["TransferLabels"]])
predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = ref.object@meta.data[["TransferLabels"]],
  prediction.assay = TRUE,
  weight.reduction = query.object[["pca"]],
  dims = 1:50
)
# save predictions as dataframe
#out.dir <- file.path(results.dir, "TransferAnchors")
#dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
file.name <- file.path(results.dir, glue("08-3-TransferAnchors-healthy-PredictionsAssay.rds"))
saveRDS(predictions.assay, file = file.name)

str(predictions.assay)
predictions.assay <- readRDS(file.name)
# look at predictions
rownames(predictions.assay)
predictions.df <- predictions.assay@data
u <- which(rownames(predictions.assay)=="max")# for each spot the cell type with the highest score is the predicted cell type
predictions.df <- predictions.df[-u,]
table(rownames(predictions.assay)[-u][apply(predictions.df[-u,], 2, which.max)])
str(predictions.df)
rownames(predictions.df)
predictions.df <- as.data.frame(t(predictions.df))
str(predictions.df)

threshold1 <- 0.35
threshold2 <- 0.35

M <- as.matrix(predictions.df)  # rows = spots, cols = labels

# Top-1 per row
max.idx <- apply(M, 1, which.max)
max.val <- M[cbind(seq_len(nrow(M)), max.idx)]
pred1   <- colnames(M)[max.idx]
pred1[max.val < threshold1] <- NA

# Store (align by rownames to be safe)
df <- data.frame(
  barcode = rownames(M),
  predicted_celltype = NA_character_,
  prediction_score = NA_real_,
  stringsAsFactors = FALSE
)
rownames(df) <- df$barcode
df[rownames(M), "prediction_score"]   <- max.val
df[rownames(M), "predicted_celltype"] <- pred1

# Top-2 per row (mask the max first)
M2 <- M
M2[cbind(seq_len(nrow(M2)), max.idx)] <- -Inf
second.idx <- apply(M2, 1, which.max)
second.val <- M2[cbind(seq_len(nrow(M2)), second.idx)]
pred2      <- colnames(M)[second.idx]

keep <- !is.na(df[rownames(M), "predicted_celltype"]) & (second.val >= threshold2)
pred2[!keep] <- NA
df$prediction_score2   <- NA_real_
df$predicted_celltype2 <- NA_character_
df[rownames(M), "prediction_score2"]   <- second.val
df[rownames(M), "predicted_celltype2"] <- pred2

# Combine
df$predicted_combo <- ifelse(
  is.na(df$predicted_celltype2),
  df$predicted_celltype,
  paste(df$predicted_celltype, df$predicted_celltype2, sep = " | ")
)

# Checks
table(is.na(df$predicted_celltype))
table(is.na(df$predicted_celltype2))
table(df$predicted_celltype,  useNA = "ifany")
table(df$predicted_celltype2, useNA = "ifany")
table(df$predicted_combo,    useNA = "ifany")

head(df)
predictions.df$barcode <- rownames(predictions.df)
str(df)
table(df$predicted_celltype)
file.name <- file.path(results.dir, "08-3-TransferAnchors-healthy-Predictions-DataFrame.rds")
saveRDS(df, file = file.name)




## ----  Tumor sample  ----
# query.object <- qs::qread("./data/processed/06-1-visium-hd-tumor-16um.qs")
query.object <- readRDS("./data/processed/seurat/seurat_visium-hd-tumor-8um_v1.0.rds")
analysis.name <- glue::glue("TransferAnchors-visium.hd.tp")
query.object # 19061 features across 140278 samples within 1 assay

### --- Prep Ref and query data objects ----
### Quickfix so that RFP probe is named similaly on both platforms
rownames(query.object)[rownames(query.object) == "CUSTOMPROBE-RFP"] <- "Tg-RFP-1"
# Process query object (visium data set)
query.object <- query.object %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

### ---- Integration with RNA seq ref.object ----
seurat.object # full rnaseq object
### Filter so that it matches sample setup in spatial
ref.object <- seurat.object # rnaseq
#ref.object <- subset(ref.object, subset = sample_type %in% c("Healthy")) # for tumor sample, run with all cell types, since its a whole section with both tumoe and non tumor parts
ref.object <- subset(ref.object, subset = Level_1 != c("Ambiguous"))
ref.object <- subset(ref.object, subset = Level_4ACM != c("Neopl-RNA-low"))
ref.object # 19073 features across 19494 samples within 1 assay
table(ref.object[["Level_4"]])

# For transferlabels, adopt a simpler taxonomy. Merge the Cell cycle clusters

ref.object@meta.data[["TransferLabels"]] <- ref.object@meta.data[["Level_4"]] %>%
  forcats::fct_collapse(
    "Endothelial/Mural" = c("Endothelial","Mural"),
    "Neopl-CC" = c("Neopl-CC-I", "Neopl-CC-II","Neopl-CC-III"),
    "Neopl-OPC-COP"      = c("Neopl-COP", "Neopl-OPC"),
    .other_level = NULL  # keep others unchanged
  )
table(ref.object@meta.data[["TransferLabels"]] )
ref.object <- .seuratFactorizeMdata(ref.object)


base::intersect(rownames(ref.object), rownames(visium.hd.tp))
# subsetlect ref cells to predict in spatial. Use only cells from Healthy brain (and exclude ambient if any)
Idents(ref.object) <- "TransferLabels"

table(Idents(ref.object))
# idents.keep.ref <- 25 # discard if an identity has less than this number of cells (none has here)
# cell_counts <- table(Idents(ref.object))
# idents_to_keep <- names(cell_counts[cell_counts >= idents.keep.ref])
# ref.object <- subset(ref.object, idents = idents_to_keep)
table(Idents(ref.object))
ref.object <- ref.object %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
  NormalizeData() %>%
  ScaleData() %>%
  .seuratFactorizeMdata()

anchors <- FindTransferAnchors(
  reference = ref.object,
  query = query.object,
  normalization.method = "LogNormalize",
  npcs = 50
)
table(ref.object[["TransferLabels"]])
predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = ref.object@meta.data[["TransferLabels"]],
  prediction.assay = TRUE,
  weight.reduction = query.object[["pca"]],
  dims = 1:50
)

# save predictions as dataframe
# out.dir <- file.path(results.dir, "TransferAnchors")
#dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
file.name <- file.path(results.dir, glue("08-3-TransferAnchors-tumor-PredictionsAssay.rds"))
saveRDS(predictions.assay, file = file.name)
predictions.assay <- readRDS(file.name)


# look at predictions
rownames(predictions.assay)
predictions.df <- predictions.assay@data
u <- which(rownames(predictions.assay)=="max")# for each spot the cell type with the highest score is the predicted cell type
predictions.df <- predictions.df[-u,]
table(rownames(predictions.assay)[-u][apply(predictions.df[-u,], 2, which.max)])
str(predictions.df)
rownames(predictions.df)
predictions.df <- as.data.frame(t(predictions.df))
str(predictions.df)



threshold1 <- 0.35
threshold2 <- 0.35

M <- as.matrix(predictions.df)  # rows = spots, cols = labels

# Top-1 per row
max.idx <- apply(M, 1, which.max)
max.val <- M[cbind(seq_len(nrow(M)), max.idx)]
pred1   <- colnames(M)[max.idx]
pred1[max.val < threshold1] <- NA

# Store (align by rownames to be safe)
df <- data.frame(
  barcode = rownames(M),
  predicted_celltype = NA_character_,
  prediction_score = NA_real_,
  stringsAsFactors = FALSE
)
rownames(df) <- df$barcode
df[rownames(M), "prediction_score"]   <- max.val
df[rownames(M), "predicted_celltype"] <- pred1

# Top-2 per row (mask the max first)
M2 <- M
M2[cbind(seq_len(nrow(M2)), max.idx)] <- -Inf
second.idx <- apply(M2, 1, which.max)
second.val <- M2[cbind(seq_len(nrow(M2)), second.idx)]
pred2      <- colnames(M)[second.idx]

keep <- !is.na(df[rownames(M), "predicted_celltype"]) & (second.val >= threshold2)
pred2[!keep] <- NA
df$prediction_score2   <- NA_real_
df$predicted_celltype2 <- NA_character_
df[rownames(M), "prediction_score2"]   <- second.val
df[rownames(M), "predicted_celltype2"] <- pred2

# Combine
df$predicted_combo <- ifelse(
  is.na(df$predicted_celltype2),
  df$predicted_celltype,
  paste(df$predicted_celltype, df$predicted_celltype2, sep = " | ")
)

# Checks
table(is.na(df$predicted_celltype))
table(is.na(df$predicted_celltype2))
table(df$predicted_celltype,  useNA = "ifany")
table(df$predicted_celltype2, useNA = "ifany")
table(df$predicted_combo,    useNA = "ifany")

head(df)
str(df)
file.name <- file.path(results.dir, "08-3-TransferAnchors-tumor-Predictions-DataFrame.rds")
saveRDS(df, file = file.name)





