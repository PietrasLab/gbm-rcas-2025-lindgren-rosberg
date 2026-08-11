
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
library(dplyr)
library(patchwork)

set.seed(169999)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table
results.dir <- "./results/07-deg-findmarkers/"
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)


## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features.")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")



# ---- FindAllMarkers ----

## ---- Level_3: NonNeaoplastic ----
annot <- "Level_1"
seurat.subset <- subset(seurat.object, subset = Level_1 != "Ambiguous")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df[[annot]] <- deg_df[["cluster"]]
deg_df <- .FactorizeMdata(deg_df)


tibble(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-Level_1.csv")
readr::write_csv(deg_df, csv.file)


## ---- Level_3: all ----
annot <- "Level_3"
seurat.subset <- subset(seurat.object, subset = Level_1 != "Ambiguous")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df[[annot]] <- deg_df[["cluster"]]
deg_df <- .FactorizeMdata(deg_df)


tibble(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-Level_3.csv")
readr::write_csv(deg_df, csv.file)


## ---- Level_3: Neoplastic TP vs TR ----
annot <- "sample_type"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  #min.pct = 0.5
  # logfc.threshold = 0.5
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df[[annot]] <- deg_df[["cluster"]]
deg_df <- .FactorizeMdata(deg_df)


tibble(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-SampleType-Neoplastic.csv")
readr::write_csv(deg_df, csv.file)


## ---- Level_3: NonNeaoplastic ----
annot <- "Level_3"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Non-Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df[[annot]] <- deg_df[["cluster"]]
deg_df <- .FactorizeMdata(deg_df)


tibble(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-NonNeoplastic-Level_3.csv")
readr::write_csv(deg_df, csv.file)

## ---- Level_3ACM: NonNeaoplastic ----
# astrocyte TE/NT merged
annot <- "Level_3ACM"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Non-Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df[[annot]] <- deg_df[["cluster"]]
deg_df <- .FactorizeMdata(deg_df)


tibble(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-NonNeoplastic-Level_3ACM.csv")
readr::write_csv(deg_df, csv.file)


## ---- Level_4: Neaoplastic Sybtypes ----
annot <- "Level_4"
seurat.subset <- subset(seurat.object, subset = Level_1 == "Neoplastic")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)

tibble(deg_df)
head(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-Neoplastic-Level_4.csv")
readr::write_csv(deg_df, csv.file)
readr::read_csv(csv.file)


## ---- Level_4: All,  ----
annot <- "Level_4ACM"
seurat.subset <- subset(seurat.object, subset = Level_1 != "Ambiguous")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)

tibble(deg_df)
head(deg_df)
csv.file <- file.path(results.dir, "07-1-FindAllMarkers-Level_4ACM.csv")
readr::write_csv(deg_df, csv.file)
readr::read_csv(csv.file)


## ---- Level_4: All ----
annot <- "Level_4"
seurat.subset <- subset(seurat.object, subset = Level_1 != "Ambiguous")
seurat.subset # 19073 features across 30345 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)
Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)

tibble(deg_df)
head(deg_df)

csv.file <- file.path(results.dir, "07-1-FindAllMarkers-Level_4.csv")
readr::write_csv(deg_df, csv.file)
readr::read_csv(csv.file)


# ---- AC Subtypes ----
## ---- Astrocyte TENTR Markers ----
annot <- "Level_4"
Idents(seurat.object) <- annot
table(Idents(seurat.object))
annotations <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R", "Neopl-ACR") #

seurat.subset <- subset(seurat.object, subset = Level_4 %in% annotations)
seurat.subset # 19073 features across 5144 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)
#. 19073 features across 5553 samples within 1 assay
Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()


deg_df <- FindMarkers(
  seurat.subset,
  ident.1 = "Neopl-ACR",
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)
deg_df$cluster <- "Neopl-ACR"
deg_df$gene <- rownames(deg_df)
deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df <- .FactorizeMdata(deg_df)

readr::write_csv(deg_df,  "./results/07-deg-findmarkers/07-1-FindMarkers-NeoplACR_vs_Astrocyte.csv")

topN <- 10
topGenes <- deg_df %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=50)




# ---- AC Subtypes ----
## ---- Tumor ACR Vs Astrocyte TENTR Markers ----
annot <- "Level_3AC"
annotations <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R") #

seurat.subset <- subset(seurat.object, subset = Level_3AC %in% annotations)
seurat.subset # 19073 features across 5144 samples within 1 assay
seurat.subset <- .seuratFactorizeMdata(seurat.subset)

Idents(seurat.subset) <- annot
table(Idents(seurat.subset))

seurat.subset <- seurat.subset %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()


deg_df <- FindAllMarkers(
  seurat.subset,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1
)

deg_df <- deg_df %>%
  arrange(cluster, p_val_adj)
deg_df <- .FactorizeMdata(deg_df)

readr::write_csv(deg_df,  "./results/07-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv")



## ---- Astrocyte R: Primary vs Tumor----

seurat.sub <- subset(seurat.object, subset = Level_3AC %in%
    c("Astrocyte R"))
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
seurat.sub <- seurat.sub %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()
seurat.sub # 713 samples
Idents(seurat.sub) <- "disease_state"
table(Idents(seurat.sub))

deg_df <- Seurat::FindMarkers(
  seurat.sub,
  ident.1 = "Tumor",
  ident.2 = "Healthy",
  only.pos = TRUE,
  # min.pct = 0.5,
  logfc.threshold = 1
  )
deg_df
write_csv(deg_df, file = file.path(results.dir, glue("07-1-FindAllMarkers-AstroR-TvH.csv")))




