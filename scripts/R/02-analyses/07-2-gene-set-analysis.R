
# ---- Prepare Env ----
require(SeuratExtend)
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
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::first)
conflicts_prefer(matrixStats::colRanks)
set.seed(169)
## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table
results.dir <- "./results/07-2-gene-set-analysis/"
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)


## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features.")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")

seurat.object <-  seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData(verbose = T)
colnames(seurat.object@meta.data)

# ---- Load Gene Signature Lists ----
richards_sigs <- readRDS(file.path("./references/genesets/Richards_NatCancer_2021_GeneSets_Mouse.rds"))
# hall50 <- list(
#   human = msigdbr(species = "Homo sapiens", category = "H"),
#   mouse = msigdbr(species = "Mus musculus", category = "H")
# )



# ---- SeuratExtend::GeneSetAnalysis ----
## ---- Astrocyte subsets (TE NT R)  ----
annot <- "Level_3AC"
Idents(seurat.object) <- annot
table(Idents(seurat.object))

seurat.sub <- subset(seurat.object, subset = Level_3AC %in%
    c("Astrocyte TE","Astrocyte NT", "Astrocyte R"))
seurat.sub #  5110 samples within 1 assay

seurat.sub <- seurat.sub %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()
table(Idents(seurat.sub))

##  ---- GO HALLMARKS  Astrocyte reactive Tum vs Healthy  ----
seurat.sub <- SeuratExtend::GeneSetAnalysis(seurat.sub, genesets = hall50$mouse)
matr <- seurat.sub@misc$AUCell$genesets
str(matr)
saveRDS(matr, file.path(results.dir, "07-2-GeneSetAnalysis-AUCell-Hall50-Astrocytes.rds"))



## ---- Astrocyte R (Tumor vs Healthy)  ----
annot <- "Level_3AC"
Idents(seurat.object) <- annot
table(Idents(seurat.object))
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
# Healthy   Tumor
# 63     650
seurat.sub <- SeuratExtend::GeneSetAnalysis(seurat.sub, genesets = hall50$mouse)
matr <- seurat.sub@misc$AUCell$genesets
str(matr)
saveRDS(matr, file.path(results.dir, "07-2-GeneSetAnalysis-AUCell-Hall50-AstrocyteR.rds"))


## ---- Neoplastic cells   ----
seurat.subset <- subset(seurat.object, subset = Level_1 %in% c("Neoplastic"))
annot <- "Level_4"
Idents(seurat.subset) <- annot
table(Idents(seurat.subset))


### ---- Neftel and Richards  ----
u <- grep("Neftel|InHouse_BulkRNAseq_2019", names(richards_sigs))
names(richards_sigs)[u]

seurat.subset <- SeuratExtend::GeneSetAnalysis(seurat.subset, genesets = richards_sigs[u])
matr <- seurat.subset@misc$AUCell$genesets
str(matr)
rownames(matr) <- gsub("\\.", "", rownames(matr))
rownames(matr) <- gsub("InHouse_BulkRNAseq_2019", "Richards_2021", rownames(matr))
saveRDS(matr, file.path(results.dir, "07-2-GeneSetAnalysis-AUCell-NeftelRichads-Neoplastic.rds"))



## ---- Astro/Neopl Astro vs Rest: Richards signatures ----


### ----- Astro/Astro-neoplastic focused annotation ------

annot <- "Level_4AC"
pathways.object <- seurat.object
pathways.object <- subset(pathways.object, subset = Level_1 != "Ambiguous")
pathways.object <- subset(pathways.object, subset = Level_2 %in% c("Neural","Neoplastic"))
pathways.object <-  pathways.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = Inf)
colnames(pathways.object@meta.data)
pathways.object <- .seuratFactorizeMdata(pathways.object)
pathways.object <- SetIdent(pathways.object, value = annot)
table(Idents(pathways.object))

# make a Astro/Astro-neoplastic focused annotation
astro.vec <- c("Astrocyte TE","Astrocyte NT","Astrocyte R")
mdata <- pathways.object@meta.data %>%
  mutate(
    annot = Level_3AC,                               # start with Level_3AC
    # merge TE + NT before other recoding
    annot = if_else(annot %in% c("Astrocyte TE","Astrocyte NT"),
      "Astrocyte TE/NT", annot),
    # recode non-neoplastic cells not in astro.vec → Non-Neoplastic (other)
    annot = case_when(
      Level_1 == "Non-Neoplastic" & !annot %in% c("Astrocyte TE/NT","Astrocyte R")
      ~ "Non-Neoplastic (other)",
      TRUE ~ annot
    )
  ) %>%
  mutate(
    annot = if_else(Level_4 == "Neopl-ACR", "Neoplastic ACR", annot)
  ) %>%
  # rename Neoplastic to Neoplastic (other)
  mutate(
    annot = if_else(annot == "Neoplastic",
      "Neoplastic (other)", annot)
  )
pathways.object$annot <- mdata$annot
table(pathways.object$annot)
# turn into factor with desired ordering
base_levels <- c("Non-Neoplastic (other)", "Astrocyte TE/NT","Astrocyte R","Neoplastic ACR","Neoplastic (other)")
other_levels <- base::setdiff(sort(unique(pathways.object$annot)), base_levels)
pathways.object$annot <- factor(pathways.object$annot, levels = c(base_levels, other_levels))
Idents(pathways.object) <- "annot"
table(Idents(pathways.object))

# Non-Neoplastic (other)        Astrocyte TE/NT            Astrocyte R
# 25235                   4397                    713
# Neoplastic ACR     Neoplastic (other)
# 443                  33284

# Non-Neoplastic (other)        Astrocyte TE/NT            Astrocyte R         Neoplastic ACR
# 6580                   4397                    713                    443
# Neoplastic (other)
# 33284




### ---- Run geneset analysis ----
richards_sigs <- readRDS(file.path("~/Resources/Genesets/RichardsWhitley_GeneSets_Mouse.rds"))
# https://www.nature.com/articles/s43018-020-00154-9
# Supp table 6. https://static-content.springer.com/esm/art%3A10.1038%2Fs43018-020-00154-9/MediaObjects/43018_2020_154_MOESM8_ESM.xlsx
# Curated neural cell and GBM gene signatures. FirstAuthor_Journal_PublicationYear_SignatureName’.
# remove redundant (duplicated) Richards Injury/Differentiation signatures
names(richards_sigs)

### --- select richard sigs for astro/neoplastic analysis ----
richards_is_neural_neoplastic <- c(
  "Cahoy_JNeurosci_2008_astrocyte"               = TRUE,
  "Cahoy_JNeurosci_2008_oligodendrocyte"         = TRUE,
  "Cahoy_JNeurosci_2009_neuron"                  = TRUE,
  "Cahoy_JNeurosci_2010_astro_young"             = TRUE,
  "Cahoy_JNeurosci_2011_astro_mature"            = TRUE,
  "Cahoy_JNeurosci_2012_OPC"                     = TRUE,
  "Cahoy_JNeurosci_2013_OL_myel"                 = TRUE,
  "Cahoy_JNeurosci_2014_MOG_pos"                 = TRUE,
  "Cahoy_JNeurosci_2015_astro_in_vivo"           = TRUE,
  "Cahoy_JNeurosci_2016_cultured_astroglia"      = TRUE,
  "InHouse_BulkRNAseq_2019_DevelopmentalGSC"     = TRUE,
  "InHouse_BulkRNAseq_2019_InjuryResponseGSC"    = TRUE,
  "Liddelow_Nature_2017_A1_reactive_astrocytes"  = TRUE,
  "Liddelow_Nature_2017_A2_reactive_astrocytes"  = TRUE,
  "Neftel_Cell_2019_AC"                          = TRUE,
  "Neftel_Cell_2019_G1.S"                        = TRUE,
  "Neftel_Cell_2019_G2.M"                        = TRUE,
  "Neftel_Cell_2019_MES1"                        = TRUE,
  "Neftel_Cell_2019_MES2"                        = TRUE,
  "Neftel_Cell_2019_NPC1"                        = TRUE,
  "Neftel_Cell_2019_NPC2"                        = TRUE,
  "Neftel_Cell_2019_OPC"                         = TRUE,
  "Nowakowski_Science_2017_Astrocyte_upreg"      = TRUE,
  "Nowakowski_Science_2017_Choroid_upreg"        = FALSE,
  "Nowakowski_Science_2017_EN-PFC1_upreg"        = TRUE,
  "Nowakowski_Science_2017_EN-PFC2_upreg"        = TRUE,
  "Nowakowski_Science_2017_EN-PFC3_upreg"        = TRUE,
  "Nowakowski_Science_2017_EN-V1-1_upreg"        = TRUE,
  "Nowakowski_Science_2017_EN-V1-2_upreg"        = TRUE,
  "Nowakowski_Science_2017_EN-V1-3_upreg"        = TRUE,
  "Nowakowski_Science_2017_Endothelial_upreg"    = FALSE,
  "Nowakowski_Science_2017_Glyc_upreg"           = FALSE,
  "Nowakowski_Science_2017_IN-CTX-CGE1_upreg"    = TRUE,
  "Nowakowski_Science_2017_IN-CTX-CGE2_upreg"    = TRUE,
  "Nowakowski_Science_2017_IN-CTX-MGE1_upreg"    = TRUE,
  "Nowakowski_Science_2017_IN-CTX-MGE2_upreg"    = TRUE,
  "Nowakowski_Science_2017_IN-STR_upreg"         = TRUE,
  "Nowakowski_Science_2017_IPC-div1_upreg"       = TRUE,
  "Nowakowski_Science_2017_IPC-div2_upreg"       = TRUE,
  "Nowakowski_Science_2017_IPC-nEN1_upreg...21"  = TRUE,
  "Nowakowski_Science_2017_IPC-nEN1_upreg...22"  = TRUE,
  "Nowakowski_Science_2017_IPC-nEN3_upreg"       = TRUE,
  "Nowakowski_Science_2017_MGE-IPC1_upreg"       = TRUE,
  "Nowakowski_Science_2017_MGE-IPC2_upreg"       = TRUE,
  "Nowakowski_Science_2017_MGE-IPC3_upreg"       = TRUE,
  "Nowakowski_Science_2017_MGE-RG1_upreg"        = TRUE,
  "Nowakowski_Science_2017_MGE-RG2_upreg"        = TRUE,
  "Nowakowski_Science_2017_MGE-div_upreg"        = TRUE,
  "Nowakowski_Science_2017_Microglia_upreg"      = FALSE,
  "Nowakowski_Science_2017_Mural_upreg"          = FALSE,
  "Nowakowski_Science_2017_OPC_upreg"            = TRUE,
  "Nowakowski_Science_2017_RG-div1_upreg"        = TRUE,
  "Nowakowski_Science_2017_RG-div2_upreg"        = TRUE,
  "Nowakowski_Science_2017_RG-early_upreg"       = TRUE,
  "Nowakowski_Science_2017_U1_upreg"             = FALSE,
  "Nowakowski_Science_2017_U2_upreg"             = FALSE,
  "Nowakowski_Science_2017_U3_upreg"             = FALSE,
  "Nowakowski_Science_2017_U4_upreg"             = FALSE,
  "Nowakowski_Science_2017_nEN-early1_upreg"     = TRUE,
  "Nowakowski_Science_2017_nEN-early2_upreg"     = TRUE,
  "Nowakowski_Science_2017_nEN-late_upreg"       = TRUE,
  "Nowakowski_Science_2017_nIN1_upreg"           = TRUE,
  "Nowakowski_Science_2017_nIN2_upreg"           = TRUE,
  "Nowakowski_Science_2017_nIN3_upreg"           = TRUE,
  "Nowakowski_Science_2017_nIN4_upreg"           = TRUE,
  "Nowakowski_Science_2017_nIN5_upreg"           = TRUE,
  "Nowakowski_Science_2017_oRG_upreg"            = TRUE,
  "Nowakowski_Science_2017_tRG_upreg"            = TRUE,
  "Nowakowski_Science_2017_vRG_upreg"            = TRUE,
  "Verhaak_CancerCell_2010_Classical"            = TRUE,
  "Verhaak_CancerCell_2010_Mesenchymal"          = TRUE,
  "Verhaak_CancerCell_2010_Neural"               = TRUE,
  "Verhaak_CancerCell_2010_Proneural"            = TRUE,
  "Zhong_Nature_2018_Astrocytes_upreg"           = TRUE,
  "Zhong_Nature_2018_Excitatory_neurons_upreg"   = TRUE,
  "Zhong_Nature_2018_Interneurons_upreg"         = TRUE,
  "Zhong_Nature_2018_Microglia_upreg"            = TRUE,
  "Zhong_Nature_2018_OPC_upreg"                  = TRUE,
  "Developmental_Richards"                       = TRUE,
  "Injury_Response_Richards"                     = TRUE
)
all(names(richards_is_neural_neoplastic) %in% names(richards_sigs))
richards_selected <- names(richards_is_neural_neoplastic)[unlist(richards_is_neural_neoplastic)]
richards_sigs <- richards_sigs[richards_selected]
length(richards_sigs)
require(SeuratExtend)

pathways.object <- pathways.object %>%
  droplevels(.) %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::NormalizeData() %>%
  Seurat::ScaleData()

Idents(pathways.object) <- factor(
  pathways.object$annot,
  levels = c(
    "Non-Neoplastic (other)",
    "Astrocyte TE/NT",
    "Astrocyte R",
    "Neoplastic ACR",
    "Neoplastic (other)"
  )
)

# Run geneset analysis
conflicts_prefer(DelayedMatrixStats::colRanks)
pathways.object <- SeuratExtend::GeneSetAnalysis(pathways.object, genesets = richards_sigs)
matr <- pathways.object@misc$AUCell$genesets
saveRDS(matr, file.path(results.dir, "7-2-GeneSetAnalysis-AUCell-Richards-AC-NACR.rds"))


