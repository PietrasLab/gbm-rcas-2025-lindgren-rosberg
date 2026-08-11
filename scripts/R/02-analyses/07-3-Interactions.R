reticulate::use_python(Sys.getenv("GBM_RCAS_PYMINER_PYTHON",
  unset = path.expand("~/mambaforge/envs/gbm-rcas-pyminer/bin/python")), required = TRUE)
# Manually import umap to make sure it's registered
umap <- reticulate::import("umap", delay_load = FALSE)
future::plan("sequential")

# ---- Prepare Env ----
require(qs2)
require(cli)
require(glue)
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
#conflicts_prefer(dplyr::filter)
require(kableExtra)
library(dplyr)
library(patchwork)
set.seed(169)
library(patchwork)
library(CellChat)
library(NMF)
library(future)
future::plan("sequential")  #
NMF::nmf.options(shared.memory = FALSE)
require(ggalluvial)



require(Seurat)
require(dplyr)
require(ggplot2)
require(Matrix)
require(ggalluvial)
require(glue)
require(cli)
set.seed(169999)
# conflicts_prefer(dplyr::filter)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
# remotes::install_github("huayc09/SeuratExtend")
require(SeuratExtend)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

# ---- Create results dir ----
results.dir <- glue("./results/07-3-interactions")
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)



## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features.")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")


# create annotations with AC TE and NT in one cluster
mdata <- seurat.object@meta.data
table(mdata$Level_3ACM)
table(mdata$Level_4ACM)
seurat.object@meta.data <- mdata
colnames(seurat.object@meta.data)
seurat.object <-  seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData(verbose = T)
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)
table(seurat.object[["Level_3ACM"]])




# ---- CellChat: Level_3 ACM ----
# use only tumor samples, remove ambiguous cells
# Merge the AC TE and NT clusters
seurat.sub <- subset(seurat.object, subset = Level_1 != "Ambiguous" & sample_type != "Healthy" )
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
seurat.sub # 44676 samples within 1 assay
annotation <- "Level_3ACM"
Idents(seurat.sub) <- annotation
table(seurat.sub[[annotation]])
# set min cells per cluster to 25
min.cells <- 25
cell.n <- table(seurat.sub[[annotation]])
# get the groups with >100 cells in all tumor data
keep.idents <- names(cell.n)[cell.n > min.cells]
seurat.sub <- subset(seurat.sub, subset = !!sym(annotation) %in% keep.idents)
seurat.sub # 43830 samples within 1 assay
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
table(seurat.sub[[annotation]])


###  ---- Run & Save data ----
### run for both ECM and SC (secreted signaling)
for (db.i in c("Secreted Signaling", "ECM-Receptor")) { # db.i <- "ECM-Receptor"
  cli::cli_alert_info("Running CellChat for {db.i}")
  cellchat <- CellChat::createCellChat(object = seurat.sub, group.by = "ident", assay = "RNA")
  # get database
  CellChatDB <- CellChat::CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  CellChat::showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChat::subsetDB(CellChatDB, search = db.i, key = "annotation")
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database
  # future::plan("multisession", workers = 4) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  options(future.globals.maxSize = 20 * 1024^3)  # Set to 20 GB
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", nboot = 1000)

  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  ## save
  cellchat.object.name = glue("07-3-CellChat_object_{annotation}_", gsub(" ", "", db.i), ".rds")
  cli::cli_alert_info("Saving CellChat object to {cellchat.object.name}")
  dir.create(file.path(results.dir, annotation), showWarnings = FALSE, recursive = TRUE)
  saveRDS(cellchat, file = file.path(results.dir, annotation, cellchat.object.name))
}



# ---- CellChat: Level_4 AC ----
# use only neoplastic groups and Astrocyte Reactive samples, remove ambiguous cells
# Merge the AC TE and NT clusters
seurat.sub <- subset(seurat.object, subset = Level_3ACM %in% c("Astrocyte R","Neoplastic") )
seurat.sub <- subset(seurat.sub, subset = sample_type != c("Healthy") )
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
seurat.sub #  34377 samples within 1 assay
annotation <- "Level_4"
Idents(seurat.sub) <- annotation
table(seurat.sub[[annotation]])
# set min cells per cluster to 25
min.cells <- 25
cell.n <- table(seurat.sub[[annotation]])
# get the groups with >100 cells in all tumor data
keep.idents <- names(cell.n)[cell.n > min.cells]
seurat.sub <- subset(seurat.sub, subset = !!sym(annotation) %in% keep.idents)
seurat.sub # 34377 samples within 1 assay
seurat.sub <-  seurat.sub %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::ScaleData(verbose = T)
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
table(seurat.sub[[annotation]])
Idents(seurat.sub)

###  ---- Run & Save data ----
### run for both ECM and SC (secreted signaling)
for (db.i in c("Secreted Signaling", "ECM-Receptor")) { # db.i <- "Secreted Signalling"
  cli::cli_alert_info("Running CellChat for {db.i}")
  cellchat <- CellChat::createCellChat(object = seurat.sub, group.by = "ident", assay = "RNA")
  # get database
  CellChatDB <- CellChat::CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  CellChat::showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChat::subsetDB(CellChatDB, search = db.i, key = "annotation")
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database

  # future::plan("multisession", workers = 4) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  options(future.globals.maxSize = 20 * 1024^3)  # Set to 20 GB
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", nboot = 1000)

  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  ## save
  cellchat.object.name = glue("07-3-CellChat_object_{annotation}_", gsub(" ", "", db.i), ".rds")
  cli::cli_alert_info("Saving CellChat object to {cellchat.object.name}")
  saveRDS(cellchat, file = file.path(results.dir, cellchat.object.name))
  rm(cellchat)
}





# ---- CellChat: Level_3ACM TP ----
# Primary tumors
seurat.sub <- subset(seurat.object, subset = Level_1 != "Ambiguous" & sample_type != "Healthy")
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
min.cells <- 25
annotation <- "Level_3ACM"
cell.n <- table(seurat.sub[[annotation]])
# get the groups with >100 cells in all tumor data
keep.idents <- names(cell.n)[cell.n > min.cells]
seurat.sub <- subset(seurat.sub, subset = !!sym(annotation) %in% keep.idents)
seurat.sub <- subset(seurat.sub, sample_type == "Primary")
seurat.sub #
seurat.sub <-  seurat.sub %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::ScaleData(verbose = T)
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
annotation <- "Level_3ACM"
Idents(seurat.sub) <- annotation
table(Idents(seurat.sub))
table(seurat.sub[[annotation]])
Idents(seurat.sub)

###  ---- Run & Save data ----
### run for both ECM and SC (secreted signaling)
for (db.i in c("Secreted Signaling", "ECM-Receptor")) { # db.i <- "Secreted Signalling"
  cli::cli_alert_info("Running CellChat for {db.i}")
  cellchat <- CellChat::createCellChat(object = seurat.sub, group.by = "ident", assay = "RNA")
  # get database
  CellChatDB <- CellChat::CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  CellChat::showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChat::subsetDB(CellChatDB, search = db.i, key = "annotation")
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database

  # future::plan("multisession", workers = 4) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  options(future.globals.maxSize = 20 * 1024^3)  # Set to 20 GB
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", nboot = 1000)

  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  cellchat.object.name = glue("07-3-CellChat_object_{annotation}_PrimaryTumors_", gsub(" ", "", db.i), ".rds")

  # cellchat <- .cellchat_plots(
  #   cellchat,
  #   cellchat.object.name = cellchat.object.name,
  #   annotation = annotation,
  #   analysis.id = "{annotation}_PrimaryTumor",
  #   main.dir = results.dir,
  #   cell_types_to_plot = c("Neoplastic", "Astrocyte R")
  # )
  ## save
  cli::cli_alert_info("Saving CellChat object to {cellchat.object.name}")
  saveRDS(cellchat, file = file.path(results.dir, cellchat.object.name))
  rm(cellchat)
}

# ---- CellChat: Level_3ACM TR ----
# Primary tumors
seurat.sub <- subset(seurat.object, subset = Level_1 != "Ambiguous" & sample_type != "Healthy")
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
min.cells <- 25
cell.n <- table(seurat.sub[[annotation]])
# get the groups with >100 cells in all tumor data
keep.idents <- names(cell.n)[cell.n > min.cells]
seurat.sub <- subset(seurat.sub, subset = !!sym(annotation) %in% keep.idents)
seurat.sub <- subset(seurat.sub, sample_type == "Recurrent")
seurat.sub #  22305 samples within 1 assay
seurat.sub <-  seurat.sub %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 5000) %>%
  Seurat::ScaleData(verbose = T)
seurat.sub <- .seuratFactorizeMdata(seurat.sub)
annotation <- "Level_3ACM"
Idents(seurat.sub) <- annotation
table(Idents(seurat.sub))
table(seurat.sub[[annotation]])
Idents(seurat.sub)

###  ---- Run & Save data ----
### run for both ECM and SC (secreted signaling)
for (db.i in c("Secreted Signaling", "ECM-Receptor")) { # db.i <- "Secreted Signalling"
  cli::cli_alert_info("Running CellChat for {db.i}")
  cellchat <- CellChat::createCellChat(object = seurat.sub, group.by = "ident", assay = "RNA")
  # get database
  CellChatDB <- CellChat::CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  CellChat::showDatabaseCategory(CellChatDB)
  cellchat@DB <- CellChat::subsetDB(CellChatDB, search = db.i, key = "annotation")
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database

  # future::plan("multisession", workers = 4) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  options(future.globals.maxSize = 20 * 1024^3)  # Set to 20 GB
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean", nboot = 1000)

  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  ## save
  cellchat.object.name = glue("07-3-CellChat_object_{annotation}_RecurrentTumors_", gsub(" ", "", db.i), ".rds")

  # cellchat <- .cellchat_plots(
  #   cellchat,
  #   cellchat.object.name = cellchat.object.name,
  #   annotation = annotation,
  #   analysis.id = "{annotation}_PrimaryTumor",
  #   main.dir = results.dir,
  #   cell_types_to_plot = c("Neoplastic", "Astrocyte R")
  # )

  cli::cli_alert_info("Saving CellChat object to {cellchat.object.name}")
  saveRDS(cellchat, file = file.path(results.dir, cellchat.object.name))
  rm(cellchat)
}




# ---- LIANA ----
annot <- "Level_3ACM"
seurat.sub <- seurat.object
Idents(seurat.sub) <- annot
table(Idents(seurat.sub))

## ---- Lev3_AC Healthy -----
annotation <- "Level_3ACM"
liana.object <- subset( seurat.object, subset = Level_1 != "Ambiguous" &
    sample_type == "Healthy" )
liana.object # 19396
liana.object <- .seuratFactorizeMdata(liana.object)
Idents(liana.object) <- annotation
Idents(liana.object)
table(Idents(liana.object))

.liana.wrapper(
  analysis.id = "Level3AC_Healthy",
  analysis.main.dir = results.dir,
  seurat.object = liana.object,
  liana.annot = annotation,
  n_varFeatures = 5000,
  liana.downsample.n = 1000,
  liana.min.cells = 15
)

## ---- Lev3_AC Primary -----
annotation <- "Level_3ACM"
liana.object <- subset( seurat.object, subset = Level_1 != "Ambiguous"
  & sample_type == "Primary" )
liana.object # 22343
liana.object <- .seuratFactorizeMdata(liana.object)
Idents(liana.object) <- annotation
table(Idents(liana.object))

.liana.wrapper(
  analysis.id = "Level3AC_Primary",
  analysis.main.dir = results.dir,
  seurat.object = liana.object,
  liana.annot = annotation,
  n_varFeatures = 5000,
  liana.downsample.n = 1000,
  liana.min.cells = 15
)


#### ---- Lev3_AC Recurrent -----
annotation <- "Level_3AC"
liana.object <- subset( seurat.object, subset = Level_1 != "Ambiguous" &
    sample_type == "Recurrent" )
liana.object # 19396
liana.object <- .seuratFactorizeMdata(liana.object)
Idents(liana.object) <- annotation
table(Idents(liana.object))

.liana.wrapper(
  analysis.id = "Level3AC_Recurrent",
  analysis.main.dir = results.dir,
  seurat.object = liana.object,
  liana.annot = annotation,
  n_varFeatures = 5000,
  liana.downsample.n = 1000,
  liana.min.cells = 15
)
gc()
