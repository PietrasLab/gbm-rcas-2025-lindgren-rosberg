
# ---- Prepare Env ----
library(Seurat)
library(SeuratExtend)
library(scCustomize)
library(CellChat)
library(CellMetaVerse)
library(scDblFinder)
library(Matrix)
library(DelayedMatrixStats)
library(BiocParallel)

library(tidyverse)
library(patchwork)
library(gridExtra)
library(kableExtra)
library(glue)
library(cli)
library(qs2)
library(reticulate)
library(NMF)

future::plan("sequential")  #
NMF::nmf.options(shared.memory = FALSE)

set.seed(169)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table

plot.dir <- glue("./manuscript-figures/figure-6")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)


## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features with added metadata")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
# create annotations with AC TE and NT in one cluster

colnames(seurat.object@meta.data)
seurat.object <-  seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData(verbose = T)
colnames(seurat.object@meta.data)
table(seurat.object[["Level_3"]])
table(seurat.object@meta.data$Level_3ACM, seurat.object@meta.data$sample_type)




# ---- Figure 6A. CellChat Astrocyte R & Neoplastic ----


## ---- Barplots Level_3ACM, Astrocyte R. Secreted ----

annot <- "Level_3ACM"
#cellchat <- readRDS("./results/07-4-interactions/CellChat_object_Level_3ACM_SecretedSignaling.rds")
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_3ACM_SecretedSignaling.rds")

# sizes & colors for the identities present
include.idents <- names(table(cellchat@idents))
color.use <- .color_pal[["Level_3AC"]][include.idents]

cell_type <- "Astrocyte R"
cellchat.df <- list()

for (role.i in c("receiver","sender")) {
  if (role.i == "receiver") {
    weights <- cellchat@net$weight[, cell_type]
    counts  <- cellchat@net$count[, cell_type]
    xlab <- "Sender Cell Type"
    title_suffix <- "as Receiver"
  } else if (role.i == "sender") {
    weights <- cellchat@net$weight[cell_type, ]
    counts  <- cellchat@net$count[cell_type, ]
    xlab <- "Receiver Cell Type"
    title_suffix <- "as Sender"
  }
  cellchat.df[[role.i]] <- data.frame(
    CellType = names(weights),
    Weight   = as.numeric(weights),
    Count    = as.numeric(counts),
    Size     = as.numeric(table(cellchat@idents)[names(weights)]),
    Color    = color.use[names(weights)]
  )
}

for (role.i in c("receiver","sender")){
  g <- ggplot(cellchat.df[[role.i]],
    aes(x = reorder(CellType, -Weight), y = Weight, fill = Color)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_identity() +   # use the colors you provided in Color column
    labs(title = glue("Interactions {role.i} ({cell_type}, Weight)"),
      x = xlab, y = "Interaction Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  ggsave(glue("{plot.dir}/Figure-6-CellChat_AstrocyteR_SecretedSignaling_{role.i}.pdf"),
    plot = g, width = 8, height = 6)
  # save the dataframe
  write.csv(cellchat.df[[role.i]],
    file = glue("{plot.dir}/Figure-6-CellChat_AstrocyteR_SecretedSignaling_{role.i}.csv"),
    row.names = FALSE)
  readr::write_csv(g$data,
    glue("{plot.dir}/Figure-6-CellChat_AstrocyteR_SecretedSignaling_{role.i}-sourcedata.csv"))
}



## ---- Barplots Level_3ACM, Neoplastic. Secreted ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_3ACM_SecretedSignaling.rds")

# sizes & colors for the identities present
include.idents <- names(table(cellchat@idents))
color.use <- .color_pal[["Level_3AC"]][include.idents]

cell_type <- "Neoplastic"
cellchat.df <- list()

for (role.i in c("receiver","sender")) {
  if (role.i == "receiver") {
    weights <- cellchat@net$weight[, cell_type]
    counts  <- cellchat@net$count[, cell_type]
    xlab <- "Sender Cell Type"
    title_suffix <- "as Receiver"
  } else if (role.i == "sender") {
    weights <- cellchat@net$weight[cell_type, ]
    counts  <- cellchat@net$count[cell_type, ]
    xlab <- "Receiver Cell Type"
    title_suffix <- "as Sender"
  }
  cellchat.df[[role.i]] <- data.frame(
    CellType = names(weights),
    Weight   = as.numeric(weights),
    Count    = as.numeric(counts),
    Size     = as.numeric(table(cellchat@idents)[names(weights)]),
    Color    = color.use[names(weights)]
  )
}

for (role.i in c("receiver","sender")){
  g <- ggplot(cellchat.df[[role.i]],
    aes(x = reorder(CellType, -Weight), y = Weight, fill = Color)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_identity() +   # use the colors you provided in Color column
    labs(title = glue("Interactions {role.i} ({cell_type}, Weight)"),
      x = xlab, y = "Interaction Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  ggsave(glue("{plot.dir}/Figure-6-CellChat_Neoplastic_SecretedSignaling_{role.i}.pdf"),
    plot = g, width = 8, height = 6)
  # save the dataframe
  write.csv(cellchat.df[[role.i]],
    file = glue("{plot.dir}/Figure-6-CellChat_Neoplastic_SecretedSignaling_{role.i}.csv"),
    row.names = FALSE)
  readr::write_csv(g$data,
    glue("{plot.dir}/Figure-6-CellChat_Neoplastic_SecretedSignaling_{role.i}-sourcedata.csv"))
}


## ---- Barplots Level_4, Neoplastic. Secreted ----
#cellchat <- readRDS("./results/07-4-interactions/CellChat_object_Level_3ACM_SecretedSignaling.rds")
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_SecretedSignaling.rds")
table(cellchat@idents)
# sizes & colors for the identities present
include.idents <- names(table(cellchat@idents))
include.idents %in% names(.color_pal[["Level_4ACM"]])
color.use <- .color_pal[["Level_4ACM"]][include.idents]

cell_type <- "Astrocyte R"
cellchat.df <- list()

for (role.i in c("receiver","sender")) {
  if (role.i == "receiver") {
    weights <- cellchat@net$weight[, cell_type]
    counts  <- cellchat@net$count[, cell_type]
    xlab <- "Sender Cell Type"
    title_suffix <- "as Receiver"
  } else if (role.i == "sender") {
    weights <- cellchat@net$weight[cell_type, ]
    counts  <- cellchat@net$count[cell_type, ]
    xlab <- "Receiver Cell Type"
    title_suffix <- "as Sender"
  }
  cellchat.df[[role.i]] <- data.frame(
    CellType = names(weights),
    Weight   = as.numeric(weights),
    Count    = as.numeric(counts),
    Size     = as.numeric(table(cellchat@idents)[names(weights)]),
    Color    = color.use[names(weights)]
  )
}

for (role.i in c("receiver","sender")){
  g <- ggplot(cellchat.df[[role.i]],
    aes(x = reorder(CellType, -Weight), y = Weight, fill = Color)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_identity() +   # use the colors you provided in Color column
    labs(title = glue("Interactions {role.i} ({cell_type}, Weight)"),
      x = xlab, y = "Interaction Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  ggsave(glue("{plot.dir}/Figure-6-CellChat_ACR-Neoplastic_SecretedSignaling_{role.i}.pdf"),
    plot = g, width = 8, height = 6)
  # save the dataframe
  write.csv(cellchat.df[[role.i]],
    file = glue("{plot.dir}/Figure-6-CellChat_ACR-Neoplastic_SecretedSignaling_{role.i}.csv"),
    row.names = FALSE)
  readr::write_csv(g$data,
    glue("{plot.dir}/Figure-6-CellChat_ACR-Neoplastic_SecretedSignaling_{role.i}-sourcedata.csv"))
}


## ---- Barplots Level_4, Neoplastic. ECM ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_ECM-Receptor.rds")
table(cellchat@idents)
# sizes & colors for the identities present
include.idents <- names(table(cellchat@idents))
color.use <- .color_pal[["Level_4ACM"]][include.idents]

cell_type <- "Astrocyte R"
cellchat.df <- list()

for (role.i in c("receiver","sender")) {
  if (role.i == "receiver") {
    weights <- cellchat@net$weight[, cell_type]
    counts  <- cellchat@net$count[, cell_type]
    xlab <- "Sender Cell Type"
    title_suffix <- "as Receiver"
  } else if (role.i == "sender") {
    weights <- cellchat@net$weight[cell_type, ]
    counts  <- cellchat@net$count[cell_type, ]
    xlab <- "Receiver Cell Type"
    title_suffix <- "as Sender"
  }
  cellchat.df[[role.i]] <- data.frame(
    CellType = names(weights),
    Weight   = as.numeric(weights),
    Count    = as.numeric(counts),
    Size     = as.numeric(table(cellchat@idents)[names(weights)]),
    Color    = color.use[names(weights)]
  )
}

for (role.i in c("receiver","sender")){
  g <- ggplot(cellchat.df[[role.i]],
    aes(x = reorder(CellType, -Weight), y = Weight, fill = Color)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_identity() +   # use the colors you provided in Color column
    labs(title = glue("Interactions {role.i} ({cell_type}, Weight)"),
      x = xlab, y = "Interaction Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  ggsave(glue("{plot.dir}/Figure-6-CellChat_ACR-Neoplastic_ECMSignaling_{role.i}.pdf"),
    plot = g, width = 8, height = 6)
  # save the dataframe
  write.csv(cellchat.df[[role.i]],
    file = glue("{plot.dir}/Figure-6-CellChat_ACR-Neoplastic_ECMSignaling_{role.i}.csv"),
    row.names = FALSE)
  readr::write_csv(g$data,
    glue("{plot.dir}/Figure-6-CellChat_ACR-Neoplastic_ECMSignaling_{role.i}-sourcedata.csv"))
}


## ---- Barplots Level_4, Neoplastic. Secreted ----
#cellchat <- readRDS("./results/07-4-interactions/CellChat_object_Level_3ACM_SecretedSignaling.rds")
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_SecretedSignaling.rds")
table(cellchat@idents)
# sizes & colors for the identities present
include.idents <- names(table(cellchat@idents))
include.idents %in% names(.color_pal[["Level_4ACM"]])
color.use <- .color_pal[["Level_4ACM"]][include.idents]

cell_type <- "Neopl-ACR"
cellchat.df <- list()

for (role.i in c("receiver","sender")) {
  if (role.i == "receiver") {
    weights <- cellchat@net$weight[, cell_type]
    counts  <- cellchat@net$count[, cell_type]
    xlab <- "Sender Cell Type"
    title_suffix <- "as Receiver"
  } else if (role.i == "sender") {
    weights <- cellchat@net$weight[cell_type, ]
    counts  <- cellchat@net$count[cell_type, ]
    xlab <- "Receiver Cell Type"
    title_suffix <- "as Sender"
  }
  cellchat.df[[role.i]] <- data.frame(
    CellType = names(weights),
    Weight   = as.numeric(weights),
    Count    = as.numeric(counts),
    Size     = as.numeric(table(cellchat@idents)[names(weights)]),
    Color    = color.use[names(weights)]
  )
}

for (role.i in c("receiver","sender")){
  g <- ggplot(cellchat.df[[role.i]],
    aes(x = reorder(CellType, -Weight), y = Weight, fill = Color)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_identity() +   # use the colors you provided in Color column
    labs(title = glue("Interactions {role.i} ({cell_type}, Weight)"),
      x = xlab, y = "Interaction Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  ggsave(glue("{plot.dir}/Figure-6-CellChat_Neopl-ACR_SecretedSignaling_{role.i}.pdf"),
    plot = g, width = 8, height = 6)
  # save the dataframe
  write.csv(cellchat.df[[role.i]],
    file = glue("{plot.dir}/Figure-6-CellChat_Neopl-ACR_SecretedSignaling_{role.i}.csv"),
    row.names = FALSE)
  readr::write_csv(g$data,
    glue("{plot.dir}/Figure-6-CellChat_Neopl-ACR_SecretedSignaling_{role.i}-sourcedata.csv"))
}


## ---- Barplots Level_4, Neoplastic. ECM ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_ECM-Receptor.rds")
table(cellchat@idents)
# sizes & colors for the identities present
include.idents <- names(table(cellchat@idents))
color.use <- .color_pal[["Level_4ACM"]][include.idents]

cell_type <- "Neopl-ACR"
cellchat.df <- list()

for (role.i in c("receiver","sender")) {
  if (role.i == "receiver") {
    weights <- cellchat@net$weight[, cell_type]
    counts  <- cellchat@net$count[, cell_type]
    xlab <- "Sender Cell Type"
    title_suffix <- "as Receiver"
  } else if (role.i == "sender") {
    weights <- cellchat@net$weight[cell_type, ]
    counts  <- cellchat@net$count[cell_type, ]
    xlab <- "Receiver Cell Type"
    title_suffix <- "as Sender"
  }
  cellchat.df[[role.i]] <- data.frame(
    CellType = names(weights),
    Weight   = as.numeric(weights),
    Count    = as.numeric(counts),
    Size     = as.numeric(table(cellchat@idents)[names(weights)]),
    Color    = color.use[names(weights)]
  )
}

for (role.i in c("receiver","sender")){
  g <- ggplot(cellchat.df[[role.i]],
    aes(x = reorder(CellType, -Weight), y = Weight, fill = Color)) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_identity() +   # use the colors you provided in Color column
    labs(title = glue("Interactions {role.i} ({cell_type}, Weight)"),
      x = xlab, y = "Interaction Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(g)
  ggsave(glue("{plot.dir}/Figure-6-CellChat_Neopl-ACR_ECMSignaling_{role.i}.pdf"),
    plot = g, width = 8, height = 6)
  # save the dataframe
  write.csv(cellchat.df[[role.i]],
    file = glue("{plot.dir}/Figure-6-CellChat_Neopl-ACR_ECMSignaling_{role.i}.csv"),
    row.names = FALSE)
  readr::write_csv(g$data,
    glue("{plot.dir}/Figure-6-CellChat_Neopl-ACR_ECMSignaling_{role.i}-sourcedata.csv"))
}


# ---- Figure 6B: individual weights ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_SecretedSignaling.rds")

# min.cells <- 25
# cell.n <- table(suerat.sub[[annotation]])

mat <- cellchat@net$weight
include.idents <- names(table(cellchat@idents))
include.idents
color.use <- .color_pal[["Level_4"]][include.idents]
groupSize <- as.numeric(table(cellchat@idents))


pdf(file.path(plot.dir, "Figure-6-AstroR-Neopl-InteractionCircles_individual_Weight_Secreted.pdf"), width = 12, height = 7)
cli::cli_alert_info("Plotting individual interaction weights")
mat <- cellchat@net$weight
str(mat)
#par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:nrow(mat)) {
i <- "Astrocyte R"
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2, color.use = color.use, vertex.weight = groupSize,
  weight.scale = TRUE, arrow.size = 0.75, edge.weight.max = max(mat),
  title.name = rownames(mat)[i])
#}
dev.off()
readr::write_csv(
  data.frame(source = i, target = colnames(mat), weight = as.numeric(mat[i, ]),
             vertex_size = groupSize),
  file.path(plot.dir, "Figure-6-AstroR-Neopl-InteractionCircles_individual_Weight_Secreted-sourcedata.csv")
)


### ---- ECM  ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_ECM-Receptor.rds")
# min.cells <- 25
# cell.n <- table(suerat.sub[[annotation]])

mat <- cellchat@net$weight
include.idents <- names(table(cellchat@idents))
include.idents
color.use <- .color_pal[["Level_4"]][include.idents]
groupSize <- as.numeric(table(cellchat@idents))


pdf(file.path(plot.dir, "Figure-6-AstroR-Neopl-InteractionCircles_individual_Weight_ECM.pdf"), width = 6, height = 7)
cli::cli_alert_info("Plotting individual interaction weights")
mat <- cellchat@net$weight
str(mat)
# par(mfrow = c(4,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
i <- "Astrocyte R"
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2, color.use = color.use, vertex.weight = groupSize,
  weight.scale = TRUE, arrow.size = 0.75, edge.weight.max = max(mat),
  title.name = rownames(mat)[i])
#}
dev.off()
readr::write_csv(
  data.frame(source = i, target = colnames(mat), weight = as.numeric(mat[i, ]),
             vertex_size = groupSize),
  file.path(plot.dir, "Figure-6-AstroR-Neopl-InteractionCircles_individual_Weight_ECM-sourcedata.csv")
)




# ---- Figure 6: individual weights ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_SecretedSignaling.rds")

# min.cells <- 25
# cell.n <- table(suerat.sub[[annotation]])

mat <- cellchat@net$weight
include.idents <- names(table(cellchat@idents))
color.use <- .color_pal[["Level_4"]][include.idents]
groupSize <- as.numeric(table(cellchat@idents))



# network plots
pdf(file.path(plot.dir, "Figure-6-Neopl-ACR-InteractionCircles_individual_Weight_Secreted.pdf"), width = 12, height = 7)
cli::cli_alert_info("Plotting individual interaction weights")
mat <- cellchat@net$weight
str(mat)
#par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:nrow(mat)) {
i <- "Neopl-ACR"
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2, color.use = color.use, vertex.weight = groupSize,
  weight.scale = TRUE, arrow.size = 0.75, edge.weight.max = max(mat),
  title.name = rownames(mat)[i])
#}
dev.off()
readr::write_csv(
  data.frame(source = i, target = colnames(mat), weight = as.numeric(mat[i, ]),
             vertex_size = groupSize),
  file.path(plot.dir, "Figure-6-Neopl-ACR-InteractionCircles_individual_Weight_Secreted-sourcedata.csv")
)



### ---- ECM  ----
cellchat <- readRDS("./results/07-3-interactions/07-3-CellChat_object_Level_4_ECM-Receptor.rds")

mat <- cellchat@net$weight
include.idents <- names(table(cellchat@idents))
color.use <- .color_pal[["Level_4"]][include.idents]
groupSize <- as.numeric(table(cellchat@idents))


pdf(file.path(plot.dir, "Figure-6-Neopl-ACR-InteractionCircles_individual_Weight_ECM.pdf"), width = 6, height = 7)
cli::cli_alert_info("Plotting individual interaction weights")
mat <- cellchat@net$weight
str(mat)
# par(mfrow = c(4,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
i <- "Neopl-ACR"
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2, color.use = color.use, vertex.weight = groupSize,
  weight.scale = TRUE, arrow.size = 0.75, edge.weight.max = max(mat),
  title.name = rownames(mat)[i])
#}
dev.off()
readr::write_csv(
  data.frame(source = i, target = colnames(mat), weight = as.numeric(mat[i, ]),
             vertex_size = groupSize),
  file.path(plot.dir, "Figure-6-Neopl-ACR-InteractionCircles_individual_Weight_ECM-sourcedata.csv")
)




# ---- Figure 6C ----

for (db.i in c("ECM-Receptor","SecretedSignaling")) { # db.i = "SecretedSignaling"

  # cli::cli_alert_info("Loading CellChat object from {.file {file.name}}")
  # cellchat.object.name = paste0("CellChat_object_", gsub(" ", "", db.i), ".rds")

  cellchat_tp <- readRDS(glue("./results/07-3-interactions/07-3-CellChat_object_Level_3ACM_PrimaryTumors_{db.i}.rds"))
  cellchat_tr <- readRDS(glue("./results/07-3-interactions/07-3-CellChat_object_Level_3ACM_RecurrentTumors_{db.i}.rds"))
  gc()


  # Define the common cell groups
  levels(cellchat_tp@idents)
  levels(cellchat_tr@idents)
  common_groups <- intersect(levels(cellchat_tp@idents), levels(cellchat_tr@idents))

  # Subset each CellChat object to include only the common groups
  cellchat_tp <- CellChat::netAnalysis_computeCentrality(cellchat_tp)
  cellchat_tr <- CellChat::netAnalysis_computeCentrality(cellchat_tr)

  object.list <- list(primary=cellchat_tp, recurrent=cellchat_tr)
  cellchat_merged <- CellChat::mergeCellChat(object.list, add.names = names(object.list))

  # cellchat_merged <- CellChat::computeCommunProb(cellchat_merged)
  # cellchat_merged <- CellChat::computeCommunProbPathway(cellchat_merged)
  str(cellchat_merged@meta)
  str(cellchat_merged@idents)
  table(cellchat_merged@idents$joint)
  cellchat_merged@meta$labels
  str(cellchat_merged@netP )
  str(cellchat_merged@netP )

  ### some selected interactions
  pdf.name <- file.path(plot.dir, glue("Figure-6C-CellChat_rankNets_Selected_{db.i}.pdf"))
  pdf(file=pdf.name, width = 7, height = 12)

  p <- CellChat::rankNet(cellchat_merged, mode = "comparison", slot.name = "netP",
    sources.use = "Astrocyte R",
    targets.use = "Neoplastic", stacked = TRUE, do.stat = TRUE,
    title = glue("Source: Astrocyte R, Target: Neoplastic"))
  plot(p)

  dev.off()
  # rankNet $data: per-pathway information-flow contribution (primary vs recurrent)
  readr::write_csv(p$data, file.path(plot.dir, glue("Figure-6C-CellChat_rankNets_Selected_{db.i}-sourcedata.csv")))

}









# ---- Figure 6E: Dotplots Compare TP TR -----
annotation <- "Level_3AC"
liana.object <- subset( seurat.object, subset = Level_3AC %in% c("Astrocyte R", "Neoplastic") &
    sample_type != "Healthy" )
liana.object # 44676
Idents(liana.object) <- annotation
table(Idents(liana.object))

### ---- Load Liana results -----
liana_list <- list(
  healthy = readRDS("./results/07-3-interactions/Level3AC_Healthy/liana_res.rds"),
  primary = readRDS("./results/07-3-interactions/Level3AC_Primary/liana_res.rds"),
  recurrent = readRDS("./results/07-3-interactions/Level3AC_Recurrent/liana_res.rds")
)
table(liana_list[[1]]$logfc$target)

str(liana_list[[2]])

# set IDs for all three lists
for ( x in names(liana_list) ){
  for ( y in names(liana_list[[x]]) ){
    liana_list[[x]][[y]]$ligand_receptor_pair <- paste0(liana_list[[x]][[y]]$ligand, "_", liana_list[[x]][[y]]$receptor)
    liana_list[[x]][[y]]$row_id <- 1:nrow(liana_list[[x]][[y]])
    liana_list[[x]][[y]] <- liana_list[[x]][[y]] %>%
      dplyr::select(row_id, ligand_receptor_pair, everything())
  }
}
liana_list
colnames(liana_list$healthy[[1]])

liana_list_trunc <- lapply( liana_list,
  function(x) {
    y <- x %>% liana::liana_aggregate(.)  %>%
      dplyr::arrange(aggregate_rank) # %>%
    #  dplyr::filter(aggregate_rank <= 0.01)
    #  y$ligand_receptor_pair <- paste0(x$ligand, "_", x$receptor)
    return(y)
  })
colnames(liana_list_trunc[[1]])
head(as.data.frame(liana_list_trunc[[2]]))

#

### ---- Create dotplot_combined.5 ----
# After running dotplot_combined.4
annotation <- "Level_3AC"
liana_in <- liana_list_trunc[c("primary", "recurrent")]
names(liana_in) <- c("Primary", "Recurrent")
analysis.name="Figure-6E-Dotplot-LIANA-TRvTP-diff-ACRvNeop"

liana_aggr_results_list = liana_in
plot.object = liana.object
annotation = "Level_3AC"
source_groups = "Astrocyte R"
target_groups = "Neoplastic"
plot.n.interactions = 25
analysis.name=analysis.name
alternative.plot = TRUE
highlight_differential = TRUE
expression_palette = rev(viridis::plasma(100))


# Define analysis ID
analysis.id <- glue("{analysis.name}_UnionTop{plot.n.interactions}")

# Initialize variables
sample_types <- names(liana_aggr_results_list)
mean_expr_list <- list()

# Add interaction.complex and interaction_label early
liana_aggr_results_list <- lapply(liana_aggr_results_list, function(x) {
  x %>% mutate(
    interaction.complex = paste(ligand.complex, receptor.complex, sep = "_"),
    interaction_label = paste(source, target, interaction.complex, sep = "__")
  )
})

# Prepare interaction label selection
if (highlight_differential) {
  combined_interactions <- liana_aggr_results_list %>%
    purrr::imap_dfr(~ .x %>%
        dplyr::filter(source %in% source_groups, target %in% target_groups) %>%
        mutate(sample_type = .y))

  all_labels <- unique(combined_interactions$interaction_label)

  complete_data <- expand.grid(interaction_label = all_labels, sample_type = sample_types) %>%
    left_join(combined_interactions, by = c("interaction_label", "sample_type")) %>%
    mutate(
      aggregate_rank = ifelse(is.na(aggregate_rank), 1, aggregate_rank),
      source = ifelse(is.na(source), sub("__.*", "", interaction_label), source),
      target = ifelse(is.na(target), sub(".*__(.*?)__.*", "\\1", interaction_label), target),
      interaction.complex = ifelse(is.na(interaction.complex), sub(".*__", "", interaction_label), interaction.complex),
      ligand.complex = ifelse(is.na(ligand.complex), NA, ligand.complex),
      receptor.complex = ifelse(is.na(receptor.complex), NA, receptor.complex),
      neg_log10_rank = -log10(aggregate_rank)
    )

  significant_labels <- complete_data %>%
    dplyr::group_by(interaction_label) %>%
    dplyr::filter(any(aggregate_rank < 0.05)) %>%
    pull(interaction_label) %>%
    unique()

  diff_data <- complete_data %>%
    dplyr::filter(interaction_label %in% significant_labels) %>%
    group_by(interaction_label) %>%
    summarize(
      max_diff = max(neg_log10_rank) - min(neg_log10_rank),
      n_samples = n_distinct(sample_type)
    ) %>%
    dplyr::filter(n_samples >= 2) %>%
    arrange(desc(max_diff)) %>%
    dplyr::slice_head(n = plot.n.interactions)

  liana_selected <- complete_data %>%
    filter(interaction_label %in% diff_data$interaction_label)

  interaction_labels <- unique(liana_selected$interaction_label)

} else {
  extract_top_interactions <- function(data, source_groups, target_groups, top_n = 50) {
    data %>%
      dplyr::filter(source %in% source_groups, target %in% target_groups) %>%
      group_by(source) %>%
      arrange(aggregate_rank) %>%
      dplyr::slice_head(n = top_n) %>%
      ungroup()
  }

  interactions_list <- lapply(names(liana_aggr_results_list), function(name) {
    extract_top_interactions(liana_aggr_results_list[[name]], source_groups, target_groups, plot.n.interactions)
  })
  names(interactions_list) <- names(liana_aggr_results_list)

  interactions_list <- lapply(interactions_list, function(x) {
    x %>% mutate(interaction_label = paste(source, target, interaction.complex, sep = "__"))
  })
  interaction_labels <- unique(unlist(lapply(interactions_list, function(x) x$interaction_label)))
}

# Combine selected interactions
filtered_list <- lapply(names(liana_aggr_results_list), function(name) {
  liana_aggr_results_list[[name]] %>%
    dplyr::filter(source %in% source_groups, target %in% target_groups) %>%
    dplyr::filter(interaction_label %in% interaction_labels) %>%
    mutate(sample_type = name) %>%
    group_by(interaction_label) %>%
    dplyr::slice(1) %>%
    ungroup()
})

liana_selected <- bind_rows(filtered_list) %>%
  mutate(interaction_label = factor(interaction_label, levels = unique(interaction_label))) %>%
  arrange(aggregate_rank)

# Compute cell type sizes
source_sizes <- plot.object@meta.data %>%
  dplyr::filter(!!sym(annotation) %in% source_groups) %>%
  group_by(!!sym(annotation), sample_type) %>%
  summarize(n_cells_source = n(), .groups = "drop") %>%
  dplyr::rename(source = !!sym(annotation))

target_sizes <- plot.object@meta.data %>%
  dplyr::filter(!!sym(annotation) %in% target_groups) %>%
  group_by(!!sym(annotation), sample_type) %>%
  summarize(n_cells_target = n(), .groups = "drop") %>%
  dplyr::rename(target = !!sym(annotation))

# Merge interaction data with cell sizes
liana_plot <- liana_selected %>%
  left_join(source_sizes, by = c("source", "sample_type")) %>%
  left_join(target_sizes, by = c("target", "sample_type")) %>%
  group_by(interaction_label) %>%
  arrange(desc(n_cells_target), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    neg_log10_rank = -log10(aggregate_rank),
    source_order = match(source, source_groups),
    target_order = match(target, target_groups)
  ) %>%
  group_by(interaction_label) %>%
  mutate(max_neg_log10_rank = max(neg_log10_rank, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(source_order, target_order, desc(max_neg_log10_rank), desc(n_cells_target)) %>%
  mutate(interaction_label = factor(interaction_label, levels = unique(interaction_label)))

# Global size range
size_range <- range(c(liana_plot$n_cells_source, liana_plot$n_cells_target), na.rm = TRUE)

# Main dotplot
p_main <- ggplot(liana_plot, aes(x = neg_log10_rank, y = interaction_label, color = sample_type, size = n_cells_target)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = .color_pal[["sample_type"]]) +
  scale_size(range = c(2, 10), limits = size_range, name = "Target group size (cells)") +
  scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 13) +
  labs(
    x = "-log10(Aggregate Rank)",
    y = "Source_Target_Interaction",
    size = "Group size (cells)",
    color = "Sample Type"
  ) +
  theme(
    panel.grid = element_line(color = "grey90"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 30)
  ) +
  ggtitle("Top ranked interactions: ligand → receptor") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

# Expression calculation
ligand_genes <- unique(unlist(strsplit(unique(liana_plot$ligand.complex), "_")))
receptor_genes <- unique(unlist(strsplit(unique(liana_plot$receptor.complex), "_")))
all_genes <- unique(c(ligand_genes, receptor_genes))

for (st in sample_types) {
  cells_st <- colnames(plot.object)[plot.object$sample_type == st]
  seurat_st <- subset(plot.object, cells = cells_st)

  for (src in source_groups) {
    cells_src <- colnames(seurat_st)[seurat_st[[annotation]] == src]
    if (length(cells_src) > 0) {
      mean_expr_list[[paste(st, src, "source", sep = "_")]] <- colMeans(FetchData(seurat_st, vars = all_genes, cells = cells_src, layer = "data"))
    }
  }

  for (tgt in target_groups) {
    cells_tgt <- colnames(seurat_st)[seurat_st[[annotation]] == tgt]
    if (length(cells_tgt) > 0) {
      mean_expr_list[[paste(st, tgt, "target", sep = "_")]] <- colMeans(FetchData(seurat_st, vars = all_genes, cells = cells_tgt, layer = "data"))
    }
  }
}

# Expression range
all_expr_values <- unlist(lapply(mean_expr_list, function(x) unname(unlist(x))))
expr_range <- range(all_expr_values, na.rm = TRUE)


# ---- Extract returned objects ----
# liana_plot <- result$liana_plot
# mean_expr_list <- result$mean_expr_list
# expr_range <- result$expr_range

# ---- Recreate p_source_target ----
p_source_target <- ggplot2::ggplot(liana_plot) +
  geom_tile(aes(x = -5, y = as.numeric(interaction_label), fill = source), width = 2.5, height = 0.9) +
  geom_tile(aes(x = 5, y = as.numeric(interaction_label), fill = target), width = 2.5, height = 0.9) +
  geom_segment(aes(x = -2, xend = 2, y = as.numeric(interaction_label), yend = as.numeric(interaction_label)),
    arrow = arrow(length = unit(0.2, "cm")), color = "gray") +
  geom_text(aes(x = -6.5, y = as.numeric(interaction_label), label = ligand.complex), hjust = 1.2, size = 3) +
  geom_text(aes(x = 6.5, y = as.numeric(interaction_label), label = receptor.complex), hjust = 0, size = 3) +
  scale_fill_manual(values = .color_pal[[annotation]]) +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_reverse() +
  theme_void() +
  ggtitle("Ligand → Receptor") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key = element_rect(fill = "white", color = "black"),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 40)
  )

# ---- Generate per-sample source-target plots ----
sample_types <- unique(liana_plot$sample_type)

per_sample_plots <- .liana_generate_all_per_sample_plots(
  sample_types = sample_types,
  liana_plot = liana_plot,
  mean_expr_list = mean_expr_list,
  expr_range = expr_range,
  source_groups = source_groups,
  target_groups = target_groups
)
per_sample_plots
# ---- Final combined figure ----
# all_plots <- c(list(result$combined_plot, p_source_target), per_sample_plots)
all_plots <- c(list(p_main, p_source_target), per_sample_plots)
widths_vector <- c(3, 1, rep(1, length(per_sample_plots)))

final_combined_plot <- patchwork::wrap_plots(all_plots, ncol = length(all_plots)) +
  patchwork::plot_layout(widths = widths_vector, guides = "collect") +
  patchwork::plot_annotation(
    title = glue("MyAnalysis_UnionTop25"),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  )

# ---- Display the final figure ----
final_combined_plot
pdf.name <- file.path(plot.dir, glue("{analysis.name}.pdf"))
pdf(file = pdf.name, width = 22, height= 12)
print(final_combined_plot)
dev.off()
# Source data: the ranked ligand-receptor interaction table underlying every panel
# (p_main dots, source/target strip, per-sample expression), incl. rank + cell sizes.
readr::write_csv(
  liana_plot %>% dplyr::select(dplyr::any_of(c(
    "interaction_label", "source", "target", "ligand.complex", "receptor.complex",
    "sample_type", "aggregate_rank", "neg_log10_rank",
    "n_cells_source", "n_cells_target"))),
  file.path(plot.dir, glue("{analysis.name}-sourcedata.csv"))
)

