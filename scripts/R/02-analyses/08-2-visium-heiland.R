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
results.dir <- glue("./results/08-2-visium-correlations")
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)


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




# ---- Heiland Visium v1, n=20 ----

### Load spatial from h5,
raw.data.path.main <- "./data/external/10x-visium-heiland/"
sample.paths <- list.dirs(path = raw.data.path.main, recursive = FALSE, full.names = TRUE)
sample.paths <- sample.paths[grepl("#UK", basename(sample.paths))]
# Remove '#' and create clean names
spatial.assay  <- "Spatial"
sample.ids <- basename(sample.paths)

out.dir <- file.path("./results/08-2-Visium-Heiland")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
#out.dir <- results.dir
for (i in seq_along(sample.paths)) { # i <- 1
  cli::cli_h2("Processing Heiland data: sample {i} of {length(sample.paths)}: {sample.ids[i]}")
  sample.path <- sample.paths[i]
  # sample.id <- sample.ids[i]
  sample.id <- basename(sample.path)


  cli::cli_h1("Processing sample: {sample.id}")
  # Set variables for this sample
  raw.data.path <- file.path(sample.path, "outs")
  analysis.name.i <- glue::glue("Visium_{sample.id}")

  #  Load raw data
  h5.file <- list.files(raw.data.path, pattern = "filtered_feature_bc_matrix.h5", full.names = FALSE)
  object <- Seurat::Load10X_Spatial(
    data.dir = raw.data.path,
    filename = h5.file
  )
  object
  head(rownames(object))

  ###  Normalize data
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)

  #module.scores <- lapply(deg.df.list, function(deg.df) {
  object.tmp <- AddModuleScore(
    object = object,
    features = gene_sigs.hs,
    # ctrl = deg.df$ctrl,
    assay = spatial.assay,
    name = glue("ModuleScore__")
  )
  colnames(object.tmp@meta.data)
  module.score.df <- object.tmp@meta.data[, grep(glue("ModuleScore__"), colnames(object.tmp@meta.data) )]
  #module_names <- paste0(glue(names(gene_sigs.hs)))
  module_names <- names(gene_sigs.hs)
  module_names <- gsub(" |-", "_", module_names)
  colnames(module.score.df) <- module_names

  # save module scores
  file.name <- file.path(out.dir, glue("08-2-Heiland-ModuleScores-SigCollection-{sample.id}.rds"))
  saveRDS(module.score.df, file = file.name)

  # Add metadata to object
  object.tmp <- AddMetaData(object.tmp, module.score.df)
  colnames(object.tmp@meta.data)


  # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  require(corrplot)
  library(Hmisc)

  cor_matrix <- cor(module.score.df,  method = "spearman")
  head(cor_matrix)
  # save tje correlation matrix
  file.name <- file.path(out.dir, glue("08-2-Heiland-ModuleScores-Correlation-{sample.id}.rds"))
  saveRDS(cor_matrix, file = file.name)


  # Compute the correlation matrix
  str(cor_matrix)
  pdf.name <- file.path(out.dir, glue("08-2-Heiland-ModuleScores-Correlation-{sample.id}.pdf"))
  pdf(pdf.name, width = 10, height = 11)
  corrplot(cor_matrix, tl.col = "black", title = sample.id,
    method = 'square', diag = FALSE, order = 'hclust',
    addrect = 12, rect.col = 'black', rect.lwd = 3,
    # col = colorRampPalette(c("blue", "white", "red"))(200),
    col = rev(COL2('BrBG')),
    insig = "n"          # blank out insignifica
    # tl.pos = 'd'
  )
  dev.off()


  # Extract color scale
  brbg_colors <- rev(COL2("BrBG"))
  #
  # Prepare long-format data for each target
  plot_cor_for <- function(target_sig) {
    cor_df <- tibble(
      signature = rownames(cor_matrix),
      correlation = cor_matrix[, target_sig]
    ) %>%
      filter(signature != target_sig) %>%
      mutate(signature = fct_reorder(signature, correlation))

    ggplot(cor_df, aes(x = signature, y = correlation, fill = correlation)) +
      geom_col() +
      geom_hline(yintercept = 0, color = "black") +
      scale_fill_gradientn(colors = brbg_colors, limits = c(-1, 1)) +
      coord_cartesian(ylim = c(-1, 1)) +  # clamp y-axis to [-1, 1]
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid.major = element_line(color = "gray80", size = 0.4),
        panel.grid.minor = element_line(color = "gray90", size = 0.2)
      ) +
      labs(
        title = glue("Correlations with {target_sig} in {sample.id} data"),
        x = "Module Score",
        y = "Pearson Correlation",
        fill = "Correlation"
      )
  }
  # Generate and print plots
  target_sigs <- c("Neopl_ACR", "Astrocyte_R")

  pdf.name <- file.path(out.dir, glue("08-2-Heiland-ModuleScores-Correlation-Rankplots-{sample.id}.pdf"))
  pdf(pdf.name, width = 10, height = 7)
  for( ts in target_sigs){
    cli::cli_alert("Plotting correlations for: {ts}")
    p <- plot_cor_for(ts)
    plot(p)
  }
  dev.off()

  #  plot module.scores
  pdf.name <- file.path(out.dir, glue("08-2-Heiland-ModuleScores-SpatialPlot-{sample.id}.pdf"))
  pdf(pdf.name, width = 14, height = 10)
  for ( module_name in module_names) { # module_name <- module_names[1]
    cli::cli_alert("Plotting: {module_name}")
    p <- SpatialFeaturePlot(
      object.tmp,
      features = module_name,
      pt.size.factor = 2.0,
      image.alpha = 0
    )
    p <- p + theme_void() +
      theme(
        panel.background = element_rect(fill = "gray15", color = NA),
        plot.background = element_rect(fill = "gray15", color = NA),
        legend.background = element_rect(fill = "gray15"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, color = "white"),
        legend.text = element_text(size = 8, color = "white"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "white")
      ) +
      ggtitle(glue("Module score: {module_name}"))
    print(p)
  }
  dev.off()

  rm(object.tmp)
} # end loop all samples

# Get correlation files



