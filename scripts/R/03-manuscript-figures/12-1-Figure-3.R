

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
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::first)
set.seed(169)


## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
# save all DEGs for Supp Table

plot.dir <- glue("./manuscript-figures/figure-3")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

# ---- Prepare objects ----

## ---- Load seurat object -----
cli::cli_alert_info("Loading UMI-downsampled and cell QC-Filtered Seurat object with 64804 cells and 19073 features with added metadata")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")
colnames(seurat.object@meta.data)
seurat.object <- .seuratFactorizeMdata(seurat.object)
# add level annotations

seurat.object <-  seurat.object %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures(nfeatures = 2000) %>%
  Seurat::ScaleData(verbose = T)
colnames(seurat.object@meta.data)



# ---- Figure 3A UMAP ------
annot <- "Level_3AC"
reduction.name <- "umap_AllCells"
fig.name <- glue("Figure-3-umap-astrocyte-subtypes")
plot.object <- seurat.object
Idents(seurat.object) <- annot
table(seurat.object@meta.data[[annot]])
Idents(seurat.object) <- annot
annotations <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R") # , "Tumor ACR-sig"

plot.object@meta.data[[annot]] <- factor(
  if_else(plot.object@meta.data[[annot]] %in% annotations,
    plot.object@meta.data[[annot]],
    "other"),
  levels = c(annotations, "other")  # Preserve the order of levels
)
Idents(plot.object) <- annot
table(Idents(plot.object))

my.pal <- .color_pal[[annot]][annotations]
my.pal[["other"]] <-  "grey"

p <- DimPlot(pt.size = 1, alpha = 1,shuffle = T,
  seurat.object,
  reduction = reduction.name,
  group.by = annot,
  cols = my.pal,
  label = TRUE,
  raster = FALSE
)
p <- p + ggtitle(glue("{fig.name}"))
p

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 13, height = 9)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))


# ---- GFAP expression UMAP ----
# Extract UMAP coordinates and expression values
gene <- "Gfap"
fig.name <- glue("Figure-3-umap-gfap")

p <- .seuratFeaturePlotHexbin(
  seurat.object = plot.object,
  feature = gene,
  reduction = reduction.name
)

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
p <- p + ggtitle(glue("{fig.name}"))
p
ggsave(pdf.name, p, width = 12, height = 9)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))



# ---- Astrocyte PieCharts   -----
annot <- c( "Level_3AC")
fig.name <- glue("Figure-3-pieChart-astroSubtypes")
# subset <- "Neoplastic"
# fig.name <- glue("Barplot_{annot}_{subset}")

# Process Metadata
mdata.plot <- plot.object@meta.data %>%
  dplyr::filter(!!sym(annot) != "other")
# dplyr::filter(sample_type %in% c("Primary", "Recurrent")) %>%



p1 <- .piechart_plot(
  my.pal = my.pal,
  plot_df = mdata.plot,
  group.var = "sample_type",
  color.var = annot,
)
p1
p1 <- p1 + ggtitle(glue("{fig.name}"))

pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p1, width = 9, height = 5)
readr::write_csv(p1$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))




# ---- DotPlots AC markers ----

## ----  Pan AC markers ----
deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-NonNeoplastic-Level_3.csv")

# get top 10
topN <- 10
topGenes.pan <- deg_df %>%
  # dplyr::rename(Level_3 = cluster) %>%
  dplyr::filter(Level_3 == "Astrocyte") %>%
  .FactorizeMdata(.) %>%
  arrange(Level_3) %>%
  group_by(Level_3) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes.pan, n=50)



## ----  AC TEvNTvR markers ----
deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv")

# get top 10
topN <- 10
topGenes <- deg_df %>%
  dplyr::rename(Level_3AC = cluster) %>%
  .FactorizeMdata(.) %>%
  arrange(Level_3AC) %>%
  group_by(Level_3AC) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=50)

topGenes$gene %in% topGenes.pan$gene



p <- SCpubr::do_DotPlot(cluster = F,
  diverging.direction = 1,
  diverging.palette = "Blues",
  zscore.data = F,
  group.by = "Level_3AC", flip = F,
  seurat.object %>% subset(Level_3 == "Astrocyte"),
  features = c(topGenes.pan$gene, topGenes$gene))

p +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(
    expand = expansion(add = 0.5)  # shifts tick marks relative to grid centers
  )
p
fig.name <- file.path("Figure-3-dotPlot-findMarkers-astro-top10-alt")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 11, height = 3.3)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))



sig.list <- c(
  setNames(split(topGenes.pan, topGenes.pan$Level_3),
    paste0("pan_", names(split(topGenes.pan, topGenes.pan$Level_3)))),
  setNames(split(topGenes, topGenes$Level_3AC),
    paste0("AC_", names(split(topGenes, topGenes$Level_3AC))))
)
names(sig.list)
sig.list <- lapply(sig.list, function(df) {
  df$gene
})
#module.scores <- lapply(deg.df.list, function(deg.df) {
seurat.sub <- seurat.object %>% subset(Level_3 == "Astrocyte")
object.tmp <- AddModuleScore(
  object = seurat.sub,
  features = sig.list,
  # ctrl = deg.df$ctrl,
  # assay = "Spatial.016um",
  name = glue("ModuleScore__")
)

colnames(object.tmp@meta.data)
module.score.df <- object.tmp@meta.data[, grep(glue("ModuleScore__"), colnames(object.tmp@meta.data) )]
module_names <- paste0(glue("ModuleScore__"), names(sig.list))
module_names <- gsub(" |-", "_", module_names)
colnames(module.score.df) <- module_names
saveRDS(module.score.df, file = file.path(plot.dir, "Figure-3-ModuleScores-AstroPanTENTR.rds"))

for(col.i in module_names){
  object.tmp[[col.i]] <- module.score.df[, col.i]
}

object.acr <- object.tmp %>% subset(Level_3AC == "Astrocyte R")

p <- Seurat::VlnPlot(
  fill.by = "ident",
  flip=T,
  object.acr, pt.size = 0,
  stack = F,
  cols = .color_pal[["sample_type"]],
  features = "ModuleScore__AC_Astrocyte_R",
  group.by = "sample_type"
)
print(p)
p <- SCpubr::do_BoxPlot(
  object.acr,
  colors.use = .color_pal[["sample_type"]],
  feature = "ModuleScore__AC_Astrocyte_R",
  group.by = "sample_type"
)
print(p)

library(ggpubr)

## ---- Stats: Kruskal-Wallis omnibus + BH-corrected pairwise Wilcoxon ----
## P-values are reported on every stats panel; STAR Methods needs
## them in source data. Brackets show the BH-ADJUSTED pairwise Wilcoxon p (not the
## ggpubr `stat_compare_means(comparisons=)` default, which is UNADJUSTED). A tidy
## companion `-stats.csv` records the omnibus KW test and each pairwise contrast
## (W statistic, raw p, BH-adjusted p, significance code).
sm.score <- object.acr@meta.data$ModuleScore__AC_Astrocyte_R
sm.group <- object.acr@meta.data$sample_type
sm.n     <- table(sm.group)

kw <- kruskal.test(sm.score ~ sm.group)

.pcode <- function(p) cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
  labels = c("****", "***", "**", "*", "ns"))

pairs <- list(c("Healthy", "Primary"), c("Healthy", "Recurrent"), c("Primary", "Recurrent"))
pw.df <- purrr::map_dfr(pairs, function(pr) {
  wt <- wilcox.test(sm.score[sm.group == pr[1]], sm.score[sm.group == pr[2]])
  tibble::tibble(group1 = pr[1], group2 = pr[2],
    n1 = as.integer(sm.n[pr[1]]), n2 = as.integer(sm.n[pr[2]]),
    statistic = base::unname(wt$statistic), p = wt$p.value)
}) %>% dplyr::mutate(p.adj = p.adjust(p, method = "BH"),
  p.signif = as.character(.pcode(p.adj)))

## bracket y-positions: adjacent contrasts low, the Healthy-vs-Recurrent span highest
sm.rng <- range(sm.score, na.rm = TRUE); sm.step <- 0.12 * diff(sm.rng)
pw.df$y.position <- max(sm.rng) + sm.step * c(1, 3, 2)

## companion stats CSV (omnibus KW row + one row per BH-adjusted pairwise contrast)
stats.out <- dplyr::bind_rows(
  tibble::tibble(test = "Kruskal-Wallis", group1 = "all groups", group2 = NA_character_,
    n1 = length(sm.score), n2 = NA_integer_, statistic = base::unname(kw$statistic),
    df = base::unname(kw$parameter), p = kw$p.value, p.adj = kw$p.value,
    p.adjust.method = "none", p.signif = as.character(.pcode(kw$p.value)),
    method = "Kruskal-Wallis rank sum"),
  pw.df %>% dplyr::transmute(test = "pairwise Wilcoxon", group1, group2, n1, n2,
    statistic, df = NA_real_, p, p.adj, p.adjust.method = "BH", p.signif,
    method = "Wilcoxon rank sum")
)

p <- p +
  ggpubr::stat_pvalue_manual(pw.df, label = "p.signif", tip.length = 0.01) +
  annotate("text", x = 1, y = max(sm.rng) + sm.step * 4, hjust = 0, size = 3,
    label = sprintf("Kruskal-Wallis p = %.2g", kw$p.value))

fig.name <- file.path("Figure-3-ModuleScore-Astrocyte-R-Signature")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 5, height = 7)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))
readr::write_csv(stats.out, file.path(plot.dir, glue("{fig.name}-stats.csv")))



# ---- Visium HD Astro Signatures ----
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
require(glue)
require(cli)
require(hdf5r)
# pak::pkg_install("arrow")
require(arrow)
?SeuratData
?Seurat::Load10X_Spatial
library(SeuratWrappers)
library(Banksy)
options(future.globals.maxSize = 1000 * 1024^2)  # Set limit to 1000 MiB (1 GiB),


if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(spacexr)



# ---- Load data (Visium HD healthy and tumor) ----
data.processed.paths = list(
  "visium-hd-healthy" = "./data/processed/seurat/seurat_visium-hd-healthy-16um_v1.0.rds",
  "visium-hd-tumor" = "./data/processed/seurat/seurat_visium-hd-tumor-16um_v1.0.rds"
)

sd.list <- list()
def.assay <- "Spatial.016um"
for (data.name in names(data.processed.paths)) {
  data.file <- data.processed.paths[[data.name]]
  sd.list[[data.name]] <- readRDS(file = data.file)
  DefaultAssay(sd.list[[data.name]]) <- def.assay
}



# ---- GBM-RCAS Module Scores ----


## ----  Pan AC markers ----
deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-NonNeoplastic-Level_3.csv")

# get top 10
topN <- 10
topGenes.pan <- deg_df %>%
  # dplyr::rename(Level_3 = cluster) %>%
  dplyr::filter(Level_3 == "Astrocyte") %>%
  .FactorizeMdata(.) %>%
  arrange(Level_3) %>%
  group_by(Level_3) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes.pan, n=50)



## ----  AC TEvNTvR markers ----
deg_df <- readr::read_csv("./results/07-1-deg-findmarkers/07-1-FindAllMarkers-AstrocyteSubtypes_TE_NT_R.csv")

# get top 10
topN <- 10
topGenes <- deg_df %>%
  dplyr::rename(Level_3AC = cluster) %>%
  .FactorizeMdata(.) %>%
  arrange(Level_3AC) %>%
  group_by(Level_3AC) %>%
  arrange(p_val_adj, .by_group = TRUE) %>%
  slice_head(n=topN)
print(topGenes, n=50)

topGenes$gene %in% topGenes.pan$gene


gene.list <- c(
  setNames(split(topGenes.pan, topGenes.pan$Level_3),
    paste0("pan_", names(split(topGenes.pan, topGenes.pan$Level_3)))),
  setNames(split(topGenes, topGenes$Level_3AC),
    paste0("AC_", names(split(topGenes, topGenes$Level_3AC))))
)
names(gene.list)
gene.list <- lapply(gene.list, function(df) {
  df$gene
})


cli_h1("Extracted Genes per Cluster")
for (cluster.i in names( gene.list)) {
  cli_alert_success("Cluster {.val {cluster.i}}: extracted {.strong {length(gene.list[[cluster.i]])}} genes ")
}



## --- Healthy ----
#module.scores <- lapply(deg.df.list, function(deg.df) {
object.h <- AddModuleScore(
  object = sd.list[["visium-hd-healthy"]],
  features = gene.list,
  # ctrl = deg.df$ctrl,
  assay = "Spatial.016um",
  name = glue("ModuleScore__")
)

colnames(object.tmp@meta.data)
module.score.df <- object.h@meta.data[, grep(glue("ModuleScore__"), colnames(object.h@meta.data) )]
module_names <- paste0(glue("ModuleScore__"), names(gene.list))
module_names <- gsub(" |-", "_", module_names)
colnames(module.score.df) <- module_names

# save module scores
# dir.create("./results/06-Visium")
saveRDS(module.score.df, file = "./manuscript-figures/figure-3/Figure-3-ModuleScores-VisiumHD-H-AstroTENTR.rds")

# Add metadata to object
object.h <- AddMetaData(object.h, module.score.df)
colnames(object.h@meta.data)


## --- Tumor ----
#module.scores <- lapply(deg.df.list, function(deg.df) {
object.t <- AddModuleScore(
  object = sd.list[["visium-hd-tumor"]],
  features = gene.list,
  # ctrl = deg.df$ctrl,
  assay = "Spatial.016um",
  name = glue("ModuleScore__")
)

colnames(object.t@meta.data)
module.score.df <- object.t@meta.data[, grep(glue("ModuleScore__"), colnames(object.t@meta.data) )]
module_names <- paste0(glue("ModuleScore__"), names(gene.list))
module_names <- gsub(" |-", "_", module_names)
colnames(module.score.df) <- module_names

# save module scores
# dir.create("./results/06-Visium")
saveRDS(module.score.df, file = "./manuscript-figures/figure-3/Figure-3-ModuleScores-VisiumHD-T-AstroTENTR.rds")

# Add metadata to object
object.t <- AddMetaData(object.t, module.score.df)
colnames(object.t@meta.data)


## --- Plot module scores spatial ----
#  plot module.scores
# NOTE: this 2-sample Visium HD AstroTENTR panel is SUPERSEDED by MF4 (the 4-sample rebuild,
#   script 12-3-Figure-3-visium-astrocyte-signature.R). The tracked PDF was removed from the repo;
#   re-running this block only produces an untracked local file, not the manuscript panel.
pdf.name <- file.path(plot.dir, glue("Figure-3-ModuleScores-VisiumHD-AstroTENTR.pdf"))
pdf(pdf.name, width = 14, height = 10)
for ( module_name in module_names) { # module_name <- module_names[1]
  cli::cli_alert("Plotting: {module_name}")

  h_score   <- FetchData(object.h,   vars = module_name, slot = "data")[,1]
  t_score   <- FetchData(object.h,   vars = module_name, slot = "data")[,1]
  # Use a robust cap so a few hot spots don’t blow out the legend
  # global_max <- as.numeric(quantile(c(h_score, t_score), probs = 0.999, na.rm = TRUE))
  # global_limits <- c(0, global_max)     # keep 0 on the low end color
  global_limits <- c(0, 0.75)     # keep 0 on the low end color



  p <- .spatialFeaturePlot(
    object.h,
    gene.i = module_name,
    limits = global_limits,
    pt.size.factor = 2.0,
    image.alpha = 0
  )
  print(p + ggtitle(glue("Healthy: {module_name}")))

  p <- .spatialFeaturePlot(
    object.t,
    gene.i = module_name,
    limits = global_limits,
    pt.size.factor = 2.0,
    image.alpha = 0
  )
  print(p + ggtitle(glue("Tumor: {module_name}")))

}
dev.off()



# ---- Pathway analysis ----

## ---- WaterfallPlot Astro TE v NT ----
matr <- readRDS(file.path("./results/07-2-gene-set-analysis/07-2-GeneSetAnalysis-AUCell-Hall50-Astrocytes.rds"))
str(matr)
seurat.sub <- subset(seurat.object, subset = Level_3AC %in%
    c("Astrocyte TE","Astrocyte NT", "Astrocyte R"))
seurat.sub <- .seuratFactorizeMdata(seurat.sub)

rownames(matr) <- gsub("HALLMARK_", "", rownames(matr))

p <- SeuratExtend::WaterfallPlot(
  matr,
  f = seurat.sub$Level_3AC,
  ident.1 = "Astrocyte TE",
  ident.2 = "Astrocyte NT",
  style = "segment",
  len.threshold = 10,
  top.n = 10
)
p
fig.name <- glue("Figure-3-Waterfall_Hall50_Astro_TEvNT")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 6, height = 4)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))




## ---- WaterfallPlot Astrocyte R vs Astrocyte TE+NT  ----
p <- SeuratExtend::WaterfallPlot(
  matr,
  f = seurat.sub$Level_3AC,
  ident.1 = "Astrocyte R",
  style = "segment",
  len.threshold = 10,
  top.n = 10
)
p
fig.name <- glue("Figure-3-Waterfall_Hall50_Astro_RvTENT")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 5, height = 3.8)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))




## ----  WaterfallPlot Astrocyte Reactive -----
## Tumor vs Healty
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

matr <- readRDS("./results/07-2-gene-set-analysis/07-2-GeneSetAnalysis-AUCell-Hall50-AstrocyteR.rds")
str(matr)
rownames(matr) <- gsub("HALLMARK_", "", rownames(matr))
p <- SeuratExtend::WaterfallPlot(
  matr, f = seurat.sub$disease_state,
  ident.1 = "Tumor",
  ident.2 = "Healthy",
  style = "segment",
  len.threshold = 5,
  top.n = 10
)
p
fig.name <- glue("Figure-3-Waterfall_Hall50_AstroR_TvH")
pdf.name <- file.path(plot.dir, glue("{fig.name}.pdf"))
ggsave(pdf.name, p, width = 5, height = 3.8)
readr::write_csv(p$data, sub("\\.pdf$", "-sourcedata.csv", pdf.name))



### ---- Volcanoes: Astrocyte Reactive DEG (Tumor-vs-Healthy + post-RT Recurrent-vs-Primary) ----
## Same seurat.sub + SeuratExtend::VolcanoPlot (companion to Fig 3H). Two things
## harmonised across the two panels so the side-by-side is a fair apples-to-apples:
##  (1) common p-cutoff. VolcanoPlot's default y.threshold is a per-panel 99th-percentile
##      of -log10(p), which for T-vs-H lands at p<9e-17 but for P->R at p<5e-4 — so each
##      panel would be scored on a different statistical bar. We instead FIX a common
##      y.threshold = -log10(0.01) (p<0.01) for BOTH; any gene passing in one panel would
##      pass in the other. (P->R genes are genuinely significant — Nr4a1 p=6e-13 etc.,
##      just fewer and smaller-fold than T-vs-H.)
##  (2) shared logFC (x) axis; two-tier x colouring + guide lines at BOTH cutoffs:
##      |logFC|>1 red (strong), 0.5<|logFC|<=1 paler red (modest), else grey.
## P->R y-axis is expanded to 25 (its data max ~12) to visually downplay it vs the much
## stronger T-vs-H panel. A source-data CSV (incl. the colour tier) is written per panel.
Y_CUT    <- -log10(0.01)          # common p-cutoff for both panels (p < 0.01)
X_MOD    <- 0.5; X_STRONG <- 1     # logFC guide lines / colour tiers
tier_cols <- c(ns = "grey72", modest = "#E8807B", strong = "#C0161B")  # grey / paler red / red

styled_volcano <- function(vp, ylim_max = NULL) {
  d <- vp$data
  d$tier <- dplyr::case_when(
    d$significant & abs(d$logFC) > X_STRONG ~ "strong",   # sig & |logFC|>1
    d$significant                           ~ "modest",   # sig & 0.5<|logFC|<=1
    TRUE                                    ~ "ns")
  d$tier <- factor(d$tier, levels = c("ns", "modest", "strong"))
  lab <- d[!is.na(d$label) & d$label != "", ]
  ggplot2::ggplot(d, ggplot2::aes(logFC, p)) +
    ggplot2::geom_point(ggplot2::aes(colour = tier), size = 0.7) +
    ggplot2::geom_vline(xintercept = c(-X_MOD, X_MOD),       linetype = "dotted", colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-X_STRONG, X_STRONG), linetype = "dashed", colour = "grey35", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = Y_CUT,                  linetype = "dashed", colour = "grey35", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = label, colour = tier),
                             size = 2.6, max.overlaps = Inf, min.segment.length = 0,
                             segment.size = 0.2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = tier_cols, guide = "none") +
    ggplot2::coord_cartesian(xlim = xr, ylim = if (!is.null(ylim_max)) c(0, ylim_max)) +
    ggplot2::labs(title = vp$labels$title, x = vp$labels$x, y = vp$labels$y) +
    ggplot2::theme_bw(base_size = 10)
}

Idents(seurat.sub) <- "disease_state"
vp_tvh  <- SeuratExtend::VolcanoPlot(seurat.sub, top.n = 20, ident.1 = "Tumor",     ident.2 = "Healthy",  x.threshold = X_MOD, y.threshold = Y_CUT)
Idents(seurat.sub) <- "sample_type"
vp_trtp <- SeuratExtend::VolcanoPlot(seurat.sub, top.n = 20, ident.1 = "Recurrent", ident.2 = "Primary",  x.threshold = X_MOD, y.threshold = Y_CUT)
Idents(seurat.sub) <- "disease_state"   # restore

xr <- range(c(vp_tvh$data$logFC, vp_trtp$data$logFC), na.rm = TRUE)
p_tvh  <- styled_volcano(vp_tvh)                 # T-vs-H keeps its own tall y-axis (~85)
p_trtp <- styled_volcano(vp_trtp, ylim_max = 25) # P->R expanded to 25 to downplay vs T-vs-H

ggsave(file.path(plot.dir, "Figure-3-FindMarkers-Volcano_AstroR_TvH.pdf"),    p_tvh,  width = 7, height = 6)
ggsave(file.path(plot.dir, "Figure-3-FindMarkers-Volcano_AstroR_TRvsTP.pdf"), p_trtp, width = 7, height = 6)
readr::write_csv(p_tvh$data,  file.path(plot.dir, "Figure-3-Volcano_AstroR_TvH-sourcedata.csv"))
readr::write_csv(p_trtp$data, file.path(plot.dir, "Figure-3-Volcano_AstroR_TRvsTP-sourcedata.csv"))




# ---- Supp figures ----

### ---- violin rcas call ----
sd.sub <- subset(seurat.object, subset = Level_3 == "Astrocyte" & Level_1 == "Non-Neoplastic")
gene.vec <- c("Gfap","Vim","Cd44","Sulf2",
  "C4b","Cebpd","S100a11","Ifitm3","S1pr3",
  "S100a6", "Nes", "Fabp7")
pdf.name <- file.path(plot.dir, glue("Figure-S2-ViolinPlots_Astrocyte_RCAS_Genes.pdf"))
pdf(pdf.name, width = 4, height = 5)
for (gene in gene.vec){
  p <- Seurat::VlnPlot(
    sd.sub,
    fill.by = "Level_3AC",
    flip=T, pt.size = 0,
    stack = F,
    cols = .color_pal[["Level_3AC"]],
    features =  gene,
    group.by = "Level_3AC"
  )
  plot(p)
}
dev.off()



