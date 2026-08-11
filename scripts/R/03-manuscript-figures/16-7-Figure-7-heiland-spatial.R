# 16-7 — Figure 7F (MF5): Heiland human GBM spatial module-score maps (UKF260_T_ST)
# ============================================================================
# The published Fig 7F / supplement S6D, reworked with the harmonised + REFINED signatures. Representative
# Heiland tumour section UKF260_T_ST, rotated to the published orientation, the four modules one per page:
# Richards Developmental, Richards Injury-Response, Neftel MES1, and Neopl-ACR (refined). gray15 near-black
# background, default Seurat Spectral palette (matching the original 08-2 styling).
#
# Reads the refined-signature panel builder (_signature-panel.R) and the pre-normalised Heiland
# .qs Seurat object from data/external/10x-visium-heiland (env-overridable). The image
# is hidden, so the spots are replotted directly from Seurat coords to apply the rotation (px = -x, py = -y).
#
# NOTE: this repo's Seurat cannot run Load10X_Spatial on this data (FOV object-model
# mismatch), so this script reads a pre-normalised .qs object exported upstream; with
# that object it runs and reproduces the committed panel directly in this environment.
# ============================================================================
suppressPackageStartupMessages({ library(Seurat); library(qs); library(dplyr); library(stringr); library(ggplot2); library(RColorBrewer); library(readr); library(cli) })
set.seed(169)
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
source(file.path(proj_root, "scripts/R/00-0-source-functions.R"))
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
HEILAND   <- Sys.getenv("HEILAND_DIR", unset = file.path(proj_root, "data/external/10x-visium-heiland"))
source(file.path(proj_root, "scripts/R/03-manuscript-figures/_signature-panel.R"))   # build_panel_full (harmonised + refined signatures)

SPATIAL_COLS <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)
panel   <- build_panel_full("human", acr_targets = "both")
# the 4 published 7F modules
MODULES <- intersect(c("Richards Developmental","Richards Injury-Response","Neftel MES1","Neopl-ACR (ours)"), names(panel$sigs))
sigs    <- panel$sigs[MODULES]

SID <- "UKF260_T_ST"   # representative Heiland tumour section (the published 7F / S6D example)

one_map <- function(o, feat) {
  co <- GetTissueCoordinates(o); cells <- if ("cell" %in% names(co)) co$cell else colnames(o)
  df <- data.frame(px = -co$x, py = -co$y, score = o@meta.data[cells, feat])
  ggplot(df, aes(px, py, colour = score)) + geom_point(size = 1.4) +
    scale_colour_gradientn(colours = SPATIAL_COLS, name = NULL) + coord_fixed() + theme_void() +
    theme(panel.background=element_rect(fill="gray15",color=NA), plot.background=element_rect(fill="gray15",color=NA),
          legend.background=element_rect(fill="gray15"), legend.position="bottom",
          legend.title=element_text(size=10,color="white"), legend.text=element_text(size=8,color="white"),
          plot.title=element_text(size=12,face="bold",hjust=0.5,color="white")) +
    ggtitle(sub("\\(ours\\)","(refined)",feat))
}

cli::cli_alert("Fig 7F: {SID}")
# This repo's Seurat cannot Load10X_Spatial this data (FOV object-model mismatch), so read the
# pre-normalised .qs Seurat object exported upstream. GetTissueCoordinates() -> x,y.
o <- qs::qread(file.path(HEILAND, sprintf("00-seurat-Heiland_%s-norm-object.qs", SID)))
assay <- grep("Spatial", Assays(o), value=TRUE)[1]; DefaultAssay(o) <- assay
o <- AddModuleScore(o, features=lapply(sigs, function(g) intersect(g, rownames(o))), name="M_", assay=assay)
for (k in seq_along(MODULES)) o[[MODULES[k]]] <- o[[paste0("M_",k)]]
present <- MODULES[vapply(sigs, function(g) length(intersect(g, rownames(o)))>=5, logical(1))]

# four modules, one per page, at full per-panel resolution
pdf(file.path(plot_dir, "Figure-7-heiland-spatial-UKF260-4sigs.pdf"), width=6, height=6)
for (m in intersect(MODULES, present)) print(one_map(o, m))
invisible(dev.off())

# source data: per-spot rotated coords + module score, one row per spot x module
co <- GetTissueCoordinates(o); cells <- if ("cell" %in% names(co)) co$cell else colnames(o)
sd <- bind_rows(lapply(intersect(MODULES, present), function(m)
  data.frame(module = sub("\\(ours\\)","(refined)",m), px = -co$x, py = -co$y, score = o@meta.data[cells, m])))
write_csv(sd, file.path(plot_dir, "Figure-7-heiland-spatial-UKF260-4sigs-sourcedata.csv"))
cli::cli_alert_success("16-7 Fig 7F Heiland spatial (MF5) -> {plot_dir}")
