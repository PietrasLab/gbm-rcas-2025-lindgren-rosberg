# 08-6 — Heiland human GBM Visium: per-section signature-panel correlations
# ============================================================================
# For each TUMOR-group Heiland/Ravi section (Tumor = _T_/_TC_/_TI_; excludes
# Control _C_ and IDHMutant), scores the harmonized HUMAN signature panel
# (build_panel_full, refined + original ACR/Astrocyte-R targets) with
# AddModuleScore and saves the per-section Spearman correlation matrices.
# Consumed by the Figure 7E ranked-correlation script (16-5).
#
# INPUTS: raw Heiland sections in data/external/10x-visium-heiland
#   (env HEILAND_DIR; per-section <id>/outs/filtered_feature_bc_matrix.h5).
#   Loading falls back to plain Read10X_h5 (no image/coords needed here).
# OUTPUT (results/08-visium-hd-derived/): heiland-per-section-cor-list.rds
# ============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(stringr); library(cli) })
set.seed(169)
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
source(file.path(proj_root, "scripts/R/03-manuscript-figures/_signature-panel.R"))
out_dir   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
HEILAND   <- Sys.getenv("HEILAND_DIR", unset = file.path(proj_root, "data/external/10x-visium-heiland"))

panel <- build_panel_full("human", acr_targets = "both")

group_of_id <- function(sid) dplyr::case_when(
  stringr::str_detect(sid, "_IDHMutant_") ~ "IDHMutant",
  stringr::str_detect(sid, "_T[CI]?_")     ~ "Tumor",
  stringr::str_detect(sid, "_C_")          ~ "Control", TRUE ~ "Other")
dirs <- list.dirs(HEILAND, recursive = FALSE, full.names = TRUE)
sids <- sub("^#", "", basename(dirs)); keep <- group_of_id(sids) == "Tumor"
dirs <- dirs[keep]; sids <- sids[keep]
cli::cli_h1("Heiland TUMOR-group sections: {length(dirs)}  (panel {length(panel$sigs)} sigs, refined + original targets)")

load_heiland <- function(dpath) {   # Load10X_Spatial; fallback to plain h5 (no coords needed)
  outs <- file.path(dpath, "outs"); h5 <- list.files(outs, pattern = "filtered_feature_bc_matrix.h5$", full.names = FALSE)
  tryCatch(Load10X_Spatial(outs, filename = h5), error = function(e)
    tryCatch(CreateSeuratObject(Read10X_h5(file.path(outs, h5)), assay = "Spatial"), error = function(e2) NULL))
}
cor_list <- list()
for (i in seq_along(dirs)) {
  sid <- sids[i]; cli::cli_alert("{sid} ({i}/{length(dirs)})")
  o <- load_heiland(dirs[i]); if (is.null(o) || ncol(o) < 100) { cli::cli_alert_warning("{sid}: load failed/small"); next }
  assay <- grep("Spatial", Assays(o), value = TRUE)[1]; DefaultAssay(o) <- assay
  o <- NormalizeData(o, verbose = FALSE)
  o <- AddModuleScore(o, features = lapply(panel$sigs, function(g) intersect(g, rownames(o))), name = "S_", assay = assay)
  msdf <- o@meta.data[, paste0("S_", seq_along(panel$sigs))]; colnames(msdf) <- names(panel$sigs)
  cor_list[[sid]] <- cor(msdf, method = "spearman"); rm(o); gc(FALSE)
}
saveRDS(cor_list, file.path(out_dir, "heiland-per-section-cor-list.rds"))
cli::cli_alert_success("08-6: scored {length(cor_list)} Heiland tumor sections -> heiland-per-section-cor-list.rds")
