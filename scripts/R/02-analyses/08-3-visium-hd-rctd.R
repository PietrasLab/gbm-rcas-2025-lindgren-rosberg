# 08-3 — Visium HD RCTD (spacexr) bin deconvolution, 14-type reference
# ============================================================================
# Per-bin cell-type deconvolution of the 4 publication Visium HD samples
# (HB_hm01 healthy, TP_987 + TP_1083 primary, TR_03d_990 post-RT recurrent),
# full mode, 16um bins, in-tissue bins only. Produces the per-bin weight tables
# consumed by the Figure 7 scripts (16-2 / 16-3 / 16-6) and by 08-4 / 08-5.
#
# REFERENCE ANNOTATION — 14-type curated collapse of Level_4ACM:
#   Neoplastic (5): Neopl-ACR, Neopl-ECM, Neopl-OPC, Neopl-COP, Neopl-Bulk
#     - Neopl-OPC and Neopl-COP kept SEPARATE (they split on the reference UMAP,
#       and this is a proneural/OPC-biased PDGFB model).
#     - Neopl-Bulk = Bulk + NC + RNA-low + CC-I/II/III. NC folded because it is
#       transcriptionally degenerate with normal Neural (spurious healthy-neopl
#       calls) and Bulk already shares NC markers (Fig 4B). CC folded because the
#       cell-cycle program is cross-cutting, not a readable spatial state.
#   Non-neoplastic (9): Astrocyte R, Astrocyte TE-NT, OPC-COP-OLG, Neuron,
#     Ependymal, Choroid, Myeloid, Immune-other, Vascular.
#
# INPUTS:
#   - Flex scRNA reference (Zenodo): data/processed/seurat/seurat_flex_filtered_v1.0.rds
#   - Per-sample 16um Seurat objects (Zenodo): data/processed/seurat/visium-hd
#     (env VISIUM_HD_OBJ_ROOT); regenerate from SpaceRanger via SPACERANGER_ROOT.
#   - Tissue masks (tracked): references/visium-hd-masks/masks-workspace
#     (env VISIUM_HD_MASK_DIR)
# OUTPUT (results/08-visium-hd-derived/, gitignored, regenerable):
#   rctd-weights-<sample>-16um.csv.gz + weight-summary / validation CSVs
# NB: RCTD is compute-heavy (~1e5 bins x 4 samples; hours, multicore).
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(spacexr); library(Matrix)
  library(dplyr); library(readr); library(tidyr); library(cli)
})

proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
source(file.path(proj_root, "scripts/R/00-2-mask-utils.R"))
out_dir <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- config ----------------------------------------------------------------
OBJ_ROOT  <- Sys.getenv("VISIUM_HD_OBJ_ROOT",
  unset = file.path(proj_root, "data/processed/seurat/visium-hd"))
SPACERANGER_ROOT <- Sys.getenv("SPACERANGER_ROOT", unset = "")  # optional: {ROOT}/{sample}/outs
mask_dir  <- Sys.getenv("VISIUM_HD_MASK_DIR",
  unset = file.path(proj_root, "references/visium-hd-masks/masks-workspace"))
ref_path  <- file.path(proj_root, "data/processed/seurat/seurat_flex_filtered_v1.0.rds")
BIN_SIZE  <- 16L
REF_LABEL_COL <- "Level_4ACM"
RCTD_MODE     <- "full"
RCTD_MAX_CORES <- max(1L, parallel::detectCores() - 1L)

is_neoplastic <- function(types) grepl("^Neopl", types)
REACTIVE_ASTRO_TYPE <- "Astrocyte R"
sanitize_ct <- function(x) gsub("/", "-", x)   # RCTD forbids '/'

SAMPLES <- c("TP_987", "TR_03d_990", "HB_hm01", "TP_1083")
.cond <- function(s) ifelse(grepl("^HB", s), "Healthy",
                     ifelse(grepl("^TP", s), "Primary", "PostRT"))

set.seed(1)
cli::cli_h1("RCTD 14-type deconvolution — {length(SAMPLES)} sample(s): {paste(SAMPLES, collapse=', ')}")

# ============================================================================
# 1. Build the RCTD Reference (Flex scRNA raw counts + Level_4ACM, 14-type collapse)
# ============================================================================
cli::cli_h2("Building RCTD reference from {basename(ref_path)} [{REF_LABEL_COL}]")
ref_seu <- readRDS(ref_path)
ref_counts <- LayerData(ref_seu, layer = "counts", assay = "RNA")
ref_labels <- ref_seu@meta.data[[REF_LABEL_COL]]
stopifnot(!is.null(ref_labels), length(ref_labels) == ncol(ref_counts))
keep_cells <- !is.na(ref_labels) & !(ref_labels %in% c("Ambiguous"))
ref_counts <- ref_counts[, keep_cells]
sane <- sanitize_ct(ref_labels[keep_cells])

# ---- 14-type collapse (keys are SANITIZED labels; '/' already -> '-') -------
collapse_map <- c(
  # neoplastic (5): OPC and COP kept SEPARATE; NC + RNA-low + CC -> Neopl-Bulk
  "Neopl-ACR" = "Neopl-ACR", "Neopl-ECM" = "Neopl-ECM",
  "Neopl-OPC" = "Neopl-OPC", "Neopl-COP" = "Neopl-COP",
  "Neopl-Bulk" = "Neopl-Bulk", "Neopl-NC" = "Neopl-Bulk", "Neopl-RNA-low" = "Neopl-Bulk",
  "Neopl-CC-I" = "Neopl-Bulk", "Neopl-CC-II" = "Neopl-Bulk", "Neopl-CC-III" = "Neopl-Bulk",
  # non-neoplastic (9)
  "Astrocyte R" = "Astrocyte R", "Astrocyte TE-NT" = "Astrocyte TE-NT",
  "OPC-COP-OLG" = "OPC-COP-OLG", "Neural" = "Neuron",
  "Ependymal" = "Ependymal", "Choroid" = "Choroid",
  "Macrophage" = "Myeloid", "Microglia" = "Myeloid",
  "Dendritic" = "Immune-other", "Neutrophil" = "Immune-other", "NKTB" = "Immune-other",
  "Endothelial" = "Vascular", "Mural" = "Vascular", "Fibroblast" = "Vascular"
)
unmapped <- setdiff(unique(sane), names(collapse_map))
if (length(unmapped)) cli::cli_alert_warning("unmapped ref labels kept as-is: {paste(unmapped, collapse=', ')}")
curated <- ifelse(sane %in% names(collapse_map), unname(collapse_map[sane]), sane)
cell_types <- factor(curated)
names(cell_types) <- colnames(ref_counts)
ref_nUMI <- colSums(ref_counts)

cli::cli_alert_info("Reference: {ncol(ref_counts)} cells x {nrow(ref_counts)} genes, {nlevels(cell_types)} cell types")
print(sort(table(cell_types), decreasing = TRUE))

reference <- spacexr::Reference(counts = ref_counts, cell_types = cell_types,
                                nUMI = ref_nUMI, require_int = TRUE)
rm(ref_seu, ref_counts); gc(FALSE)

neopl_levels <- levels(cell_types)[is_neoplastic(levels(cell_types))]
cli::cli_alert_info("Neoplastic types ({length(neopl_levels)}): {paste(neopl_levels, collapse=', ')}")

# ============================================================================
# 2. Per-sample: load spatial bins, tissue-mask, run RCTD full-mode
# ============================================================================
run_one <- function(sid, cond) {
  cli::cli_h2("{sid} [{cond}]")
  obj_path <- file.path(OBJ_ROOT, sprintf("seurat_visium-hd_%s_%dum_v1.1.rds", sid, BIN_SIZE))
  obj <- if (file.exists(obj_path)) {
    cli::cli_alert_info("reading {basename(obj_path)}"); readRDS(obj_path)
  } else {
    stopifnot(nzchar(SPACERANGER_ROOT))
    cli::cli_alert_info("object absent - Load10X_Spatial from SPACERANGER_ROOT")
    Load10X_Spatial(file.path(SPACERANGER_ROOT, sid, "outs"), bin.size = BIN_SIZE)
  }
  assay <- grep("Spatial", Assays(obj), value = TRUE)[1]
  obj <- apply_tissue_mask(obj, sid, mask_dir, bin_size = BIN_SIZE)
  obj <- attach_region_labels(obj, sid, mask_dir, bin_size = BIN_SIZE)

  sp_counts <- LayerData(obj, layer = "counts", assay = assay)
  coords    <- GetTissueCoordinates(obj)[, c("x", "y")]
  coords    <- coords[colnames(sp_counts), , drop = FALSE]
  sp_nUMI   <- colSums(sp_counts)
  ok <- sp_nUMI > 0
  sp_counts <- sp_counts[, ok]; coords <- coords[ok, ]; sp_nUMI <- sp_nUMI[ok]

  query <- spacexr::SpatialRNA(coords = coords, counts = sp_counts, nUMI = sp_nUMI)
  cli::cli_alert_info("RCTD create + run [{RCTD_MODE}] on {ncol(sp_counts)} in-tissue bins")
  rctd <- spacexr::create.RCTD(query, reference, max_cores = RCTD_MAX_CORES, test_mode = FALSE)
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = RCTD_MODE)

  W  <- as.matrix(rctd@results$weights)
  Wn <- as.matrix(spacexr::normalize_weights(W))
  bc <- rownames(Wn)

  neopl_frac <- if (length(neopl_levels)) rowSums(Wn[, intersect(neopl_levels, colnames(Wn)), drop = FALSE]) else rep(0, nrow(Wn))
  ra_frac    <- if (REACTIVE_ASTRO_TYPE %in% colnames(Wn)) Wn[, REACTIVE_ASTRO_TYPE] else rep(NA_real_, nrow(Wn))

  md <- obj@meta.data[bc, , drop = FALSE]
  in_tumor <- if ("in_tumor" %in% colnames(md)) md$in_tumor else NA
  weights_df <- cbind(
    tibble(sample_id = sid, condition = cond, barcode = bc,
           neopl_frac = neopl_frac, reactive_astro_frac = ra_frac,
           in_tumor = in_tumor, x = coords[bc, "x"], y = coords[bc, "y"]),
    as.data.frame(Wn))
  wpath <- file.path(out_dir, sprintf("rctd-weights-%s-%dum.csv.gz", sid, BIN_SIZE))
  write_csv(weights_df, wpath)
  cli::cli_alert_success("weights -> {basename(wpath)}")
  cli::cli_alert("median neopl_frac = {round(median(neopl_frac),3)} | median reactive_astro_frac = {round(median(ra_frac, na.rm=TRUE),3)}")
  list(sid = sid, cond = cond, Wn = Wn, neopl_frac = neopl_frac, ra_frac = ra_frac)
}

results <- list()
for (sid in SAMPLES) { results[[sid]] <- run_one(sid, .cond(sid)); gc(FALSE) }

# ============================================================================
# 3. QC / sanity outputs
# ============================================================================
cli::cli_h2("QC / sanity outputs")
summ <- lapply(results, function(r) {
  tibble(sample_id = r$sid, condition = r$cond, cell_type = colnames(r$Wn),
         mean_weight = colMeans(r$Wn), median_weight = apply(r$Wn, 2, median),
         pct_bins_gt10 = colMeans(r$Wn > 0.10) * 100, max_weight = apply(r$Wn, 2, max))
}) |> bind_rows()
write_csv(summ, file.path(out_dir, "rctd-weight-summary-by-celltype.csv"))

# ============================================================================
# 4. VALIDATION — neoplastic weight healthy-vs-tumor (+ OPC/COP degeneracy check)
# ============================================================================
cli::cli_h2("VALIDATION: neoplastic weight healthy-vs-tumor")
val <- lapply(results, function(r) {
  tibble(sample_id = r$sid, condition = r$cond, n_bins = length(r$neopl_frac),
         mean_neopl_frac = mean(r$neopl_frac), median_neopl_frac = median(r$neopl_frac),
         pct_bins_neopl_gt50 = mean(r$neopl_frac > 0.5) * 100,
         mean_reactive_astro = mean(r$ra_frac, na.rm = TRUE))
}) |> bind_rows()
write_csv(val, file.path(out_dir, "rctd-validation-neoplastic-by-sample.csv"))
print(as.data.frame(val))
for (i in seq_len(nrow(val))) {
  s <- val[i, ]
  if (s$condition == "Healthy" && s$mean_neopl_frac > 0.10)
    cli::cli_alert_danger("FLAG {s$sample_id} (healthy): mean neopl_frac = {round(s$mean_neopl_frac,3)} (expected ~0)")
  if (s$condition == "Healthy" && s$mean_neopl_frac <= 0.10)
    cli::cli_alert_success("OK {s$sample_id} (healthy): mean neopl_frac = {round(s$mean_neopl_frac,3)} (~0)")
  if (s$condition %in% c("Primary","PostRT") && s$mean_neopl_frac >= 0.15)
    cli::cli_alert_success("OK {s$sample_id} (tumor): mean neopl_frac = {round(s$mean_neopl_frac,3)} (substantial)")
}

# OPC vs COP degeneracy diagnostic (per tumor sample): correlation of the two
# neoplastic weights across bins. Strong NEGATIVE r = the pair may be degenerate.
if (all(c("Neopl-OPC","Neopl-COP") %in% colnames(results[[1]]$Wn))) {
  cli::cli_h2("Neopl-OPC vs Neopl-COP per-bin correlation (degeneracy check)")
  oc <- lapply(results, function(r) tibble(sample_id = r$sid, condition = r$cond,
        r_opc_cop = suppressWarnings(cor(r$Wn[,"Neopl-OPC"], r$Wn[,"Neopl-COP"], method = "spearman")))) |> bind_rows()
  write_csv(oc, file.path(out_dir, "rctd-opc-cop-degeneracy.csv"))
  print(as.data.frame(oc))
}

cli::cli_alert_success("RCTD 14-type deconvolution complete -> {out_dir}")
