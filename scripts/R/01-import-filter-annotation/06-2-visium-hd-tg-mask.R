#!/usr/bin/env Rscript
# =============================================================================
# 06-2 — Visium HD transgene (Tg+) overlay mask -> QuPath GeoJSON (per sample)
# =============================================================================
# Generates the spatial TRANSGENE mask used to define tumour (Tg+) area in the
# Visium HD figures (Fig 1 transgene maps; Fig 7 spatial series).
#
# PROVENANCE (see STAR Methods "Spatial masks and histological annotation"):
#   1. Tissue region: HPC image-qc workflow (Otsu on the H&E), one candidate chosen
#      + polished in QuPath. Fold/blur EXCLUSIONS and histological regions annotated
#      MANUALLY in QuPath. (Not code; described in STAR.)
#   2. Tg+ overlay (THIS script): sum the RCAS transgene probes (Tg-RFP, Tg-hPDGFB)
#      per 16um bin from SpaceRanger binned counts, call a bin Tg+ at >= TG_THRESHOLD
#      summed counts, union the Tg+ bin squares into one MultiPolygon, write a
#      QuPath-loadable GeoJSON in full-res pixel coordinates.
#   3. In QuPath the Tg+ overlay is imported and (24um) expanded to the
#      `Tg_positive_expanded` annotation used downstream for tumour membership.
#
# RERUN: regeneration needs SpaceRanger `binned_outputs/square_016um/` (h5 +
# tissue_positions + scalefactors) AND the geo R packages (sf/arrow/Matrix); square-
# bin counts are run-invariant, so the DEFAULT run is canonical. When SpaceRanger
# output is not present locally, this script instead STAGES the tracked canonical
# overlay from references/visium-hd-masks/tg-overlays into results/ so the
# figures are reproducible (the staging path needs NO extra packages). Set
# SPACERANGER_ROOT to force regeneration.
# Output (gitignored): results/06-2-tg-masks/{sample}_tg_overlay.geojson
# =============================================================================
say <- function(...) cat(sprintf(...), "\n")
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
out_dir   <- file.path(proj_root, "results", "06-2-tg-masks"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- params (mirror spatial-masks run_01_default.conf) ----------------------
SAMPLES      <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")   # publication cohort
BIN_SIZE_UM  <- 16L
TG_THRESHOLD <- 2L
TG_PATTERNS  <- "Tg-RFP|Tg-hPDGFB|CUSTOMPROBE-RFP"
SPACERANGER_ROOT <- Sys.getenv("SPACERANGER_ROOT", unset = "")   # {ROOT}/{sample}/outs/binned_outputs/square_016um/...
CACHE_TG <- file.path(proj_root, "references/visium-hd-masks/tg-overlays")  # tracked canonical overlays

# ---- Tg+ overlay generation (ported algo; heavy geo deps loaded lazily) ------
generate_overlay <- function(sid, sr_outs) {
  pkgs <- c("arrow", "sf", "jsonlite", "Matrix", "dplyr", "Seurat")
  if (!all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) {
    say("  [%s] geo packages (sf/arrow/...) not installed -> cannot regenerate", sid); return(NULL) }
  sf::sf_use_s2(FALSE)
  bin_dir  <- file.path(sr_outs, "binned_outputs", sprintf("square_%03dum", BIN_SIZE_UM))
  pos_path <- file.path(bin_dir, "spatial", "tissue_positions.parquet")
  h5_path  <- file.path(bin_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(pos_path) || !file.exists(h5_path)) return(NULL)
  mat <- Seurat::Read10X_h5(h5_path, use.names = TRUE); if (is.list(mat)) mat <- mat[["Gene Expression"]]
  tg_idx <- grep(TG_PATTERNS, rownames(mat))
  tg_sum <- if (length(tg_idx)) as.integer(Matrix::colSums(mat[tg_idx, , drop = FALSE])) else 0L
  tg  <- data.frame(barcode = colnames(mat), tg_sum = tg_sum)
  pos <- arrow::read_parquet(pos_path); if ("in_tissue" %in% names(pos)) pos <- dplyr::filter(pos, in_tissue == 1)
  merged <- dplyr::filter(dplyr::inner_join(pos, tg, by = "barcode"), tg_sum >= TG_THRESHOLD)
  say("  [%s] Tg+ bins (sum>=%d): %d / %d", sid, TG_THRESHOLD, nrow(merged), nrow(pos))
  if (nrow(merged) == 0) return("no-tg")   # healthy: legitimately no Tg+ region
  sfj <- jsonlite::fromJSON(file.path(bin_dir, "spatial", "scalefactors_json.json"))
  bin_px <- if (!is.null(sfj$bin_size_um) && !is.null(sfj$microns_per_pixel)) as.numeric(sfj$bin_size_um)/as.numeric(sfj$microns_per_pixel) else as.numeric(sfj$spot_diameter_fullres)
  half <- bin_px / 2
  polys <- lapply(seq_len(nrow(merged)), function(j) { x <- merged$pxl_col_in_fullres[j]; y <- merged$pxl_row_in_fullres[j]
    sf::st_polygon(list(matrix(c(x-half,y-half, x+half,y-half, x+half,y+half, x-half,y+half, x-half,y-half), ncol = 2, byrow = TRUE))) })
  geom_sf <- sf::st_sf(geometry = sf::st_union(sf::st_sfc(polys)))
  out_path <- file.path(out_dir, paste0(sid, "_tg_overlay.geojson"))
  tmp <- tempfile(fileext = ".geojson")
  sf::st_write(geom_sf, tmp, driver = "GeoJSON", quiet = TRUE, append = FALSE, layer_options = "COORDINATE_PRECISION=1")
  fc <- jsonlite::fromJSON(tmp, simplifyVector = FALSE); feat <- fc$features[[1]]
  hx <- function(n) paste(sample(c(0:9, letters[1:6]), n, replace = TRUE), collapse = "")
  feat[["id"]] <- sprintf("%s-%s-4%s-%s%s-%s", hx(8), hx(4), hx(3), sample(c("8","9","a","b"),1), hx(3), hx(12))
  feat[["properties"]] <- list(name = "Tg_positive", objectType = "annotation")
  jsonlite::write_json(list(feat), out_path, auto_unbox = TRUE, digits = 1, pretty = FALSE); unlink(tmp)
  out_path
}

# ---- per sample: regenerate if SpaceRanger present, else stage from cache ----
say("== 06-2 Tg overlay masks -> %s", out_dir)
for (sid in SAMPLES) {
  sr_outs <- if (nzchar(SPACERANGER_ROOT)) file.path(SPACERANGER_ROOT, sid, "outs") else ""
  res <- if (nzchar(sr_outs) && dir.exists(sr_outs)) generate_overlay(sid, sr_outs) else NULL
  if (!is.null(res) && !identical(res, "no-tg")) { say("  [%s] regenerated from SpaceRanger", sid); next }
  if (identical(res, "no-tg")) { say("  [%s] no Tg+ bins (healthy) - no overlay", sid); next }
  cached <- file.path(CACHE_TG, paste0(sid, "_tg_overlay.geojson"))
  if (file.exists(cached)) { file.copy(cached, file.path(out_dir, basename(cached)), overwrite = TRUE)
    say("  [%s] SpaceRanger input absent -> STAGED tracked canonical overlay (set SPACERANGER_ROOT to regenerate)", sid)
  } else say("  [%s] no SpaceRanger input AND no tracked overlay (%s) - healthy sample has none", sid, cached)
}
# The MANUAL QuPath annotation exports (tissue / fold / histology; not code, see STAR)
# are tracked directly in references/visium-hd-masks/qupath-annotations/ (per-sample
# GeoJSON) and read from there by the figure layer - no staging step needed.
say("== 06-2 done. Tg masks in results/06-2-tg-masks/ (gitignored)")
