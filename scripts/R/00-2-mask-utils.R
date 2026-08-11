# 00-2-mask-utils.R
# Utility: apply tissue masks to Visium HD Seurat objects
#
# Parquet files (produced externally, one per sample) contain per-bin boolean
# mask columns. mask_selections.csv specifies which mask column to use per sample.
#
# Usage:
#   source("scripts/R/00-2-mask-utils.R")
#   mask_bcs <- get_mask_barcodes(sid, mask_dir)   # → character vector
#   obj_masked <- subset(obj, cells = mask_bcs)
#
# mask_name → parquet column mapping (letter schedule):
#   A  in_tissue              → col "in_tissue"             SpaceRanger default
#   B  sat_in_tissue          → col "sat_in_tissue"         Saturation Otsu
#   C  otsu_in_tissue         → col "otsu_in_tissue"        Grayscale Otsu
#   D  hue_in_tissue          → col "hue_in_tissue"         Hue-gated
#   E  imgqc_in_tissue        → col "imgqc_in_tissue"       sat|otsu|hue union
#   F  umi_in_tissue_refined  → col "umi_in_tissue_refined" UMI Otsu refined
#   G  tissue_joint           → col "tissue_joint"          SpaceRanger & UMI refined
#   H  image_union_and_umi    → COMPUTED: imgqc & umi_in_tissue_refined
#   M  manual                 → NPZ-based (not yet implemented)
#   S  skip                   → no masking applied
#
# Note: parquet contains ALL bins; Seurat objects are pre-filtered to
# spaceranger in_tissue bins. Masking further subsets to the selected mask.

# Letter → parquet column (or special handler tag)
# mask_letter is the authoritative field; mask_name is descriptive only.
.LETTER_TO_COL <- c(
  A = "in_tissue",
  B = "sat_in_tissue",
  C = "otsu_in_tissue",
  D = "hue_in_tissue",
  E = "imgqc_in_tissue",
  F = "umi_in_tissue_refined",
  G = "tissue_joint",
  H = "__image_union_and_umi__",  # computed: imgqc & umi_in_tissue_refined
  M = "__manual__",               # NPZ-based (not yet implemented)
  Q = "qupath_in_tissue",         # QuPath GeoJSON annotation (rasterized upstream)
  S = "__skip__"                  # no masking
)


# ---- get_mask_barcodes() -----------------------------------------------
# Returns a character vector of barcodes that pass the selected mask for `sid`.
# Returns NULL (no masking) if: parquet missing, letter S, or letter M.
#
# Args:
#   sid       — sample ID matching mask_selections.csv and parquet filename
#   mask_dir  — directory containing mask_selections.csv and *_umi_per_bin_016um.parquet
#   bin_size  — bin size in µm (default 16)
#
get_mask_barcodes <- function(sid, mask_dir, bin_size = 16) {
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' is required. Install with: install.packages('arrow')")

  sel_path <- file.path(mask_dir, "mask_selections.csv")
  if (!file.exists(sel_path))
    stop("mask_selections.csv not found in: ", mask_dir)

  sel <- read.csv(sel_path, stringsAsFactors = FALSE)
  row <- sel[sel$sample_id == sid, ]
  if (nrow(row) == 0) {
    warning("No mask entry for sample: ", sid, " — returning NULL")
    return(NULL)
  }

  letter    <- toupper(trimws(row$mask_letter[1]))
  mask_name <- row$mask_name[1]   # human-readable label only

  col_tag <- .LETTER_TO_COL[letter]
  if (is.na(col_tag)) {
    warning("Unknown mask_letter '", letter, "' for ", sid,
            " — returning NULL.")
    return(NULL)
  }

  # S — skip
  if (col_tag == "__skip__") {
    cat(sprintf("  [mask] %s: S (skip) — no mask applied\n", sid))
    return(NULL)
  }

  # M — manual (NPZ not yet implemented)
  if (col_tag == "__manual__") {
    warning("Manual mask (M) for ", sid,
            " requires NPZ → coordinate mapping (not yet implemented). ",
            "Returning NULL — object will not be masked.")
    return(NULL)
  }

  # Q — QuPath GeoJSON-derived mask (column rasterized upstream into the parquet)
  if (col_tag == "qupath_in_tissue") {
    pq_path_q <- file.path(mask_dir,
                   paste0(sid, "_umi_per_bin_", sprintf("%03d", bin_size),
                          "um.parquet"))
    if (!file.exists(pq_path_q)) {
      warning("Parquet not found for Q mask, ", sid, ": ", pq_path_q,
              "\nReturning NULL — object will not be masked.")
      return(NULL)
    }
    pq <- arrow::read_parquet(pq_path_q, col_select = c("barcode", "qupath_in_tissue"))
    if (!"qupath_in_tissue" %in% names(pq)) {
      warning("Column 'qupath_in_tissue' missing from parquet for ", sid,
              ". Returning NULL.")
      return(NULL)
    }
    keep <- as.logical(pq$qupath_in_tissue)
    keep[is.na(keep)] <- FALSE
    cat(sprintf("  [mask] %s: Q/qupath_in_tissue → %d / %d bins kept\n",
                sid, sum(keep), nrow(pq)))
    return(pq$barcode[keep])
  }

  pq_path <- file.path(mask_dir,
               paste0(sid, "_umi_per_bin_", sprintf("%03d", bin_size),
                      "um.parquet"))
  if (!file.exists(pq_path)) {
    warning("Parquet not found for ", sid, ": ", pq_path,
            "\nReturning NULL — object will not be masked.")
    return(NULL)
  }

  # H — image_union_and_umi: imgqc_in_tissue AND umi_in_tissue_refined
  if (col_tag == "__image_union_and_umi__") {
    pq   <- arrow::read_parquet(pq_path,
                                 col_select = c("barcode", "imgqc_in_tissue",
                                                "umi_in_tissue_refined"))
    keep <- as.logical(pq$imgqc_in_tissue) & as.logical(pq$umi_in_tissue_refined)
    keep[is.na(keep)] <- FALSE
    cat(sprintf("  [mask] %s: H/image_union_and_umi (imgqc & umi_refined) → %d / %d bins\n",
                sid, sum(keep), nrow(pq)))
    return(pq$barcode[keep])
  }

  # A–G — direct parquet column
  pq   <- arrow::read_parquet(pq_path, col_select = c("barcode", unname(col_tag)))
  keep <- as.logical(pq[[col_tag]])
  keep[is.na(keep)] <- FALSE
  cat(sprintf("  [mask] %s: %s/%s ('%s') → %d / %d bins kept\n",
              sid, letter, mask_name, col_tag, sum(keep), nrow(pq)))
  pq$barcode[keep]
}


# ---- apply_tissue_mask() -----------------------------------------------
# Convenience wrapper: loads barcodes and subsets a Seurat object.
# If masking is unavailable (NULL barcodes), returns the object unchanged.
#
# Args:
#   obj      — Seurat object
#   sid      — sample ID
#   mask_dir — mask directory
#   bin_size — bin size in µm (default 16)
#
apply_tissue_mask <- function(obj, sid, mask_dir, bin_size = 16) {
  bcs <- get_mask_barcodes(sid, mask_dir, bin_size)
  if (is.null(bcs)) {
    cat(sprintf("  [mask] %s: no mask applied\n", sid))
    return(obj)
  }
  keep_bcs <- intersect(bcs, colnames(obj))
  dropped  <- ncol(obj) - length(keep_bcs)
  cat(sprintf("  [mask] %s: dropping %d non-tissue bins (%d → %d)\n",
              sid, dropped, ncol(obj), length(keep_bcs)))
  subset(obj, cells = keep_bcs)
}


# ---- attach_region_labels() --------------------------------------------
# Add QuPath region columns (neoplastic / hypoxic / ventricular / region
# label) to a Seurat object's metadata, matched by barcode. Does NOT subset —
# it annotates, so downstream code can stratify (e.g. tissue-wide vs
# neoplastic-only) from the same object.
#
# Adds meta.data columns:
#   in_tumor          logical — qupath_in_tumor
#   in_tumor_hypoxic  logical — qupath_in_tumor_hypoxic
#   in_ventricular    logical — qupath_in_ventricular
#   region_label      character — qupath_region_label
#
# Samples without QuPath annotations (no parquet, or columns absent) get
# in_tumor = NA and a warning — callers should treat NA as "region unknown"
# and fall back to tissue-wide analysis only.
#
# Args:
#   obj      — Seurat object (already tissue-masked or not)
#   sid      — sample ID
#   mask_dir — mask directory (holds *_umi_per_bin_{bin}um.parquet)
#   bin_size — bin size in µm (default 16)
#
attach_region_labels <- function(obj, sid, mask_dir, bin_size = 16) {
  region_cols <- c("qupath_in_tumor", "qupath_in_tumor_hypoxic",
                   "qupath_in_ventricular", "qupath_in_hippocampus",
                   "qupath_in_choroid_plexus", "qupath_in_corpus_callosum",
                   "qupath_in_fissure_any", "qupath_region_label")
  out_names   <- c(in_tumor = "qupath_in_tumor",
                   in_tumor_hypoxic = "qupath_in_tumor_hypoxic",
                   in_ventricular = "qupath_in_ventricular",
                   in_hippocampus = "qupath_in_hippocampus",
                   in_choroid_plexus = "qupath_in_choroid_plexus",
                   in_corpus_callosum = "qupath_in_corpus_callosum",
                   in_fissure = "qupath_in_fissure_any",
                   region_label = "qupath_region_label")

  set_na <- function(o) {
    o$in_tumor            <- NA
    o$in_tumor_hypoxic    <- NA
    o$in_ventricular      <- NA
    o$in_hippocampus      <- NA
    o$in_choroid_plexus   <- NA
    o$in_corpus_callosum  <- NA
    o$in_fissure          <- NA
    o$region_label        <- NA_character_
    o
  }

  pq_path <- file.path(mask_dir,
               paste0(sid, "_umi_per_bin_", sprintf("%03d", bin_size),
                      "um.parquet"))
  if (!file.exists(pq_path)) {
    warning("attach_region_labels: parquet not found for ", sid,
            " — region columns set to NA (tissue-wide only).")
    return(set_na(obj))
  }

  avail <- intersect(region_cols,
                     names(arrow::read_parquet(pq_path, as_data_frame = FALSE)))
  if (length(avail) == 0) {
    warning("attach_region_labels: no qupath region columns for ", sid,
            " — set to NA.")
    return(set_na(obj))
  }

  pq  <- arrow::read_parquet(pq_path, col_select = c("barcode", avail))
  idx <- match(colnames(obj), pq$barcode)

  for (meta_name in names(out_names)) {
    src <- out_names[[meta_name]]
    if (src %in% avail) {
      vals <- pq[[src]][idx]
      if (meta_name != "region_label") {
        vals <- as.logical(vals); vals[is.na(vals)] <- FALSE
      }
      obj[[meta_name]] <- vals
    } else {
      obj[[meta_name]] <- if (meta_name == "region_label") NA_character_ else NA
    }
  }
  cat(sprintf("  [region] %s: in_tumor = %s bins | hypoxic = %s | annotated cols: %s\n",
              sid,
              if (all(is.na(obj$in_tumor))) "NA" else sum(obj$in_tumor, na.rm = TRUE),
              if (all(is.na(obj$in_tumor_hypoxic))) "NA" else sum(obj$in_tumor_hypoxic, na.rm = TRUE),
              paste(avail, collapse = ", ")))
  obj
}
