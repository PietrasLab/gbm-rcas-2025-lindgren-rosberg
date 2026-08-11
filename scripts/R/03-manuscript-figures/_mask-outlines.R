# Mask outlines for spatial overlays, all in the bin pixel frame (== RCTD map x,y,
# validated cor==1 vs parquet pxl_*_in_fullres). Returns a tidy df (x,y,grp,kind):
#   kind = "Tg"      : Tg_positive_expanded (spatial-masks tg_overlay + 24um buffer)
#   kind = "Tissue"  : QuPath Tissue polygon (a5 export) -> tissue border
#   kind = "Exclude" : QuPath Fold/Exclude_HE/Blur -> excluded regions
# a5 geojson coords are in the same pixel frame (rasterized upstream onto
# pxl_*_in_fullres), so they overlay the RCTD maps directly.
suppressPackageStartupMessages({ library(jsonlite); library(sf); library(dplyr); library(arrow); library(RANN) })
sf::sf_use_s2(FALSE)
# Publication repo: the Tg overlay is STAGED into results/06-2-tg-masks by
# scripts/R/01-import-filter-annotation/06-2-visium-hd-tg-mask.R (gitignored);
# the manual QuPath annotations + bin capture-footprint parquet are tracked in
# references/visium-hd-masks/.
.MO_PROJ <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
.MO_TG   <- file.path(.MO_PROJ, "results/06-2-tg-masks")
.MO_ANN  <- file.path(.MO_PROJ, "references/visium-hd-masks/qupath-annotations")
.MO_MASK <- Sys.getenv("VISIUM_HD_MASK_DIR",
  unset = file.path(.MO_PROJ, "references/visium-hd-masks/masks-workspace"))

.mo_read_qp <- function(path) {
  raw <- paste(readLines(path, warn = FALSE), collapse = "")
  if (startsWith(trimws(raw), "[")) raw <- paste0('{"type":"FeatureCollection","features":', raw, "}")
  tmp <- tempfile(fileext = ".geojson"); writeLines(raw, tmp); on.exit(unlink(tmp))
  suppressWarnings(sf::st_read(tmp, quiet = TRUE))
}
.mo_label <- function(gj) {
  cls <- rep(NA_character_, nrow(gj))
  if ("classification" %in% names(gj)) {
    rc <- gj$classification
    if (is.list(rc)) cls <- vapply(rc, function(x) if (is.null(x$name)) NA_character_ else as.character(x$name), character(1))
    else { cls <- as.character(rc); cls <- gsub('.*"name"[[:space:]]*:[[:space:]]*"([^"]+)".*', "\\1", cls) }
  }
  nm <- if ("name" %in% names(gj)) as.character(gj$name) else rep(NA_character_, nrow(gj))
  cls[is.na(cls) | cls == "" | cls == "NA"] <- nm[is.na(cls) | cls == "" | cls == "NA"]
  sub("_R$", "", cls)
}
# drop polygon parts below a minimum area (removes tiny fold "bubbles")
.mo_drop_small <- function(g, min_area) {
  if (is.null(g) || length(g) == 0 || all(sf::st_is_empty(g))) return(g)
  pp <- suppressWarnings(sf::st_cast(sf::st_make_valid(g), "POLYGON"))
  pp <- pp[as.numeric(sf::st_area(pp)) >= min_area]
  if (length(pp) == 0) return(NULL)
  sf::st_union(pp)
}
.mo_df <- function(geom, sid, kind) {
  if (is.null(geom) || length(geom) == 0 || all(sf::st_is_empty(geom))) return(NULL)
  co <- as.data.frame(sf::st_coordinates(geom))
  gcols <- intersect(c("L3", "L2", "L1"), names(co))
  tibble(sample_id = sid, kind = kind, x = co$X, y = co$Y,
         grp = paste(kind, do.call(paste, c(co[gcols], sep = "_")), sep = "_"))
}
# clip an sfc to another, keeping only polygonal parts (intersections can yield
# GEOMETRYCOLLECTIONs); returns NULL if empty. All geoms forced to NA CRS.
.mo_clip <- function(g, by) {
  if (is.null(g) || is.null(by) || length(g) == 0 || all(sf::st_is_empty(g))) return(NULL)
  g  <- sf::st_make_valid(sf::st_set_crs(g, NA))
  by <- sf::st_make_valid(sf::st_set_crs(sf::st_union(by), NA))
  r  <- suppressWarnings(sf::st_intersection(g, by))
  if (length(r) == 0 || all(sf::st_is_empty(r))) return(NULL)
  r <- suppressWarnings(sf::st_collection_extract(r, "POLYGON"))
  if (length(r) == 0 || all(sf::st_is_empty(r))) return(NULL)
  r
}
mask_outlines <- function(sid) {
  out <- list()
  # ---- capture-area footprint = bbox of the captured bins (parquet). EVERYTHING
  # is clipped to this so no outline is drawn beyond the Visium capture area
  # (the tissue annotation can be larger than capture -> blank sketch otherwise).
  pqf <- file.path(.MO_MASK, sprintf("%s_umi_per_bin_016um.parquet", sid))
  foot <- NULL; pitch <- NA_real_
  if (file.exists(pqf)) {
    xy <- as.matrix(arrow::read_parquet(pqf, col_select = c("pxl_col_in_fullres", "pxl_row_in_fullres")))
    xy <- xy[complete.cases(xy), ]; pitch <- median(RANN::nn2(xy, xy, k = 2)$nn.dists[, 2])
    rx <- range(xy[, 1]); ry <- range(xy[, 2])
    foot <- sf::st_sfc(sf::st_polygon(list(rbind(c(rx[1], ry[1]), c(rx[2], ry[1]),
                                                 c(rx[2], ry[2]), c(rx[1], ry[2]), c(rx[1], ry[1])))))
  }
  # ---- Tissue + excluded (QuPath a5), clipped to capture footprint -----------
  cand <- Sys.glob(file.path(.MO_ANN, paste0(sid, "_qupath_annotations.geojson")))
  tissue_g <- NULL
  if (length(cand)) {
    gj <- .mo_read_qp(cand[which.max(file.info(cand)$mtime)]); gj$lab <- .mo_label(gj)
    tissue_g <- .mo_clip(sf::st_union(gj[gj$lab == "Tissue", ]), foot)
    out$tissue  <- .mo_df(tissue_g, sid, "Tissue")
    # excluded (Fold/Blur/etc): drop tiny bubbles (esp. TR_03d_990 folds) that plot
    # as "commas"; keep polygons >= ~9 bins of area.
    excl <- .mo_clip(sf::st_union(gj[gj$lab %in% c("Fold","Exclude_HE","ExludeHE","Blur","Exclude_analysis","Exclude_segmentation"), ]), foot)
    if (!is.null(excl) && !is.na(pitch)) excl <- .mo_drop_small(excl, min_area = 9 * pitch^2)
    out$exclude <- .mo_df(excl, sid, "Exclude")
  }
  # ---- Tg tumor (overlay + 24um buffer), clipped to TISSUE (never outside the
  # tissue border) and to the capture footprint ------------------------------
  gj_tg <- file.path(.MO_TG, sprintf("%s_tg_overlay.geojson", sid))
  if (file.exists(gj_tg) && !is.na(pitch)) {
    x <- fromJSON(gj_tg, simplifyVector = FALSE); feat <- if (!is.null(x$features)) x$features[[1]] else x[[1]]
    polys <- lapply(feat$geometry$coordinates, function(p) lapply(p, function(r) do.call(rbind, lapply(r, function(q) c(q[[1]], q[[2]])))))
    tg <- sf::st_buffer(sf::st_sfc(sf::st_multipolygon(polys)), dist = (24 / 16) * pitch)
    tg <- .mo_clip(tg, if (!is.null(tissue_g)) tissue_g else foot)   # within tissue (or capture if no tissue)
    out$tg <- .mo_df(tg, sid, "Tg")
  }
  bind_rows(out)
}

# per-bin logical: is bin (x,y in the fullres pixel frame) INSIDE the
# Tg_positive_expanded tumor polygon (the SAME polygon mask_outlines draws). Used
# to define 'tumor bins' independently of RCTD / the scored signatures. Samples
# with no Tg overlay (healthy) -> all FALSE.
tg_membership <- function(sid, x, y) {
  gj_tg <- file.path(.MO_TG, sprintf("%s_tg_overlay.geojson", sid))
  if (!file.exists(gj_tg)) return(rep(FALSE, length(x)))
  pqf <- file.path(.MO_MASK, sprintf("%s_umi_per_bin_016um.parquet", sid))
  xy <- as.matrix(arrow::read_parquet(pqf, col_select = c("pxl_col_in_fullres", "pxl_row_in_fullres")))
  xy <- xy[complete.cases(xy), ]; pitch <- median(RANN::nn2(xy, xy, k = 2)$nn.dists[, 2])
  jj <- jsonlite::fromJSON(gj_tg, simplifyVector = FALSE); feat <- if (!is.null(jj$features)) jj$features[[1]] else jj[[1]]
  polys <- lapply(feat$geometry$coordinates, function(p) lapply(p, function(r) do.call(rbind, lapply(r, function(q) c(q[[1]], q[[2]])))))
  tg <- sf::st_make_valid(sf::st_set_crs(sf::st_buffer(sf::st_sfc(sf::st_multipolygon(polys)), dist = (24 / 16) * pitch), NA))
  pts <- sf::st_set_crs(sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y")), NA)
  as.logical(sf::st_within(pts, tg, sparse = FALSE)[, 1])
}

# ggplot layers for the CANONICAL spatial-plot overlay. CONVENTION: EVERY spatial
# plot gets the tissue border (solid grey) + excluded regions (dashed grey). Plots
# where the tumor region is meaningful ALSO get the Tg_positive_expanded outline
# (pink = neoplastic palette colour), toggled with tg=TRUE. Faceted callers must
# give the outline df the same sample_id factor levels as the point data.
MASK_TG_COL <- "#08306B"   # dark blue for the Tg tumor outline (v32: thin, 50% alpha)
mask_layers <- function(outlines, tissue = TRUE, exclude = TRUE, tg = TRUE,
                        lwd_tissue = 0.3, lwd_tg = 0.2, lwd_exclude = 0.15, alpha_tg = 0.5) {   # Tg SOLID; excluded thinner
  ll <- list()
  if (tissue)  ll <- c(ll, ggplot2::geom_path(data = dplyr::filter(outlines, kind == "Tissue"),
                       ggplot2::aes(x, y, group = grp), color = "grey55", linewidth = lwd_tissue, inherit.aes = FALSE))
  if (exclude) ll <- c(ll, ggplot2::geom_path(data = dplyr::filter(outlines, kind == "Exclude"),
                       ggplot2::aes(x, y, group = grp), color = "grey30", linewidth = lwd_exclude, linetype = "22", inherit.aes = FALSE))
  if (tg)      ll <- c(ll, ggplot2::geom_path(data = dplyr::filter(outlines, kind == "Tg"),
                       ggplot2::aes(x, y, group = grp), color = MASK_TG_COL, linewidth = lwd_tg, alpha = alpha_tg, inherit.aes = FALSE))
  ll
}

# ============================================================================
# spatial_bin_grid() — THE canonical spatial map for the visium pub cohort.
# Draws every bin as its TRUE 16um square (geom_tile sized to the per-sample bin
# pitch) with coord_fixed per sample (patchwork), the fig1d continuous fill, the
# mask outlines (tissue/excluded/Tg), and the legend at the bottom. Use this for
# ALL continuous spatial maps so the look is identical everywhere.
#   dat       : data frame with sample_id, x, y, and `fill_col` (ALL bins, no
#               downsampling — tiles must cover the tissue).
#   outlines  : bind_rows(lapply(samp_ord, mask_outlines)), sample_id factored.
#   fill_scale: defaults to fig1d_fill(); pass ggplot2::scale_fill_manual(...) for
#               a DISCRETE map (e.g. dominant cell type from .color_pal).
# Requires ggplot2 + patchwork + RANN loaded by the caller.
# ============================================================================
spatial_bin_grid <- function(dat, fill_col, samp_ord, outlines, title = NULL,
                             fill_scale = fig1d_fill(limits = c(0, 1), name = fill_col),
                             tg = TRUE, title_col = "black") {
  if (!"bin_px" %in% names(dat)) {
    pt <- dat |> dplyr::group_by(sample_id) |>
      dplyr::group_modify(function(.x, .y) {
        m <- as.matrix(.x[, c("x", "y")]); data.frame(bin_px = median(RANN::nn2(m, m, k = 2)$nn.dists[, 2]))
      }) |> dplyr::ungroup()
    dat <- dplyr::left_join(dat, pt, by = "sample_id")
  }
  one <- function(s, keep_leg) {
    d <- dplyr::filter(dat, sample_id == s); o <- dplyr::filter(outlines, sample_id == s)
    p <- ggplot2::ggplot(d, ggplot2::aes(x, y, fill = .data[[fill_col]])) +
      ggplot2::geom_tile(ggplot2::aes(width = bin_px, height = bin_px)) +
      mask_layers(o, tg = tg) + fill_scale +
      ggplot2::scale_y_reverse() + ggplot2::coord_fixed() + ggplot2::labs(title = s) +
      theme_fig1d(10) + ggplot2::theme(aspect.ratio = NULL,
        plot.title = ggplot2::element_text(size = 9, hjust = 0.5))
    if (keep_leg && is.factor(dat[[fill_col]])) {
      # DISCRETE legend fix: the single legend is taken from this (first) panel, which may
      # not contain every level (e.g. the healthy section has no tumour niches). drop=FALSE
      # + override.aes(fill=) does NOT colour absent-level keys (verified). Seeding one
      # invisible (alpha=0) tile per level makes every key render with its scale colour;
      # override.aes(alpha=1) restores key opacity.
      lv   <- levels(dat[[fill_col]])
      seed <- data.frame(x = d$x[1], y = d$y[1], bin_px = d$bin_px[1])[rep(1L, length(lv)), , drop = FALSE]
      seed[[fill_col]] <- factor(lv, lv)
      p <- p + ggplot2::geom_tile(data = seed, ggplot2::aes(width = bin_px, height = bin_px),
                                  alpha = 0, show.legend = TRUE) +
               ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
    if (!keep_leg) p <- p + ggplot2::guides(fill = "none")   # see below
    p
  }
  # ONE legend only: each panel is filtered to its sample, so the 4 fill guides differ
  # (different niches present) and ggplot2 4.0 + patchwork 1.3 render each -> duplicated
  # legends. Keep the legend on the FIRST panel only (drop=FALSE there lists ALL levels)
  # and collect that single guide into a dedicated guide_area() row. Sizing = theme_fig1d.
  plots <- Map(function(s, i) one(s, i == 1L), samp_ord, seq_along(samp_ord))
  (patchwork::wrap_plots(plots, nrow = 1) / patchwork::guide_area()) +
    patchwork::plot_layout(guides = "collect", heights = c(1, 0.1)) +
    patchwork::plot_annotation(title = title,
      theme = ggplot2::theme(legend.position = "bottom",
                             plot.title = ggplot2::element_text(face = "bold", colour = title_col)))
}
