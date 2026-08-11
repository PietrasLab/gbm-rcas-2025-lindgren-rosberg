# 16-4 — Figure 7: neoplastic-panel split-triangle correlation (Primary vs Recurrent, 16um)
# ============================================================================
# Split-triangle Spearman corrplot of the neoplastic-compartment signatures (our 6 + Neftel /
# Richards / Nomura references) across the Tg-tumor 16um bins: upper-left triangle = Primary
# (TP_987 + TP_1083 pooled), lower-right = Recurrent (TR_03d_990); rows ordered by hclust of the
# Recurrent matrix, source annotation bar, anti-diagonal divider. Reads the per-bin signature
# scores + panel meta from results/08-visium-hd-derived (env VISIUM_HD_DERIVED_DIR);
# scoring is upstream (not re-run here).
# NB bin-level n (~1e5) -> p-values meaningless; magnitude is the signal.
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ComplexHeatmap); library(circlize); library(grid); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
DERIVED   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))
META <- c("sample_id","barcode","x","y","in_tg","neopl_frac")
SAMP <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990"); PT <- c("TP_987","TP_1083")

meta <- read_csv(file.path(DERIVED, "signature-panel-meta.csv"), show_col_types = FALSE)
scores <- setNames(lapply(SAMP, function(s)
  read_csv(file.path(DERIVED, sprintf("signature-scores-%s-16um.csv.gz", s)), show_col_types = FALSE)), SAMP)
present <- Reduce(intersect, lapply(scores, function(d) setdiff(names(d), META)))
get_bins <- function(sids, tissue = c("all","tumor","nontumor")) {
  tissue <- match.arg(tissue); d <- bind_rows(scores[sids])
  if (tissue == "tumor") d <- d[d$in_tg %in% TRUE, ]; if (tissue == "nontumor") d <- d[d$in_tg %in% FALSE, ]; d
}
SRC_PAL <- c("Ours"="#D5006D","Neftel 2019"="#1F78B4","Richards 2021"="#33A02C","Nomura 2025"="#FF7F00","Hallmark"="#984EA3")

split_corr <- function(dA, dB, labA, labB, sigs, scope, fname, width = 8.4, height = 7.6) {
  cA <- suppressWarnings(cor(as.matrix(dA[, sigs]), method="spearman", use="pairwise.complete.obs"))
  cB <- suppressWarnings(cor(as.matrix(dB[, sigs]), method="spearman", use="pairwise.complete.obs"))
  cA[is.na(cA)] <- 0; cB[is.na(cB)] <- 0
  ord <- hclust(as.dist(1 - cB), method="complete")$order
  N <- length(ord); ro <- ord; co <- rev(ord)
  D <- matrix(NA_real_, N, N, dimnames = list(sigs[ro], sigs[co]))
  for (i in 1:N) for (j in 1:N) { if (i + j == N + 1) next; src <- if (i + j < N + 1) cA else cB; D[i,j] <- src[ro[i], co[j]] }
  m <- meta[match(sigs[ro], meta$label), ]
  la <- rowAnnotation(df = data.frame(source = m$source), col = list(source = SRC_PAL),
                      show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 6),
                      simple_anno_size = unit(3, "mm"))
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
  ht <- Heatmap(D, name = "Spearman r", col = col_fun, na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", left_annotation = la,
    row_names_gp = gpar(fontsize = 6.5), column_names_gp = gpar(fontsize = 6.5),
    column_title = sprintf("UL = %s  |  LR = %s   (%s)", labA, labB, scope),
    column_title_gp = gpar(fontsize = 9, fontface = "bold"), rect_gp = gpar(col = "white", lwd = 0.3))
  pdf(file.path(plot_dir, fname), width = width, height = height)
  draw(ht, merge_legend = TRUE)
  decorate_heatmap_body("Spearman r", grid.lines(c(0,1), c(0,1), gp = gpar(col = "black", lwd = 1.2)))
  invisible(dev.off())
  write_csv(tibble::rownames_to_column(as.data.frame(D), "signature"),
            file.path(plot_dir, sub("\\.pdf$", "-sourcedata.csv", fname)))
  cli::cli_alert_success("{fname}  ({length(sigs)} sigs; nA={nrow(dA)}, nB={nrow(dB)})")
}

# NEOPLASTIC cells (Tg-tumor bins), neoplastic-compartment sigs only.
# width x1.3 vs the upstream default so the labels are legible.
Sneo <- intersect(meta$label[meta$compartment == "neoplastic"], present)
split_corr(get_bins(PT,"tumor"), get_bins("TR_03d_990","tumor"),
           "Primary (TP_987+TP_1083)", "Recurrent (TR_03d_990)", Sneo,
           "Tg-tumor bins, neoplastic panel (ours + refs)", "Figure-7-corr-neoplastic.pdf",
           width = 8.4 * 1.3)
cli::cli_alert_success("16-4 neoplastic split-triangle corrplot -> {plot_dir}")
