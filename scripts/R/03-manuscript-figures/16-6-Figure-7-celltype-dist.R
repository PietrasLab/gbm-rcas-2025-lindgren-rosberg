# 16-6 — Figure 7: RCTD cell-type weight distributions in TUMOUR tissue (Primary vs Recurrent)
# ============================================================================
# Distribution (violin + box) of each RCTD cell-type weight in the tumour compartment only.
# Tumour tissue = the Tg_positive_expanded mask (tg_membership), the same gate as the spatial maps
# and the niche-frequency panel; Healthy has no tumour bins so it is excluded. Restricting to
# tumour bins is what reveals the Neopl-ACR / reactive-astrocyte rise Primary -> Recurrent that is
# diluted out when normal tissue is included.
#
# Reads the per-bin RCTD weights from results/08-visium-hd-derived
# (env VISIUM_HD_DERIVED_DIR); no re-clustering here.
#
# NO SIGNIFICANCE TESTING ON THE PANEL: with replicate n = 2 primary + 1 recurrent, a by-condition
# test is not warranted, so NO stars/p-values are drawn. The panel is descriptive -- per-sample
# means (white dots) are the honest replicate-level view. A bin-level Cliff's-delta/Wilcoxon table
# is still written to CSV for the record only (pseudoreplicated, n = bins), never shown on the figure.
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-0-source-functions.R"))   # .color_pal
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
DERIVED   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))

SAMP <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")
.cond <- function(s) ifelse(grepl("^HB", s), "Healthy", ifelse(grepl("^TP", s), "Primary", "Recurrent"))
COND_LV  <- c("Healthy", "Primary", "Recurrent")
pal_st   <- .color_pal[["sample_type"]]
COND_PAL <- c(Healthy = pal_st[["Healthy"]], Primary = pal_st[["Primary"]], Recurrent = pal_st[["Recurrent"]])
META <- c("sample_id","condition","barcode","neopl_frac","reactive_astro_frac","in_tumor","x","y")

# Tumour mask (tg_expanded, the Tg_positive_expanded gate) is taken per-bin from the niche
# assignments rather than recomputed from the sf polygon, so this script has no sf dependency.
nb <- read_csv(file.path(DERIVED, "bin-niche-16um.csv.gz"), show_col_types = FALSE) |>
  select(sample_id, barcode, in_tg = tg_expanded)
raw <- lapply(SAMP, function(s) {
  w <- read_csv(file.path(DERIVED, sprintf("rctd-weights-%s-16um.csv.gz", s)), show_col_types = FALSE)
  w$cond <- .cond(s); w
}) |> bind_rows() |> left_join(nb, by = c("sample_id", "barcode"))
raw <- raw[!is.na(raw$in_tg), ]
CT <- setdiff(names(raw), c(META, "in_tg", "cond"))
raw$cond <- factor(raw$cond, COND_LV)
cli::cli_alert_info("{length(CT)} cell types | Tg-tumour bins: {sum(raw$in_tg)} | non-tumour: {sum(!raw$in_tg)}")

long <- raw |> select(sample_id, cond, in_tg, all_of(CT)) |>
  pivot_longer(all_of(CT), names_to = "cell_type", values_to = "weight")
long$cell_type <- factor(long$cell_type, CT)

cliffs_delta <- function(x, y) { n <- length(x) * length(y); (2 * as.numeric(wilcox.test(x, y)$statistic) / n) - 1 }
compute_stats <- function(d, tag) {
  st <- d |> group_by(cell_type) |> group_modify(function(g, k) {
    lv <- levels(droplevels(g$cond)); ng <- length(lv)
    p   <- if (ng > 2) kruskal.test(weight ~ cond, data = g)$p.value else wilcox.test(weight ~ droplevels(cond), data = g)$p.value
    eff <- if (ng == 2) cliffs_delta(g$weight[g$cond == lv[1]], g$weight[g$cond == lv[2]]) else NA_real_
    tibble(test = if (ng > 2) "Kruskal-Wallis" else "Wilcoxon", p = p, cliffs_delta = eff, n_bins = nrow(g))
  }) |> ungroup() |> mutate(star = as.character(cut(p, c(-Inf,1e-4,1e-3,1e-2,5e-2,Inf), c("****","***","**","*","ns"))))
  write_csv(st, file.path(plot_dir, sprintf("Figure-7-celltype-dist-%s-stats.csv", tag))); st
}
dist_plot <- function(d, ttl, tag) {
  compute_stats(d, tag)   # writes the descriptive stats CSV (Cliff's delta / bin-level Wilcoxon) for
                          # the record only; NOT drawn -- n = 2 primary + 1 recurrent, so stars would
                          # overstate a pseudoreplicated (bin-level) comparison.
  psm <- d |> group_by(cell_type, sample_id, cond) |> summarise(mw = mean(weight), .groups = "drop")
  ggplot(d, aes(cond, weight, fill = cond)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey40") +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.2, alpha = 0.6) +
    geom_point(data = psm, aes(cond, mw), shape = 21, size = 1.5, fill = "white", colour = "black", inherit.aes = FALSE) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = COND_PAL, name = NULL) +
    labs(title = ttl, x = NULL, y = "RCTD weight (per bin)",
         subtitle = "white dots = per-sample means. Descriptive only: 2 primary vs 1 recurrent section (no significance testing).") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top", strip.text = element_text(size = 7))
}

# TUMOUR tissue only: Primary vs Recurrent (Healthy has no tumour bins)
pt <- long |> filter(in_tg, cond %in% c("Primary","Recurrent")) |> mutate(cond = droplevels(cond))
ggsave(file.path(plot_dir, "Figure-7-celltype-dist-TUMOR-PTvTR.pdf"),
       dist_plot(pt, "Cell-type weight in tumour tissue (Tg mask): Primary vs Recurrent", "TUMOR-PTvTR"),
       width = 12, height = 9)
# Compact source data: the ~2.4M-row per-bin table goes to the Zenodo deposit, not git. Here we
# ship the box/violin summary per cell type x condition + the per-sample means (the plotted dots).
sd_summary <- pt |> group_by(cell_type, cond) |>
  summarise(n_bins = n(), mean = mean(weight), median = median(weight),
            q25 = quantile(weight, .25), q75 = quantile(weight, .75),
            min = min(weight), max = max(weight), .groups = "drop")
sd_psm <- pt |> group_by(cell_type, sample_id, cond) |> summarise(sample_mean = mean(weight), .groups = "drop")
write_csv(sd_summary, file.path(plot_dir, "Figure-7-celltype-dist-TUMOR-PTvTR-sourcedata.csv"))
write_csv(sd_psm,     file.path(plot_dir, "Figure-7-celltype-dist-TUMOR-PTvTR-persample-means.csv"))
cli::cli_alert_success("16-6 tumour cell-type distributions -> {plot_dir}")
