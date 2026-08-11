# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(glue)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
})
conflicts_prefer(dplyr::filter, dplyr::mutate, dplyr::summarise, .quiet = TRUE)

# Canonical manuscript palette (.color_pal): always use it for sample_type / Levels.
source("./scripts/R/00-0-source-functions.R")

# 11-2 — Supplementary Figure 2 (Figure 2 supplement): per-sample composition + Wilcoxon
# Outputs go in the matching main-figure folder (manuscript-figures/figure-2/), named
# Figure-S2-* per the repo convention.
#
# Consolidated supplemental: per-sample cell-type proportions with pairwise
# Wilcoxon (manual, BH-corrected across cell types) for each resolution level.
# Reads metadata only from the filtered Seurat object.
#   A — Level_1 (Neoplastic / Non-Neoplastic), % of all cells   [-> main Fig 2G]
#   B — Level_3 all cell types, % of all cells                  [supplement]
#   C — Astrocyte subclusters (TE/NT/R), % of astrocytes        [-> main Fig 3B]
#   D — Neoplastic clusters (Level_4), % of neoplastic          [-> main Fig 4F]
#   E — Neopl-ACR, % of neoplastic                              [-> main Fig 5D]
# Plot brackets show the BH-corrected significance (p.adj.signif from run_wilcoxon,
# drawn with stat_pvalue_manual; only significant brackets shown). The full pairwise
# table (raw p, BH p.adj) is in the stats-*.csv supplementary tables.

obj_path <- "data/processed/seurat/seurat_flex_filtered_v1.0.rds"
plot.dir <- "./manuscript-figures/figure-2"      # Supp Fig 2 lives with Figure 2
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
pfx <- "Figure-S2-composition"                    # filename prefix

COND_LEVELS <- c("Healthy", "Primary", "Recurrent")
COND_COLORS <- .color_pal[["sample_type"]]   # canonical manuscript palette
comps_all   <- list(c("Healthy", "Primary"), c("Primary", "Recurrent"), c("Healthy", "Recurrent"))
comps_tumor <- list(c("Primary", "Recurrent"))
l3_types    <- c("Astrocyte","OPC/COP/OLG","Neural","Ependymal","Microglia",
                 "Macrophage","Dendritic","Neutrophil","NKTB","Choroid",
                 "Endothelial","Mural","Fibroblast","Neoplastic")
astro_types <- c("Astrocyte TE", "Astrocyte NT", "Astrocyte R")
neopl_clusters <- c("Neopl-Bulk","Neopl-CC-I","Neopl-CC-II","Neopl-CC-III",
                    "Neopl-OPC","Neopl-COP","Neopl-NC","Neopl-ACR",
                    "Neopl-ECM","Neopl-RNA-low")

# ---- Load metadata ------------------------------------------------------------
cli_alert_info("Loading metadata from {obj_path}")
mdata <- readRDS(obj_path)@meta.data
mdata$condition <- factor(mdata$sample_type, levels = COND_LEVELS)

# ---- Helpers ------------------------------------------------------------------
per_sample_props <- function(df, group_col, filter_col = NULL, filter_vals = NULL) {
  if (!is.null(filter_col)) df <- df[df[[filter_col]] %in% filter_vals, ]
  raw <- df |>
    group_by(sample_id, condition, .data[[group_col]]) |>
    summarise(n = n(), .groups = "drop") |>
    group_by(sample_id, condition) |>
    mutate(pct = n / sum(n) * 100) |> ungroup()
  all_samples <- unique(df[, c("sample_id", "condition")])
  all_groups  <- unique(as.character(df[[group_col]]))
  full_grid   <- tidyr::crossing(all_samples, !!group_col := all_groups)
  colnames(full_grid)[3] <- group_col
  raw |>
    right_join(full_grid, by = c("sample_id", "condition", group_col)) |>
    mutate(n = coalesce(n, 0L), pct = coalesce(pct, 0)) |>
    mutate(!!group_col := factor(.data[[group_col]]))
}

run_wilcoxon <- function(props, group_col, comparisons_list) {
  clusters <- unique(props[[group_col]])
  res <- lapply(clusters, function(cl) {
    d <- props[props[[group_col]] == cl, ]
    lapply(comparisons_list, function(cmp) {
      x <- d$pct[d$condition == cmp[1]]; y <- d$pct[d$condition == cmp[2]]
      x <- x[!is.na(x)]; y <- y[!is.na(y)]
      if (length(x) < 2 || length(y) < 2) return(NULL)
      if (sum(x) == 0 && sum(y) == 0) return(NULL)
      wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
      r <- data.frame(group1 = cmp[1], group2 = cmp[2], n1 = length(x), n2 = length(y),
                      statistic = wt$statistic, p = wt$p.value)
      r[[group_col]] <- cl; r
    })
  })
  out <- bind_rows(Filter(Negate(is.null), unlist(res, recursive = FALSE)))
  out <- out[!is.nan(out$p), ]
  if (nrow(out) == 0) return(out)
  out |> mutate(p.adj = p.adjust(p, method = "BH"),
                p.adj.signif = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                                         p.adj < 0.05 ~ "*", p.adj < 0.1 ~ ".", TRUE ~ "ns"))
}

# BH-corrected significance brackets: place per-facet y positions from run_wilcoxon output
# (only significant comparisons shown), then draw with stat_pvalue_manual.
bracket_df <- function(stats, props, group_col) {
  if (is.null(stats) || nrow(stats) == 0) return(NULL)
  sig <- stats[stats$p.adj.signif != "ns", , drop = FALSE]
  if (nrow(sig) == 0) return(NULL)
  ymax <- props |> group_by(.data[[group_col]]) |>
    summarise(.ymax = max(pct, na.rm = TRUE), .groups = "drop")
  sig |> left_join(ymax, by = group_col) |>
    group_by(.data[[group_col]]) |>
    mutate(y.position = .ymax * (1.04 + 0.10 * (row_number() - 1))) |>
    ungroup()
}

panel <- function(props, group_col, comps, ylab, nrow = 1, stats = NULL) {
  p <- ggplot(props, aes(condition, pct, color = condition, fill = condition)) +
    geom_boxplot(alpha = 0.25, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.85) +
    scale_color_manual(values = COND_COLORS) + scale_fill_manual(values = COND_COLORS) +
    facet_wrap(vars(.data[[group_col]]), scales = "free_y", nrow = nrow) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none",
          strip.text = element_text(size = 8), panel.grid.major.x = element_blank())
  br <- bracket_df(stats, props, group_col)
  if (!is.null(br)) {
    p <- p + ggpubr::stat_pvalue_manual(br, label = "p.adj.signif", tip.length = 0.01,
                                        size = 3, bracket.size = 0.3, inherit.aes = FALSE)
  }
  p
}

# ---- Panels + stats -----------------------------------------------------------
props_L1    <- per_sample_props(mdata, "Level_1", "Level_1", c("Neoplastic","Non-Neoplastic"))
props_L3    <- per_sample_props(mdata, "Level_3", "Level_3", l3_types) |>
  mutate(Level_3 = factor(Level_3, levels = l3_types))
props_astro <- per_sample_props(mdata, "Level_3AC", "Level_3AC", astro_types) |>
  mutate(Level_3AC = factor(Level_3AC, levels = astro_types))
props_neopl <- per_sample_props(mdata, "Level_4", "Level_4", neopl_clusters) |>
  mutate(Level_4 = factor(Level_4, levels = neopl_clusters))
props_nacr  <- props_neopl |> filter(Level_4 == "Neopl-ACR")

st_L1    <- run_wilcoxon(props_L1,    "Level_1",   comps_all)
st_L3    <- run_wilcoxon(props_L3,    "Level_3",   comps_all)
st_astro <- run_wilcoxon(props_astro, "Level_3AC", comps_all)
st_neopl <- run_wilcoxon(props_neopl, "Level_4",   comps_tumor)
st_nacr  <- run_wilcoxon(props_nacr,  "Level_4",   comps_tumor)

write_csv(st_L1,    file.path(plot.dir, glue("{pfx}-stats-Level1.csv")))
write_csv(st_L3,    file.path(plot.dir, glue("{pfx}-stats-Level3-all.csv")))
write_csv(st_astro, file.path(plot.dir, glue("{pfx}-stats-astro.csv")))
write_csv(st_neopl, file.path(plot.dir, glue("{pfx}-stats-neoplastic.csv")))
write_csv(st_nacr,  file.path(plot.dir, glue("{pfx}-stats-NeoplACR.csv")))

# panels carry BH-corrected significance brackets (stats = run_wilcoxon output)
pA <- panel(props_L1,    "Level_1",   comps_all,   "% of all cells",       nrow = 1, stats = st_L1)
pB <- panel(props_L3,    "Level_3",   comps_all,   "% of all cells",       nrow = 2, stats = st_L3)
pC <- panel(props_astro, "Level_3AC", comps_all,   "% of astrocytes",      nrow = 1, stats = st_astro)
pD <- panel(props_neopl, "Level_4",   comps_tumor, "% of neoplastic",      nrow = 2, stats = st_neopl)
pE <- panel(props_nacr,  "Level_4",   comps_tumor, "% of neoplastic",      stats = st_nacr)

ggsave(file.path(plot.dir, glue("{pfx}-Level1.pdf")),     pA, width = 5,  height = 4)
ggsave(file.path(plot.dir, glue("{pfx}-Level3-all.pdf")), pB, width = 12, height = 6)
ggsave(file.path(plot.dir, glue("{pfx}-astro.pdf")),      pC, width = 6,  height = 4)
ggsave(file.path(plot.dir, glue("{pfx}-neoplastic.pdf")), pD, width = 11, height = 6)
ggsave(file.path(plot.dir, glue("{pfx}-NeoplACR.pdf")),   pE, width = 3.2, height = 4)

write_csv(pA$data, file.path(plot.dir, glue("{pfx}-Level1-sourcedata.csv")))
write_csv(pB$data, file.path(plot.dir, glue("{pfx}-Level3-all-sourcedata.csv")))
write_csv(pC$data, file.path(plot.dir, glue("{pfx}-astro-sourcedata.csv")))
write_csv(pD$data, file.path(plot.dir, glue("{pfx}-neoplastic-sourcedata.csv")))
write_csv(pE$data, file.path(plot.dir, glue("{pfx}-NeoplACR-sourcedata.csv")))

cli_alert_success("Supp Fig 2 (composition) panels + stats -> {plot.dir}/ ({pfx}-*)")
