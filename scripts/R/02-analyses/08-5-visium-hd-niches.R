# 08-5 — Visium HD de-novo spatial composition niches (joint k-means)
# ============================================================================
# Unsupervised tissue-organisation view complementing the reference-anchored
# maps. For each in-tissue 16um bin, the KNN-mean (k = 30) of the 08-3 RCTD
# weight vectors is computed per section, then a single JOINT k-means (12
# niches) over ALL bins of all 4 sections gives shared niche labels comparable
# across Healthy / Primary / Recurrent. Tg_positive_expanded (tg_membership) is
# carried as a per-bin STRATIFIER (not a clustering gate).
# Consumed by the Figure 7 niche script (16-3); plots are rendered there.
#
# INPUTS: 08-3 RCTD weights (run 08-3 first); Tg overlays via _mask-outlines.R
# OUTPUT (results/08-visium-hd-derived/):
#   bin-niche-16um.csv.gz + niche annotation / enrichment / frequency CSVs
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(RANN); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(fig_dir, "_mask-outlines.R"))   # tg_membership
out_dir <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))

BIN_SIZE <- 16L; K_NEIGH <- 30L; N_NICHES <- 12L
SAMP_ORD <- c("HB_hm01", "TP_1083", "TP_987", "TR_03d_990")
.cond <- function(s) ifelse(grepl("^HB", s), "Healthy", ifelse(grepl("^TP", s), "Primary", "Recurrent"))
COND_LV <- c("Healthy","Primary","Recurrent")
set.seed(1)

META_COLS <- c("sample_id","condition","barcode","neopl_frac","reactive_astro_frac","in_tumor","x","y","tg_expanded","section")
.probe <- read_csv(file.path(out_dir, sprintf("rctd-weights-%s-%dum.csv.gz", SAMP_ORD[2], BIN_SIZE)), n_max = 1, show_col_types = FALSE)
CT <- setdiff(names(.probe), META_COLS)
cli::cli_alert_info("cell types ({length(CT)}): {paste(CT, collapse=', ')}")

read_weights <- function(sid) {
  d <- read_csv(file.path(out_dir, sprintf("rctd-weights-%s-%dum.csv.gz", sid, BIN_SIZE)), show_col_types = FALSE)
  d$tg_expanded <- tg_membership(sid, d$x, d$y); d$cond <- .cond(sid)
  cli::cli_alert("{sid}: {nrow(d)} bins | Tg+ {sum(d$tg_expanded)} ({round(100*mean(d$tg_expanded),1)}%)")
  d
}
neighborhood_profiles <- function(df, k = K_NEIGH) {
  W <- as.matrix(df[, CT, drop = FALSE]); xy <- as.matrix(df[, c("x","y")]); k <- min(k, nrow(df))
  nn <- RANN::nn2(xy, xy, k = k)$nn.idx
  P <- matrix(0, nrow(W), ncol(W), dimnames = list(rownames(W), CT))
  for (j in seq_len(k)) P <- P + W[nn[, j], , drop = FALSE]
  P / k
}

# ---- joint composition niches over all 4 sections ---------------------------
cli::cli_h2("Composition niches, joint over {paste(SAMP_ORD, collapse=', ')}")
dat <- lapply(SAMP_ORD, function(sid) { d <- read_weights(sid); list(df = d, P = neighborhood_profiles(d)) })
Pall <- do.call(rbind, lapply(dat, `[[`, "P"))
km <- kmeans(Pall, centers = N_NICHES, iter.max = 100, nstart = 10)
ord <- order(table(km$cluster), decreasing = TRUE); relab <- setNames(seq_along(ord), ord)
niche <- factor(relab[as.character(km$cluster)], levels = seq_len(N_NICHES))
bind <- do.call(rbind, lapply(dat, `[[`, "df")); bind$niche <- niche
bind$sample_id <- factor(bind$sample_id, SAMP_ORD); bind$cond <- factor(bind$cond, COND_LV)

# niche x cell-type enrichment (z per cell type across niches)
enr <- bind |> group_by(niche) |> summarise(across(all_of(CT), mean), .groups = "drop") |>
  pivot_longer(all_of(CT), names_to = "cell_type", values_to = "mean_weight") |>
  group_by(cell_type) |> mutate(z = as.numeric(scale(mean_weight))) |> ungroup()
enr$cell_type <- factor(enr$cell_type, rev(CT))

# niche frequency by condition + Tg-stratified
bind$tg_stratum <- ifelse(bind$tg_expanded, "tumor (Tg+)", "non-tumor (Tg-)")
freq <- bind |> count(sample_id, cond, niche, tg_stratum) |> group_by(sample_id) |> mutate(frac = n/sum(n)) |> ungroup()

# niche annotation + source data
ann <- bind |> group_by(niche) |> summarise(n_bins = n(), mean_neopl_frac = mean(neopl_frac, na.rm=TRUE),
  mean_reactive_astro = mean(reactive_astro_frac, na.rm=TRUE), frac_tg_expanded = mean(tg_expanded, na.rm=TRUE),
  frac_in_tumor = mean(as.logical(in_tumor), na.rm=TRUE), .groups = "drop")
write_csv(ann,  file.path(out_dir, "niche-annotation.csv"))
write_csv(enr,  file.path(out_dir, "niche-celltype-enrichment.csv"))
write_csv(freq, file.path(out_dir, "niche-frequency.csv"))
write_csv(bind[, c("sample_id","cond","barcode","x","y","niche","tg_expanded","neopl_frac","reactive_astro_frac","in_tumor")],
          file.path(out_dir, "bin-niche-16um.csv.gz"))
print(as.data.frame(ann))
cli::cli_alert_success("08-5 composition niches -> {out_dir} (bin-niche-16um.csv.gz + niche-*.csv)")
