# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})
conflicts_prefer(dplyr::filter, dplyr::select, .quiet = TRUE)
source("./scripts/R/00-0-source-functions.R")   # canonical .color_pal

# 10-4 — Supp Fig 1 (Figure 1 supplement), REVIEWERS-ONLY panels: RCAS reporter
# transgene vs inferred CNA as a neoplastic marker (revision).
# Shows the reporter transgene is a clean, robust neoplastic
# marker and superior to inferCNV CNA in this RTK-driven (chromosomally stable) model:
# (1) transgene counts bimodal, (2) specific to neoplastic cells, (3) non-neoplastic
# transgene is mostly ambient (decontX), (4) head-to-head AUC vs inferCNV CNA.
# Labelled as Supp Fig 1 but intended for the response only (not the final manuscript).

obj_path  <- "data/processed/seurat/seurat_flex_filtered_v1.0.rds"
# inferCNV object is a revision-exploration product; referenced via data/external.
icnv_path <- "data/external/c2-infercnv-full-object.rds"
plot.dir  <- "./manuscript-figures/figure-1"
pfx       <- "Figure-S1-transgene"     # reviewers-only panels of Supp Fig 1 (tracked in docs)

cols_neopl <- c(`FALSE` = unname(.color_pal[["Level_1"]][["Non-Neoplastic"]]),
                `TRUE`  = unname(.color_pal[["Level_1"]][["Neoplastic"]]))
lab_neopl  <- c("non-neoplastic", "neoplastic")

auc <- function(score, pos) {
  ok <- !is.na(score); score <- score[ok]; pos <- pos[ok]
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- object + transgene quantification ----------------------------------------
cli_alert_info("Loading object + quantifying transgene...")
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(SeuratObject::JoinLayers(obj), error = function(e) obj)
counts <- SeuratObject::LayerData(obj, layer = "counts")

tg_feats <- grep("Tg-RFP|Tg-hPDGFB|CUSTOMPROBE-RFP|hPDGFB", rownames(obj),
                 value = TRUE, ignore.case = TRUE)
cli_alert_info("transgene features: {paste(tg_feats, collapse=', ')}")
obj$tg_sum   <- Matrix::colSums(counts[tg_feats, , drop = FALSE])
obj$tg_log   <- log1p(obj$tg_sum)
obj$is_neopl <- grepl("^Neopl", as.character(obj$Level_4))
md <- obj@meta.data %>% rownames_to_column("barcode") %>%
  select(barcode, Level_4, is_neopl, tg_sum, tg_log)

# 1) bimodality
p_hist <- ggplot(md, aes(tg_log, fill = is_neopl)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = cols_neopl, labels = lab_neopl, name = NULL) +
  labs(title = "RCAS transgene counts are bimodal", subtitle = NULL,
       x = "log1p(transgene UMI)", y = "cells") + theme_bw(base_size = 10)
ggsave(file.path(plot.dir, glue::glue("{pfx}-bimodal-histogram.pdf")), p_hist, width = 7, height = 4.5)
readr::write_csv(p_hist$data, file.path(plot.dir, glue::glue("{pfx}-bimodal-histogram-sourcedata.csv")))

# 2) specificity by cell type
thr <- 1
fp <- md %>% group_by(Level_4, is_neopl) %>%
  summarise(n = n(), pct_tg_pos = 100 * mean(tg_sum >= thr), median_tg = median(tg_sum), .groups = "drop") %>%
  arrange(is_neopl, desc(pct_tg_pos))
write_csv(fp, file.path(plot.dir, glue::glue("{pfx}-positivity-by-celltype.csv")))
fp_rate <- md %>% filter(!is_neopl) %>% summarise(fp = 100 * mean(tg_sum >= thr)) %>% pull(fp)
p_spec <- ggplot(md, aes(reorder(Level_4, tg_log, median), tg_log, fill = is_neopl)) +
  geom_violin(scale = "width", alpha = 0.85) +
  scale_fill_manual(values = cols_neopl, guide = "none") + coord_flip() +
  labs(title = "Transgene by cell type", subtitle = NULL, x = NULL, y = "log1p(transgene UMI)") +
  theme_bw(base_size = 9)
ggsave(file.path(plot.dir, glue::glue("{pfx}-by-celltype.pdf")), p_spec, width = 7, height = 7)
readr::write_csv(p_spec$data, file.path(plot.dir, glue::glue("{pfx}-by-celltype-sourcedata.csv")))

# 3) ambient correction (decontX)
cli_alert_info("Running decontX (ambient correction)...")
dx <- tryCatch(celda::decontX(x = counts, z = as.character(obj$Level_4)),
               error = function(e) { cli_alert_warning("decontX failed: {e$message}"); NULL })
if (!is.null(dx)) {
  md$tg_sum_decont <- Matrix::colSums(dx$decontXcounts[tg_feats, , drop = FALSE])
  p_amb <- md %>% select(is_neopl, tg_sum, tg_sum_decont) %>%
    pivot_longer(c(tg_sum, tg_sum_decont), names_to = "type", values_to = "tg") %>%
    mutate(type = recode(type, tg_sum = "raw", tg_sum_decont = "decontX")) %>%
    ggplot(aes(interaction(type, is_neopl), log1p(tg), fill = is_neopl)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = cols_neopl, guide = "none") +
    labs(title = "Transgene before/after ambient removal (decontX)",
         subtitle = "non-neoplastic transgene is largely ambient; neoplastic retained",
         x = NULL, y = "log1p(transgene UMI)") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(plot.dir, glue::glue("{pfx}-decontX-ambient.pdf")), p_amb, width = 7, height = 5)
  readr::write_csv(p_amb$data, file.path(plot.dir, glue::glue("{pfx}-decontX-ambient-sourcedata.csv")))
}

# 4) head-to-head: transgene vs inferCNV CNA score
metrics <- c(glue::glue("transgene features: {paste(tg_feats, collapse=', ')}"),
             glue::glue("non-neoplastic transgene+ false-positive rate: {round(fp_rate,2)}%"),
             glue::glue("AUC transgene (all cells): {round(auc(md$tg_sum, md$is_neopl),3)}"))
if (file.exists(icnv_path)) {
  cli_alert_info("Loading cached inferCNV object for CNA head-to-head...")
  icnv <- readRDS(icnv_path)
  cna <- colMeans(abs(icnv@expr.data - 1), na.rm = TRUE)
  md$cna_score <- cna[md$barcode]
  sub <- md %>% filter(!is.na(cna_score))
  auc_tg <- auc(sub$tg_sum, sub$is_neopl); auc_cna <- auc(sub$cna_score, sub$is_neopl)
  metrics <- c(metrics, glue::glue("cells with CNA score: {nrow(sub)}"),
               glue::glue("AUC transgene (matched): {round(auc_tg,3)}"),
               glue::glue("AUC inferCNV CNA score : {round(auc_cna,3)}"))
  p_sc <- ggplot(sub, aes(cna_score, tg_log, colour = is_neopl)) +
    geom_point(size = 0.3, alpha = 0.3) +
    scale_colour_manual(values = cols_neopl, labels = lab_neopl, name = NULL) +
    labs(title = glue::glue("Transgene vs inferCNV CNA score (AUC {round(auc_tg,2)} vs {round(auc_cna,2)})"),
         subtitle = "transgene separates neoplastic/reference; CNA score does not",
         x = "inferCNV CNA score (mean |dev|)", y = "log1p(transgene UMI)") +
    theme_bw(base_size = 9) + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  readr::write_csv(
    sub[, c("barcode", "Level_4", "is_neopl", "cna_score", "tg_log")],
    file.path(plot.dir, glue::glue("{pfx}-vs-CNA-scatter-sourcedata.csv"))
  )
  p_sc <- ggrastr::rasterise(p_sc, dpi = 300)
  ggsave(file.path(plot.dir, glue::glue("{pfx}-vs-CNA-scatter.pdf")), p_sc, width = 7.5, height = 5.5)
} else {
  cli_alert_warning("inferCNV object not found at {icnv_path} - skipped CNA scatter")
}
writeLines(metrics, file.path(plot.dir, glue::glue("{pfx}-metrics.txt")))
cli_alert_success("Reviewers-only transgene-vs-CNA panels -> {plot.dir}/{pfx}-*")
