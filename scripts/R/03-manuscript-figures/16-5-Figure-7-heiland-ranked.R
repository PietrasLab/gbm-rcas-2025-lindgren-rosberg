# 16-5 — Figure 7E: Heiland human-GBM ranked correlations (per-section, multi-target)
# ============================================================================
# For each of our neoplastic programs (targets), rank the external reference programs
# (Neftel / Richards / Nomura) by their mean per-section Spearman correlation across the Heiland
# human GBM tumor sections; err bar = min/max, points = individual sections.
# Reads the saved per-section correlation matrices
# (heiland-per-section-cor-list.rds) from results/08-visium-hd-derived (env VISIUM_HD_DERIVED_DIR);
# targets/references are derived from the matrix column names (Ours neoplastic vs Neftel/Richards/
# Nomura), so the heavy signature/geneset panel does not need to be re-built here.
# Per-SECTION Spearman; section = replicate unit; no per-spot significance (autocorrelation-
# controlled test = toroidal null, run upstream).
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(forcats); library(ggplot2); library(patchwork); library(readr); library(cli)
})
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
DERIVED   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))

cl <- readRDS(file.path(DERIVED, "heiland-per-section-cor-list.rds"))
allcols <- colnames(cl[[1]])
# targets = our neoplastic programs (both ACR variants); refs = external neoplastic programs.
ours_neo <- grep("^Neopl-.*\\(ours", allcols, value = TRUE)
refs     <- grep("^(Neftel|Richards|Nomura)", allcols, value = TRUE)
src_of   <- setNames(ifelse(grepl("^Neftel", refs), "Neftel 2019",
                     ifelse(grepl("^Richards", refs), "Richards 2021", "Nomura 2025")), refs)
SRC_PAL  <- c("Neftel 2019"="#1F78B4","Richards 2021"="#33A02C","Nomura 2025"="#FF7F00")
cli::cli_alert_info("{length(cl)} Heiland sections; {length(ours_neo)} targets x {length(refs)} refs")

# per-section correlation target x ref, from the saved cor matrices
cor_ps <- imap_dfr(cl, function(mat, sid) {
  tg <- intersect(ours_neo, colnames(mat)); rf <- intersect(refs, colnames(mat))
  expand.grid(target = tg, program = rf, stringsAsFactors = FALSE) |>
    mutate(section = sid, r = mapply(function(t, p) mat[t, p], target, program))
})
cor_ps$source <- src_of[cor_ps$program]
summ <- cor_ps |> group_by(target, program, source) |>
  summarise(mean_r = mean(r, na.rm=TRUE), min_r = min(r, na.rm=TRUE), max_r = max(r, na.rm=TRUE), n = n(), .groups = "drop")
write_csv(cor_ps, file.path(plot_dir, "Figure-7-heiland-ranked-persection-sourcedata.csv"))
write_csv(summ,   file.path(plot_dir, "Figure-7-heiland-ranked-summary-sourcedata.csv"))

plot_target <- function(tg) {
  s <- summ[summ$target == tg, ]; ps <- cor_ps[cor_ps$target == tg, ]
  ord <- s |> arrange(mean_r) |> pull(program); s$program <- factor(s$program, ord); ps$program <- factor(ps$program, ord)
  ggplot(s, aes(mean_r, program)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_col(aes(fill = source), width = 0.7, alpha = 0.85) +
    geom_errorbarh(aes(xmin = min_r, xmax = max_r), height = 0.25, linewidth = 0.3, colour = "grey30") +
    geom_point(data = ps, aes(r, program), size = 0.5, colour = "grey20", alpha = 0.5) +
    scale_fill_manual(values = SRC_PAL, name = "source") +
    labs(title = tg, x = "Spearman r (per-section)", y = NULL) + theme_bw(base_size = 8) + theme(legend.position = "bottom")
}
p <- wrap_plots(lapply(ours_neo, plot_target), ncol = 3, guides = "collect") +
  plot_annotation(title = sprintf("Fig 7E (Heiland, %d tumor sections): reference programs ranked by mean correlation with our neoplastic programs", length(cl)),
    subtitle = "bar=mean; err=min/max; points=per-section. Per-SECTION Spearman; section=replicate; no per-spot significance (autocorrelation-controlled test = toroidal null, run upstream).",
    theme = theme(legend.position = "bottom"))
ggsave(file.path(plot_dir, "Figure-7-fig7E-heiland-ranked.pdf"), p, width = 15, height = 9)
cli::cli_alert_success("16-5 Heiland ranked correlations -> {plot_dir}")
