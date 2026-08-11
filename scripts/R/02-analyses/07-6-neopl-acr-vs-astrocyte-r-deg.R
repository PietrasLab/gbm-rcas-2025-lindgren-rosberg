# 07-6 - Direct DEG: Neopl-ACR (malignant ACR) vs Astrocyte R (reactive astrocyte)
#
# Question: explicitly show the genes that distinguish a genuinely
# reactive astrocyte (Astrocyte R) from a tumour cell that has adopted the same
# reactive program (Neopl-ACR). Fits alongside Fig 5B.
#
# Method: single-cell Wilcoxon FindMarkers, TWO groups only (Neopl-ACR vs
# Astrocyte R), BOTH directions. This is the same DE framework used for every
# other single-cell volcano in the manuscript (SeuratExtend::VolcanoPlot in
# Figure 3), so this panel is method- and style-consistent with the rest of the
# paper (one DE method across the manuscript).
#
# Outputs:
#   results/07-6-neopl-acr-vs-astrocyte-r-deg/07-6-FindMarkers-NeoplACR-vs-AstrocyteR.csv
#   manuscript-figures/figure-5/Figure-5-Volcano-NeoplACR-vs-AstrocyteR.pdf
#   manuscript-figures/figure-5/Figure-5-Volcano-NeoplACR-vs-AstrocyteR-sourcedata.csv

# ---- Prepare Env ----
require(conflicted)
require(tidyverse)
require(cli)
require(glue)
suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
require(ggrepel)
require(SeuratExtend)
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, .quiet = TRUE)

set.seed(169999)

## ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

## ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

out.dir  <- "./results/07-6-neopl-acr-vs-astrocyte-r-deg"
fig.dir  <- "./manuscript-figures/figure-5"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# Step 1 - Load object, subset to the two populations
# ============================================================
cli::cli_alert_info("Loading Seurat object ...")
seurat.object <- readRDS("./data/processed/seurat/seurat_flex_filtered_v1.0.rds")

annot        <- "Level_4"
groups       <- c("Astrocyte R", "Neopl-ACR")
seurat.sub   <- subset(seurat.object, subset = Level_4 %in% groups)
seurat.sub   <- Seurat::NormalizeData(seurat.sub, verbose = FALSE)
Idents(seurat.sub) <- annot

cli::cli_alert_info("Astrocyte R: {sum(Idents(seurat.sub) == 'Astrocyte R')} cells")
cli::cli_alert_info("Neopl-ACR:   {sum(Idents(seurat.sub) == 'Neopl-ACR')} cells")
rm(seurat.object); gc()


# ============================================================
# Step 2 - FindMarkers (Neopl-ACR vs Astrocyte R), full table
# ============================================================
# logfc.threshold = 0 + only.pos = FALSE -> comprehensive, both-direction table
# suitable for a supplemental DEG list. Positive avg_log2FC = up in Neopl-ACR.
cli::cli_alert_info("FindMarkers: Neopl-ACR vs Astrocyte R ...")
deg_df <- Seurat::FindMarkers(
  seurat.sub,
  ident.1         = "Neopl-ACR",
  ident.2         = "Astrocyte R",
  only.pos        = FALSE,
  min.pct         = 0.1,
  logfc.threshold = 0
) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(
    comparison = "Neopl-ACR_vs_Astrocyte-R",
    direction  = case_when(
      p_val_adj < 0.05 & avg_log2FC >  1 ~ "Neopl-ACR",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Astrocyte R",
      TRUE                               ~ "Shared/NS"
    )
  )

cli::cli_alert_success(
  "DE (padj<0.05, |avg_log2FC|>1): Neopl-ACR {sum(deg_df$direction=='Neopl-ACR')} | Astrocyte R {sum(deg_df$direction=='Astrocyte R')} | Shared/NS {sum(deg_df$direction=='Shared/NS')}"
)

readr::write_csv(deg_df, file.path(out.dir, "07-6-FindMarkers-NeoplACR-vs-AstrocyteR.csv"))


# ============================================================
# Step 3 - Volcano (styled like the Figure-3 astrocyte volcanoes)
# ============================================================
Y_CUT     <- -log10(0.05)                 # p_adj < 0.05
X_MOD     <- 0.5; X_STRONG <- 1           # logFC guide lines / colour tiers
tier_cols <- c(ns = "grey72", modest = "#E8807B", strong = "#C0161B")

vp <- SeuratExtend::VolcanoPlot(
  seurat.sub,
  ident.1     = "Neopl-ACR",
  ident.2     = "Astrocyte R",
  top.n       = 20,
  x.threshold = X_MOD,
  y.threshold = Y_CUT
)

d <- vp$data
d$tier <- dplyr::case_when(
  d$significant & abs(d$logFC) > X_STRONG ~ "strong",
  d$significant                           ~ "modest",
  TRUE                                    ~ "ns"
)
d$tier <- factor(d$tier, levels = c("ns", "modest", "strong"))
lab <- d[!is.na(d$label) & d$label != "", ]
xr  <- range(d$logFC, na.rm = TRUE)

p_volcano <- ggplot(d, aes(logFC, p)) +
  geom_point(aes(colour = tier), size = 0.7) +
  geom_vline(xintercept = c(-X_MOD, X_MOD),       linetype = "dotted", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-X_STRONG, X_STRONG), linetype = "dashed", colour = "grey35", linewidth = 0.3) +
  geom_hline(yintercept = Y_CUT,                  linetype = "dashed", colour = "grey35", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = lab, aes(label = label, colour = tier),
                           size = 2.6, max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = tier_cols, guide = "none") +
  coord_cartesian(xlim = xr) +
  labs(
    title = "Neopl-ACR vs Astrocyte R",
    x     = expression(log[2]~"FC  (Neopl-ACR / Astrocyte R)"),
    y     = vp$labels$y
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(fig.dir, "Figure-5-Volcano-NeoplACR-vs-AstrocyteR.pdf"),
       p_volcano, width = 7, height = 6)
readr::write_csv(d, file.path(fig.dir, "Figure-5-Volcano-NeoplACR-vs-AstrocyteR-sourcedata.csv"))

cli::cli_alert_success("07-6 done. Table + volcano written.")
