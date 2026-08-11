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
  library(tibble)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
})
conflicts_prefer(dplyr::filter, dplyr::select, .quiet = TRUE)
source("./scripts/R/00-0-source-functions.R")   # canonical .color_pal
future::plan("sequential")   # cf. 15-1
options(future.globals.maxSize = 8 * 1024^3)  # FindMarkers' globals check runs even under sequential; raise the ~3 GiB object past the 500 MiB default

# 15-2 — Supp Fig 6 (Figure 6 supplement): metabolic / lipid programs vs L-R signalling
# (revision). Lets metabolic and signalling HALLMARK
# sets compete on one ranked list, and scores a curated metabolic/lipid panel across
# neoplastic states / conditions.
#   A) fgsea of the Neopl-ACR program (Neopl-ACR vs other neoplastic) against
#      HALLMARK; metabolic/lipid sets highlighted.
#   B) per-cell metabolic/lipid module scores per neoplastic state + by condition.

obj_path <- "data/processed/seurat/seurat_flex_filtered_v1.0.rds"
plot.dir <- "./manuscript-figures/figure-6"
pfx      <- "Figure-S6-metabolic"
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
metab_pat <- "METAB|LIPID|FATTY|CHOLESTEROL|STEROL|GLYCOLY|OXIDATIVE_PHOS|ADIPO|BILE_ACID|SPHINGO|GLYCEROPHOSPHO"

# ---- gene sets ----------------------------------------------------------------
H <- msigdbr(species = "Mus musculus", category = "H"); H_list <- split(H$gene_symbol, H$gs_name)
C2 <- msigdbr(species = "Mus musculus", category = "C2")
subcol <- if ("gs_subcat" %in% names(C2)) "gs_subcat" else "gs_subcollection"
lipid_c2 <- C2 |> filter(grepl("CP:KEGG|CP:REACTOME", .data[[subcol]]),
  grepl("LIPID|FATTY_ACID|CHOLESTEROL|STEROL|SPHINGOLIPID|GLYCEROPHOSPHOLIPID", gs_name))
panel <- list(
  Lipid_FattyAcid  = H_list[["HALLMARK_FATTY_ACID_METABOLISM"]],
  Cholesterol      = H_list[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]],
  Adipogenesis     = H_list[["HALLMARK_ADIPOGENESIS"]],
  OxPhos           = H_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
  Glycolysis       = H_list[["HALLMARK_GLYCOLYSIS"]],
  Lipid_KEGG_React = unique(lipid_c2$gene_symbol))

# ---- object: neoplastic cells -------------------------------------------------
cli_alert_info("Loading object + subsetting neoplastic...")
obj <- readRDS(obj_path)
neopl_levels <- grep("^Neopl", levels(factor(obj$Level_4)), value = TRUE)
neo <- subset(obj, subset = Level_4 %in% neopl_levels); rm(obj); gc()
neo <- NormalizeData(neo, verbose = FALSE)
neo$sample_type <- factor(neo$sample_type, levels = c("Healthy","Primary","Recurrent"))

# ---- (B) module scores across states / conditions -----------------------------
panel <- lapply(panel, function(g) intersect(g, rownames(neo)))
neo <- AddModuleScore(neo, features = panel, name = "MS_", seed = 1)
names(neo@meta.data)[match(paste0("MS_", seq_along(panel)), names(neo@meta.data))] <- names(panel)
md <- neo@meta.data |> select(Level_4, sample_type, all_of(names(panel))) |>
  mutate(Level_4 = factor(Level_4, levels = neopl_levels))

by_state <- md |> group_by(Level_4) |> summarise(across(all_of(names(panel)), mean), .groups = "drop")
write_csv(by_state, file.path(plot.dir, glue("{pfx}-score-by-state.csv")))
write_csv(md |> group_by(Level_4, sample_type) |> summarise(across(all_of(names(panel)), mean), .groups="drop"),
          file.path(plot.dir, glue("{pfx}-score-by-state-condition.csv")))

p_state <- by_state |> pivot_longer(-Level_4, names_to = "program", values_to = "score") |>
  ggplot(aes(program, Level_4, fill = score)) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(title = "Metabolic/lipid module scores per neoplastic state", x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(plot.dir, glue("{pfx}-score-by-state-heatmap.pdf")), p_state, width = 8, height = 5)
write_csv(p_state$data, file.path(plot.dir, glue("{pfx}-score-by-state-heatmap-sourcedata.csv")))

# (lipid-score-by-condition violin dropped — differences across conditions are
#  minimal; the per-state heatmap above carries the metabolic/lipid read-out.)

# ---- (A) fgsea: Neopl-ACR program vs HALLMARK ---------------------------------
acr_like <- intersect("Neopl-ACR", neopl_levels)
Idents(neo) <- ifelse(neo$Level_4 %in% acr_like, "Neopl_ACR", "other_neoplastic")
dm <- FindMarkers(neo, ident.1 = "Neopl_ACR", ident.2 = "other_neoplastic",
                  logfc.threshold = 0, min.pct = 0.05, verbose = FALSE)
ranks <- dm$avg_log2FC; names(ranks) <- rownames(dm)
ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
set.seed(1)
fg <- fgsea(pathways = H_list, stats = ranks, eps = 0, minSize = 10, maxSize = 500) |>
  as_tibble() |> arrange(desc(NES)) |>
  mutate(metabolic = grepl(metab_pat, pathway),
         leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))
write_csv(fg, file.path(plot.dir, glue("{pfx}-fgsea-NeoplACR-HALLMARK.csv")))

top <- fg |> filter(padj < 0.05) |> slice_max(abs(NES), n = 20) |>
  mutate(pathway = sub("^HALLMARK_", "", pathway), pathway = factor(pathway, levels = pathway[order(NES)]))
p_fg <- ggplot(top, aes(NES, pathway, fill = metabolic)) + geom_col() +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"),
                    labels = c("other", "metabolic/lipid"), name = NULL) +
  labs(title = "HALLMARK enrichment in the Neopl-ACR program",
       subtitle = "Positive NES = up in Neopl-ACR vs other neoplastic; metabolic/lipid sets highlighted",
       x = "NES", y = NULL) + theme_bw(base_size = 9)
ggsave(file.path(plot.dir, glue("{pfx}-fgsea-barplot.pdf")), p_fg, width = 9, height = 7)
write_csv(p_fg$data, file.path(plot.dir, glue("{pfx}-fgsea-barplot-sourcedata.csv")))

cli_alert_success("Metabolic/lipid panels -> {plot.dir}/{pfx}-* (metabolic sets among {sum(fg$padj<0.05,na.rm=TRUE)} sig HALLMARK)")
