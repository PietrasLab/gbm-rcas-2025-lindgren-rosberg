# ---- Prepare Env ----
set.seed(42)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(glue)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
})
conflicts_prefer(dplyr::filter, dplyr::select, .quiet = TRUE)
source("./scripts/R/00-0-source-functions.R")   # canonical .color_pal

# 10-5 — Supp Fig 1 (Figure 1 supplement): the neoplastic (RCAS) signature is a
# bona fide tumour program (revision).
# Parts A+B (Visium-HD cross-technique part C is a heavier follow-up):
#   A) HALLMARK fgsea on the endogenous neoplastic-marker ranking (reporter genes
#      excluded) -> tumourigenesis/proliferation/RTK programs.
#   B) Flex: the endogenous signature score alone separates neoplastic vs non (AUC),
#      a transgene-independent neoplastic call.

markers_csv <- "results/03-qc/03-6-rcas-findallmarkers.csv"
obj_path    <- "data/processed/seurat/seurat_flex_filtered_v1.0.rds"
plot.dir    <- "./manuscript-figures/figure-1"
pfx         <- "Figure-S1-neoplastic"
RCAS_N <- 50L

cols_neopl <- c(neoplastic = unname(.color_pal[["Level_1"]][["Neoplastic"]]),
                `non-neoplastic` = unname(.color_pal[["Level_1"]][["Non-Neoplastic"]]))
auc <- function(score, pos) {
  ok <- !is.na(score) & !is.na(pos); score <- score[ok]; pos <- as.logical(pos[ok])
  r <- rank(score); n1 <- as.numeric(sum(pos)); n0 <- as.numeric(sum(!pos))
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- endogenous neoplastic signature (reporter genes excluded) -----------------
mk <- read_csv(markers_csv, show_col_types = FALSE)
mk <- mk[mk$cluster == "rcas_pos", ]
mk_endo <- mk[!grepl("^Tg-", mk$gene), ]
rcas_endo_genes <- head(mk_endo$gene, RCAS_N)
cli_alert_info("rcas_pos markers {nrow(mk)} | endogenous {nrow(mk_endo)} | top-{RCAS_N} signature")

# ---- (A) HALLMARK fgsea: signature is a tumour program -------------------------
H_list <- split(msigdbr(species = "Mus musculus", category = "H")$gene_symbol,
                 msigdbr(species = "Mus musculus", category = "H")$gs_name)
ranks <- mk_endo$avg_log2FC; names(ranks) <- mk_endo$gene
ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
fg <- fgsea(pathways = H_list, stats = ranks, eps = 0, minSize = 5, maxSize = 500, scoreType = "pos") |>
  arrange(padj) |> mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))
write_csv(fg, file.path(plot.dir, glue("{pfx}-HALLMARK-fgsea.csv")))

p_fg <- fg |> filter(padj < 0.1) |> arrange(desc(NES)) |> head(20) |>
  mutate(pathway = sub("^HALLMARK_", "", pathway), pathway = factor(pathway, levels = rev(pathway))) |>
  ggplot(aes(NES, pathway, fill = -log10(padj))) +
  geom_col(colour = "grey25", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("#FDD0A2", "#FDAE6B", "#FF4A00"), name = "-log10(padj)") +  # Fig 5B scheme, darker low end
  labs(title = "HALLMARK enrichment of the neoplastic (RCAS) signature",
       subtitle = "fgsea on rcas_pos up-marker ranking (reporter genes excluded)", x = "NES", y = NULL) +
  theme_bw(base_size = 9)
ggsave(file.path(plot.dir, glue("{pfx}-HALLMARK-barplot.pdf")), p_fg, width = 9, height = 6)
write_csv(p_fg$data, file.path(plot.dir, glue("{pfx}-HALLMARK-barplot-sourcedata.csv")))
cli_alert_info("HALLMARK fgsea: {sum(fg$padj<0.05)} pathways padj<0.05")

# ---- (B) Flex: endogenous signature separates neoplastic (no transgene) --------
cli_alert_info("Loading Flex object + scoring endogenous signature...")
flex <- readRDS(obj_path); DefaultAssay(flex) <- "RNA"
flex <- tryCatch(SeuratObject::JoinLayers(flex), error = function(e) flex)
flex <- NormalizeData(flex, verbose = FALSE)
flex$is_neopl <- grepl("^Neopl", as.character(flex$Level_4))
endo_avail <- intersect(rcas_endo_genes, rownames(flex))
flex <- AddModuleScore(flex, features = list(endo_avail), name = "rcas_endo_", assay = "RNA")
flex$rcas_endo_score <- flex$rcas_endo_1
flex_auc <- auc(flex$rcas_endo_score, flex$is_neopl)
cli_alert_info("Flex endogenous-signature AUC (neoplastic vs non): {round(flex_auc,3)} [{length(endo_avail)}/{length(rcas_endo_genes)} genes]")

flex_md <- flex@meta.data |> rownames_to_column("bc") |> select(bc, Level_4, is_neopl, rcas_endo_score)
write_csv(flex_md, file.path(plot.dir, glue("{pfx}-flex-scores.csv.gz")))

p_flex <- flex_md |> mutate(grp = ifelse(is_neopl, "neoplastic", "non-neoplastic")) |>
  ggplot(aes(reorder(Level_4, rcas_endo_score, median), rcas_endo_score, fill = grp)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.2) + coord_flip() +
  scale_fill_manual(values = cols_neopl, name = NULL) +
  labs(title = "Neoplastic (RCAS) signature score by cell type (Flex)",
       subtitle = glue("top-{RCAS_N} rcas_pos markers, reporter (Tg) genes excluded; transgene-independent; neoplastic-vs-non AUC = {round(flex_auc,3)}"),
       x = NULL, y = "Neoplastic RCAS signature score") + theme_bw(base_size = 9)
ggsave(file.path(plot.dir, glue("{pfx}-signature-by-celltype.pdf")), p_flex, width = 8, height = 6)
write_csv(p_flex$data, file.path(plot.dir, glue("{pfx}-signature-by-celltype-sourcedata.csv")))

cli_alert_success("Neoplastic-signature panels (A+B) -> {plot.dir}/{pfx}-*")
