# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(glue)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(fgsea)
})
conflicts_prefer(dplyr::filter, .quiet = TRUE)

# 12-2 — bulk-seq panels from the human-astrocyte IR bulk-seq (07-5).
# (revision; GAS6 + Astrocyte-R GSEA)
#
# Inputs:  results/07-5-astrocyte-ir-bulkseq/  (run 07-5-astrocyte-ir-bulkseq.R first)
# Outputs: Astrocyte-R GSEA (panel B) -> manuscript-figures/figure-3/ (stays Figure 3)
# GAS6 (panel A) -> manuscript-figures/figure-6/ (reassigned Figure 3 -> Figure 6 per DL)

in_dir   <- "results/07-5-astrocyte-ir-bulkseq"
plot.dir <- glue("./manuscript-figures/figure-3")   # Astrocyte-R GSEA panel stays in Figure 3
gas6.dir <- glue("./manuscript-figures/figure-6")   # GAS6 panel reassigned to Figure 6 (per DL)
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gas6.dir, recursive = TRUE, showWarnings = FALSE)

cond_levels <- c("nonIR_DMSO", "IR_DMSO")
cond_labels <- c("nonIR", "IR")
cond_cols   <- c(nonIR_DMSO = "#4A90D9", IR_DMSO = "#D9534F")

# ---- Panel A: GAS6 normalised counts, nonIR vs IR (DMSO) -----------------------
gas6 <- read_csv(file.path(in_dir, "GAS6-normcounts-DMSO.csv"), show_col_types = FALSE) |>
  mutate(Condition = factor(Condition, levels = cond_levels))
deg  <- read_csv(file.path(in_dir, "DESeq2-IR-vs-nonIR-DMSO.csv"), show_col_types = FALSE)
g    <- deg |> filter(Symbol == "GAS6")
sub  <- glue("LFC = {round(g$log2FoldChange,2)}  |  padj = {signif(g$padj,2)}")

pA <- ggplot(gas6, aes(Condition, norm_count, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.45) +
  geom_line(aes(group = Rep), color = "grey55", linewidth = 0.4) +
  geom_point(size = 2.6, shape = 21, color = "grey20") +
  scale_fill_manual(values = cond_cols) +
  scale_x_discrete(labels = cond_labels) +
  labs(title = "GAS6", subtitle = sub, x = NULL, y = "Normalised counts") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold.italic", size = 13),
        plot.subtitle = element_text(size = 8, color = "grey40"))

ggsave(file.path(gas6.dir, "Figure-6-bulkseq-GAS6.pdf"), pA, width = 3.2, height = 4.3)
write_csv(pA$data, file.path(gas6.dir, "Figure-6-bulkseq-GAS6-sourcedata.csv"))

# ---- Panel B: Astrocyte-R program GSEA in IR-up genes --------------------------
ranked <- read_csv(file.path(in_dir, "ranked-stats-IR.csv"), show_col_types = FALSE)
rv     <- setNames(ranked$stat, ranked$Symbol); rv <- rv[!duplicated(names(rv))]
sig    <- read_csv(file.path(in_dir, "astrocyte_R_program_hs.csv"), show_col_types = FALSE)$gene
gsea   <- read_csv(file.path(in_dir, "GSEA-AstrocyteR-program-in-IR.csv"), show_col_types = FALSE)

pB <- fgsea::plotEnrichment(sig, rv) +
  labs(title = "Astrocyte R markers enriched in IR-upregulated genes",
       subtitle = glue("NES = {round(gsea$NES,2)}  |  padj = {signif(gsea$padj,2)}  |  n = {gsea$size} genes"),
       x = "Rank in IR vs nonIR (DESeq2 stat)", y = "Enrichment score") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, color = "grey40"))

ggsave(file.path(plot.dir, "Figure-3-bulkseq-GSEA-AstrocyteR.pdf"), pB, width = 5.4, height = 3.9)
# plotEnrichment $data is the running-enrichment-score curve (x = rank, y = ES)
write_csv(pB$data, file.path(plot.dir, "Figure-3-bulkseq-GSEA-AstrocyteR-sourcedata.csv"))

cat("Saved bulk-seq panels: GAS6 ->", gas6.dir, "| Astrocyte-R GSEA ->", plot.dir, "\n")
