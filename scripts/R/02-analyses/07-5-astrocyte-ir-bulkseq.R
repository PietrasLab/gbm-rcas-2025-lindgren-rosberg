# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(glue)
  library(DESeq2)
  library(fgsea)
  library(homologene)
  library(dplyr)
  library(readr)
  library(tibble)
})
conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, .quiet = TRUE)

# 07-5 — Bulk RNAseq validation: irradiated vs non-irradiated human astrocytes
# (revision)
#
# External data: PietrasLab/Jeannot-GBM-astrocytes-Flunarizine (data/rnaseq).
#   Clone it and symlink data/external/jeannot-gbm-astrocytes-flunarizine -> the clone.
#   Human primary astrocytes, 4 conditions x 4 replicates; we use the DMSO arm to
#   isolate the pure IR effect (IR_DMSO vs nonIR_DMSO).
#
# Scope: GAS6 is the focus gene; the Astrocyte R GSEA goes to MAIN Figure 3,
# supporting a marginal-but-real increase in the ACR program upon IR.
#
# GSEA gene set = the FULL Astrocyte R cell-state program: all significant
# Astrocyte R cluster markers (07-1 Level_4, padj<0.05, only.pos), human orthologs
# (homologene). A ranked GSEA needs a broad set; the 31-gene refined signature
# (07-4) is for module scoring / survival, not for this enrichment test.
#
# Outputs (results/07-5-astrocyte-ir-bulkseq/) — consumed by 12-2-Figure-3-bulkseq.R:
#   - DESeq2-IR-vs-nonIR-DMSO.csv          full contrast (Symbol, log2FC, stat, padj)
#   - ranked-stats-IR.csv                  Symbol, stat (for plotEnrichment)
#   - astrocyte_R_program_hs.csv           gene set used for the GSEA
#   - GSEA-AstrocyteR-program-in-IR.csv    fgsea result (ES/NES/padj/leadingEdge)
#   - GAS6-normcounts-DMSO.csv             per-sample normalised counts (boxplot)

dds_path     <- "data/external/jeannot-gbm-astrocytes-flunarizine/data/rnaseq/deseq2.filtered.rds"
markers_path <- "results/07-1-deg-findmarkers/07-1-FindAllMarkers-Level_4.csv"
out_dir      <- "results/07-5-astrocyte-ir-bulkseq"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GENE_FOCUS <- "GAS6"   # Gas6 only in the figure

# ---- Load + IR vs nonIR contrast (DMSO only) ----------------------------------
cli_alert_info("Loading DESeq2 object: {dds_path}")
dds_full <- readRDS(dds_path)

dmso <- c("nonIR_DMSO", "IR_DMSO")
dds  <- dds_full[, colData(dds_full)$Condition %in% dmso]
colData(dds)$Condition <- factor(colData(dds)$Condition, levels = dmso)
design(dds) <- ~ Condition
dds <- DESeq(dds, quiet = TRUE)

fdata <- as.data.frame(rowData(dds_full))
res <- results(dds, contrast = c("Condition", "IR_DMSO", "nonIR_DMSO")) |>
  as.data.frame() |>
  rownames_to_column("ensembl_gene_id") |>
  mutate(Symbol = fdata$external_gene_name_ensembl[match(ensembl_gene_id, fdata$ensembl_gene_id)]) |>
  filter(!is.na(Symbol), Symbol != "") |>
  arrange(padj)

write_csv(res, file.path(out_dir, "DESeq2-IR-vs-nonIR-DMSO.csv"))
g6 <- res |> filter(Symbol == "GAS6")
cli_alert_info("DESeq2 IR vs nonIR (DMSO): {nrow(res)} genes; GAS6 LFC {round(g6$log2FoldChange,2)} padj {signif(g6$padj,3)}")

# ---- Astrocyte R program gene set: all sig cluster markers -> human --------------
mk <- read_csv(markers_path, show_col_types = FALSE)
astroR_mm <- mk |> filter(cluster == "Astrocyte R", p_val_adj < 0.05, avg_log2FC > 0) |> pull(gene)
astroR_hs <- homologene::mouse2human(astroR_mm) |>
  filter(!is.na(humanGene), humanGene != "") |> distinct(humanGene) |> pull(humanGene)
write_csv(tibble(gene = astroR_hs), file.path(out_dir, "astrocyte_R_program_hs.csv"))

# ---- GSEA: Astrocyte R program among IR-up genes --------------------------------
ranked <- res |> filter(!is.na(stat)) |> arrange(desc(stat))
ranked_vec <- setNames(ranked$stat, ranked$Symbol)
ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]
write_csv(tibble(Symbol = names(ranked_vec), stat = as.numeric(ranked_vec)),
          file.path(out_dir, "ranked-stats-IR.csv"))

n_overlap <- sum(astroR_hs %in% names(ranked_vec))
cli_alert_info("Astrocyte R program: mm {length(astroR_mm)} -> hs {length(astroR_hs)}; in bulk universe {n_overlap}")

set.seed(42)
gsea <- fgsea(pathways = list(AstrocyteR_program = astroR_hs),
              stats = ranked_vec, minSize = 5, maxSize = 2000, nPermSimple = 10000)
gsea_out <- gsea |>
  mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))
write_csv(gsea_out, file.path(out_dir, "GSEA-AstrocyteR-program-in-IR.csv"))
cli_alert_info("GSEA Astrocyte R program: NES {round(gsea$NES,2)}, padj {signif(gsea$padj,3)}")

# ---- GAS6 normalised counts (DMSO samples) for the boxplot ---------------------
nc <- counts(dds_full, normalized = TRUE)
cd <- as.data.frame(colData(dds_full))
ensg <- fdata$ensembl_gene_id[match(GENE_FOCUS, fdata$external_gene_name_ensembl)]
dmso_idx <- cd$Condition %in% dmso

write_csv(
  tibble(
    sample    = colnames(nc)[dmso_idx],
    Condition = factor(cd$Condition[dmso_idx], levels = dmso),
    Rep       = cd$Rep[dmso_idx],
    gene      = GENE_FOCUS,
    norm_count = as.numeric(nc[ensg, dmso_idx])
  ),
  file.path(out_dir, "GAS6-normcounts-DMSO.csv")
)

cli_alert_success("07-5 done -> {out_dir}/ (DESeq2 + Astrocyte R program GSEA + GAS6)")
