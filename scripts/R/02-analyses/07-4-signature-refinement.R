# ---- Prepare Env ----
set.seed(169)

suppressPackageStartupMessages({
  library(conflicted)
  library(cli)
  library(glue)
  library(dplyr)
  library(readr)
  library(homologene)
})
conflicts_prefer(dplyr::filter, dplyr::slice_max, .quiet = TRUE)

# 07-4 — Signature refinement (revision)
#
# Recreates the Astrocyte R and Neopl-ACR cell-state signatures from the 07-1
# Level_4 markers and removes the genes shared between them (non-neoplastic ACR
# genes that otherwise leak into the Neopl-ACR signature).
#
# Procedure:
#   - source: results/07-1-deg-findmarkers/07-1-FindAllMarkers-Level_4.csv
#   - filter p_val_adj < 0.05, top-50 genes per cluster by avg_log2FC  ("original")
#   - refined = each signature minus the overlap with the other
#   - map mouse -> human orthologs (homologene)
#
# Also emits the HUMAN gene lists used for the GEPIA survival analysis (GBM):
# refined ACR and refined Neopl-ACR individually, and refined vs original Neopl-ACR.
#
# Outputs (results/07-4-signature-refinement/), each {mm, hs}:
#   - astrocyte_R_original / astrocyte_R_refined
#   - neopl_ACR_original   / neopl_ACR_refined
#   + C3-signature-overlap-supplementary-table.csv

TOP_N <- 50
PADJ  <- 0.05

markers_file <- "results/07-1-deg-findmarkers/07-1-FindAllMarkers-Level_4.csv"
out_dir      <- "results/07-4-signature-refinement"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mouse_to_human <- function(mouse_syms) {
  homologene::mouse2human(mouse_syms) |>
    filter(!is.na(humanGene), humanGene != "") |>
    distinct(humanGene) |>
    pull(humanGene) |>
    sort()
}
write_sig <- function(genes, stem) {
  write_csv(data.frame(gene = genes),                file.path(out_dir, glue("{stem}_mm.csv")))
  write_csv(data.frame(gene = mouse_to_human(genes)), file.path(out_dir, glue("{stem}_hs.csv")))
}

# ---- Original signatures: top-50 per cluster by avg_log2FC ----------------------
markers_df <- read_csv(markers_file, show_col_types = FALSE)

sig_top50 <- markers_df |>
  filter(p_val_adj < PADJ) |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = TOP_N) |>
  ungroup()

astR_genes   <- sig_top50 |> filter(cluster == "Astrocyte R") |> arrange(desc(avg_log2FC)) |> pull(gene)
neoACR_genes <- sig_top50 |> filter(cluster == "Neopl-ACR")   |> arrange(desc(avg_log2FC)) |> pull(gene)

# ---- Overlap removal -> refined ------------------------------------------------
overlap        <- intersect(astR_genes, neoACR_genes)
astR_refined   <- setdiff(astR_genes,   neoACR_genes)
neoACR_refined <- setdiff(neoACR_genes, astR_genes)

cli_alert_info("Astrocyte R: original {length(astR_genes)} -> refined {length(astR_refined)}")
cli_alert_info("Neopl-ACR:   original {length(neoACR_genes)} -> refined {length(neoACR_refined)}")
cli_alert_info("Overlap removed (n = {length(overlap)}): {paste(overlap, collapse = ', ')}")

# ---- Write signatures (mm + hs) ------------------------------------------------
write_sig(astR_genes,     "astrocyte_R_original")
write_sig(astR_refined,   "astrocyte_R_refined")
write_sig(neoACR_genes,   "neopl_ACR_original")
write_sig(neoACR_refined, "neopl_ACR_refined")

# ---- Supplementary overlap table ----------------------------------------------
all_genes <- union(astR_genes, neoACR_genes)
write_csv(
  data.frame(
    gene            = all_genes,
    in_AstR_top50   = all_genes %in% astR_genes,
    in_NeoACR_top50 = all_genes %in% neoACR_genes,
    in_overlap      = all_genes %in% overlap,
    retained_AstR   = all_genes %in% astR_refined,
    retained_NeoACR = all_genes %in% neoACR_refined
  ) |> arrange(desc(in_AstR_top50), desc(in_NeoACR_top50)),
  file.path(out_dir, "C3-signature-overlap-supplementary-table.csv")
)

cli_alert_success("Wrote original + refined ACR / Neopl-ACR (mm + hs) to {out_dir}/ — _hs lists are the GEPIA deliverable")
