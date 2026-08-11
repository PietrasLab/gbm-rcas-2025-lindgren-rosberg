# 09-1b — Nomura correlation CSV with COALESCED states (interim)
#
# The local datasets are re-scored by 09-1 with gene-UNION coalesced Nomura states (NEU=MP9+MP14,
# Stress=MP10+MP15; drop unassigned). Nomura itself is scored HPC-side and only the per-MP scores
# were pulled (C18c-Nomura-scores.rds). Until the HPC scoring is re-run with the coalesced states, we
# coalesce here by SCORE-AVERAGE of the member MP module scores (a proxy for the gene-union score;
# the AC/MES/OPC/NPC single-MP states are unchanged — only NEU & Stress are averaged). Clearly an
# INTERIM: replace 31-corr-Nomura-neoplastic.csv with the HPC gene-union version when available.

suppressMessages({library(dplyr); library(readr); library(tibble)})
R   <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
IN  <- file.path(R, "data/external/C18c-Nomura-scores.rds")
OUT <- file.path(R, "results/31-external-validation/31-corr-Nomura-neoplastic.csv")

s <- readRDS(IN); sm <- s$neo_score_mat; meta <- as.data.frame(s$neo_meta)
newcols <- list(); newmeta <- list()
keep_as_is <- meta$label[meta$source %in% c("Neftel 2019","Richards 2021")]
for (l in keep_as_is) { newcols[[l]] <- sm[, l]; newmeta[[l]] <- meta[meta$label==l, c("label","source","compartment")] }
# ours: merge Neopl-CC + Neopl-Bulk -> Neopl-CC/Bulk; keep the rest
ours <- meta$label[meta$source=="Ours"]
for (l in setdiff(ours, c("Neopl-CC (ours)","Neopl-Bulk (ours)"))) { newcols[[l]] <- sm[, l]; newmeta[[l]] <- tibble(label=l, source="Ours", compartment="neoplastic") }
if (all(c("Neopl-CC (ours)","Neopl-Bulk (ours)") %in% colnames(sm))) {
  newcols[["Neopl-CC/Bulk (ours)"]] <- rowMeans(sm[, c("Neopl-CC (ours)","Neopl-Bulk (ours)")])
  newmeta[["Neopl-CC/Bulk (ours)"]] <- tibble(label="Neopl-CC/Bulk (ours)", source="Ours", compartment="neoplastic") }
# Nomura: rename single-MP states; coalesce NEU (MP9+MP14) & Stress (MP10+MP15) by score-average
rename_map <- c(
  "Nomura MP2_OPC (OPC-like)"="Nomura OPC-like", "Nomura MP4_AC (AC-like)"="Nomura AC-like",
  "Nomura MP5_Hypoxia (Hypoxia)"="Nomura Hypoxia", "Nomura MP6_MES (MES-like)"="Nomura MES-like",
  "Nomura MP7_NPC (NPC-like)"="Nomura NPC-like", "Nomura MP8_GPC (GPC-like)"="Nomura GPC-like",
  "Nomura MP13_Cilia (Cilia-like)"="Nomura Cilia-like", "Nomura MP3_CC (Cell cycle)"="Nomura Cell cycle")
for (old in names(rename_map)) { newcols[[rename_map[[old]]]] <- sm[, old]
  newmeta[[rename_map[[old]]]] <- tibble(label=rename_map[[old]], source="Nomura 2025", compartment="neoplastic") }
newcols[["Nomura NEU-like"]] <- rowMeans(sm[, c("Nomura MP9_ExN (NEU-like)","Nomura MP14_NRGN (NEU-like)")])
newcols[["Nomura Stress"]]   <- rowMeans(sm[, c("Nomura MP10_Stress1 (Stress)","Nomura MP15_Stress2 (Stress)")])
newmeta[["Nomura NEU-like"]] <- tibble(label="Nomura NEU-like", source="Nomura 2025", compartment="neoplastic")
newmeta[["Nomura Stress"]]   <- tibble(label="Nomura Stress",   source="Nomura 2025", compartment="neoplastic")

M <- do.call(cbind, newcols); mt <- bind_rows(newmeta)
cm <- cor(M, method="spearman", use="pairwise.complete.obs")
df <- as.data.frame(round(cm,4)) |> rownames_to_column("signature") |>
  left_join(mt[,c("label","source","compartment")], by=c("signature"="label"))
write_csv(df, OUT); cat("wrote (interim, coalesced) 31-corr-Nomura-neoplastic.csv |", nrow(df), "sigs\n")
