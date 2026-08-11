
# Nomura et al. 2025, Nature Genetics — the revisited GBM meta-programs / cell states
# "An integrative single-cell atlas ..." Nat Genet (2025). doi:10.1038/s41588-025-02167-5
# https://www.nature.com/articles/s41588-025-02167-5
#
# These are the UPDATED (GLASS/Nomura) meta-programs that reviewers asked us to relate
# our RCAS neoplastic programs to (cross-model validation). They supersede the
# original Neftel_Cell_2019 6-state model already curated in
#   references/genesets/Richards_NatCancer_2021_GeneSets{,_Mouse}.rds .
#
# Source (local): ~/Resources/Genesets/Nomura-Tirosh/41588_2025_2167_MOESM2_ESM.xlsx
#   - TableS2 = the 15 meta-programs (MP_1..MP_15), 50 human symbols each. Header on row 5.
#   - TableS3 = 3 broad biological-process programs (BP ECM / Neuronal / Glial). Header row 4.
#
# 15 MP -> 9 cell-state assignment is from the paper's Extended Data Fig. 3 (3 MPs unassigned).
#
# Mirrors the pattern of 01-3-genesets-Richards-2021-GBM.R: build a named list of human
# symbols, save it, then translate to mouse orthologs with homologene (human2mouse).

# ---- Env ----
suppressMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(homologene)
  library(cli)
})

setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

xlsx <- "~/Resources/Genesets/Nomura-Tirosh/41588_2025_2167_MOESM2_ESM.xlsx"
out.dir <- "./references/genesets"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

PFX <- "Nomura_NatGenet_2025_"  # FirstAuthor_Journal_Year_ convention (matches Neftel_Cell_2019_*)

# ---- Read the 15 meta-programs (TableS2) ----
# Header is on sheet row 5 -> skip the 4 descriptor rows above it.
mp_wide <- suppressMessages(read_excel(xlsx, sheet = "TableS2", skip = 4))
cli_alert_info("TableS2: {nrow(mp_wide)} rows x {ncol(mp_wide)} MP columns")

mp_list <- mp_wide %>%
  pivot_longer(cols = everything(), names_to = "geneset", values_to = "gene") %>%
  filter(!is.na(gene) & gene != "") %>%
  group_by(geneset) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

# ---- Read the 3 broad BP programs (TableS3) ----
# Header is on sheet row 4 -> skip 3.
bp_wide <- suppressMessages(read_excel(xlsx, sheet = "TableS3", skip = 3))
cli_alert_info("TableS3: {nrow(bp_wide)} rows x {ncol(bp_wide)} BP columns")

bp_list <- bp_wide %>%
  pivot_longer(cols = everything(), names_to = "geneset", values_to = "gene") %>%
  filter(!is.na(gene) & gene != "") %>%
  # tidy the BP names ("BP ECM" -> "BP_ECM")
  mutate(geneset = gsub("\\s+", "_", trimws(geneset))) %>%
  group_by(geneset) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

# ---- Combine + apply the study prefix ----
nomura_human <- c(mp_list, bp_list)
names(nomura_human) <- paste0(PFX, names(nomura_human))
cli_alert_success("Assembled {length(nomura_human)} Nomura genesets (human): {sum(names(nomura_human) %in% names(nomura_human))} total")
print(vapply(nomura_human, length, integer(1)))

# ---- 15 MP -> 9 cell-state map (Extended Data Fig. 3) ----
# state = the coarse cell state; NA-state MPs (RP / MIC / LQ) are "Unassigned",
# and MP3 (CC) is the cell-cycle program (not a state). Keys are the RAW column names.
mp_state_map <- tibble::tribble(
  ~mp,             ~state,        ~note,
  "MP_1_RP",       "Unassigned",  "Ribosomal protein",
  "MP_2_OPC",      "OPC-like",    "",
  "MP_3_CC",       "Cell cycle",  "Proliferation (not a cell state)",
  "MP_4_AC",       "AC-like",     "",
  "MP_5_Hypoxia",  "Hypoxia",     "",
  "MP_6_MES",      "MES-like",    "",
  "MP_7_NPC",      "NPC-like",    "",
  "MP_8_GPC",      "GPC-like",    "",
  "MP_9_ExN",      "NEU-like",    "Excitatory-neuron-like",
  "MP_10_Stress1", "Stress",      "",
  "MP_11_MIC",     "Unassigned",  "Microglia contamination",
  "MP_12_LQ",      "Unassigned",  "Low quality",
  "MP_13_Cilia",   "Cilia-like",  "",
  "MP_14_NRGN",    "NEU-like",    "NRGN neuronal",
  "MP_15_Stress2", "Stress",      ""
)

# ---- Attach metadata as attributes so downstream code can group MP -> state ----
attr(nomura_human, "mp_state_map") <- mp_state_map
attr(nomura_human, "prefix")       <- PFX
attr(nomura_human, "source")       <- basename(xlsx)
attr(nomura_human, "citation")     <- "Nomura et al. 2025, Nat Genet, doi:10.1038/s41588-025-02167-5"

saveRDS(nomura_human, file.path(out.dir, "Nomura_NatGenet_2025_GeneSets.rds"))
cli_alert_success("Saved human genesets -> {file.path(out.dir, 'Nomura_NatGenet_2025_GeneSets.rds')}")

# ---- Translate to mouse orthologs (homologene::human2mouse) ----
cli_h1("Converting human genes to mouse orthologs (homologene)")
nomura_mouse <- map(nomura_human, ~ {
  hm <- homologene::human2mouse(.x)      # 9606 -> 10090
  unique(hm$mouseGene)
})
attributes(nomura_mouse) <- attributes(nomura_human)  # carry map/citation forward
attr(nomura_mouse, "species") <- "mouse (homologene human2mouse)"

saveRDS(nomura_mouse, file.path(out.dir, "Nomura_NatGenet_2025_GeneSets_Mouse.rds"))
cli_alert_success("Saved mouse genesets -> {file.path(out.dir, 'Nomura_NatGenet_2025_GeneSets_Mouse.rds')}")

# ---- QC: gene counts + mouse-mapping yield ----
qc <- tibble(
  geneset   = names(nomura_human),
  n_human   = vapply(nomura_human, length, integer(1)),
  n_mouse   = vapply(nomura_mouse, length, integer(1))
) %>%
  mutate(mouse_yield = round(n_mouse / n_human, 2)) %>%
  left_join(mp_state_map %>% mutate(geneset = paste0(PFX, mp)), by = "geneset")
print(as.data.frame(qc), row.names = FALSE)
readr::write_csv(qc, file.path(out.dir, "Nomura_NatGenet_2025_GeneSets-QC.csv"))
cli_alert_success("Wrote QC table -> {file.path(out.dir, 'Nomura_NatGenet_2025_GeneSets-QC.csv')}")
