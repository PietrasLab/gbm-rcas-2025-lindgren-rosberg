# Shared harmonized signature panel (species-aware) for the Visium HD / Heiland
# spatial correlation figures. Optionally swaps the Neopl-ACR / Astrocyte-R
# TARGETS to their REFINED (de-confounded) versions (refined when ACR &
# Astrocyte-R are scored side-by-side).
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble) })
.PANEL_PROJ <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
.PANEL_SIG  <- file.path(.PANEL_PROJ, "references/signatures")
.panel_orig <- function(stem, suff) { d <- read_csv(file.path(.PANEL_SIG, sprintf("%s_original_%s.csv", stem, suff)), show_col_types = FALSE); if ("original" %in% names(d)) d$gene[d$original] else d$gene }
.panel_refined <- function(stem, suff) read_csv(file.path(.PANEL_SIG, sprintf("%s_refined_%s.csv", stem, suff)), show_col_types = FALSE)$gene
.PANEL_NEFTEL   <- c("Neftel_Cell_2019_AC","Neftel_Cell_2019_MES1","Neftel_Cell_2019_MES2","Neftel_Cell_2019_NPC1","Neftel_Cell_2019_NPC2","Neftel_Cell_2019_OPC")
# FULL Richards (full correlates strongly, short is weak; Fig 4 uses full)
.PANEL_RICHARDS <- c("InHouse_BulkRNAseq_2019_DevelopmentalGSC","InHouse_BulkRNAseq_2019_InjuryResponseGSC")

build_panel_full <- function(species = c("mouse","human"), acr_targets = c("both","refined","original")) {
  species <- match.arg(species); acr_targets <- match.arg(acr_targets); suff <- if (species == "mouse") "mm" else "hs"
  gs  <- readRDS(file.path(.PANEL_PROJ, sprintf("references/genesets/Richards_NatCancer_2021_GeneSets%s.rds", if (species=="mouse") "_Mouse" else "")))
  nom <- readRDS(file.path(.PANEL_PROJ, sprintf("references/genesets/Nomura_NatGenet_2025_GeneSets%s.rds", if (species=="mouse") "_Mouse" else ""))); ns <- attr(nom, "mp_state_map")
  # ACR / Astrocyte-R TARGETS: original and/or refined (labelled explicitly). "both"
  # emits both variants so every correlation analysis can be evaluated with each.
  acr_o <- .panel_orig("neopl_ACR", suff);   acr_r <- .panel_refined("neopl_ACR", suff)
  ast_o <- .panel_orig("astrocyte_R", suff); ast_r <- .panel_refined("astrocyte_R", suff)
  acr_list <- switch(acr_targets,
    both     = list("Neopl-ACR (ours)"=acr_o, "Neopl-ACR (ours, refined)"=acr_r),
    refined  = list("Neopl-ACR (ours, refined)"=acr_r),
    original = list("Neopl-ACR (ours)"=acr_o))
  ast_list <- switch(acr_targets,
    both     = list("Astrocyte R (ours)"=ast_o, "Astrocyte R (ours, refined)"=ast_r),
    refined  = list("Astrocyte R (ours, refined)"=ast_r),
    original = list("Astrocyte R (ours)"=ast_o))
  ours_neo <- c(acr_list, list("Neopl-ECM (ours)"=.panel_orig("neopl_ECM",suff), "Neopl-COP (ours)"=.panel_orig("neopl_COP",suff),
    "Neopl-NC (ours)"=.panel_orig("neopl_NC",suff), "Neopl-OPC (ours)"=.panel_orig("neopl_OPC",suff),
    "Neopl-CC/Bulk (ours)"=unique(c(.panel_orig("neopl_CC-I",suff),.panel_orig("neopl_CC-II",suff),.panel_orig("neopl_CC-III",suff),.panel_orig("neopl_Bulk",suff)))))
  nn_stems <- c("Astrocyte","Microglia","Macrophage","OPC-COP-OLG","Neural","Endothelial","Mural","Fibroblast","Choroid","Ependymal","Dendritic","Neutrophil","NKTB")  # incl homeostatic Astrocyte
  ours_nn <- c(ast_list, setNames(lapply(nn_stems, function(s) .panel_orig(paste0("nonneopl_",s),suff)), paste0("NN-",nn_stems," (ours)")))
  neftel <- gs[.PANEL_NEFTEL]; names(neftel) <- gsub("Neftel_Cell_2019_","Neftel ",names(neftel))
  richards <- gs[.PANEL_RICHARDS]; names(richards) <- c("Richards Developmental","Richards Injury-Response")
  nom_states <- c("OPC-like","AC-like","Hypoxia","MES-like","NPC-like","GPC-like","NEU-like","Stress","Cilia-like","Cell cycle")
  nom_sel <- lapply(nom_states, function(st) unique(unlist(nom[paste0(attr(nom,"prefix"), ns$mp[ns$state == st])], use.names = FALSE))); names(nom_sel) <- paste0("Nomura ", nom_states)
  sigs <- c(ours_neo, ours_nn, neftel, richards, nom_sel)
  meta <- bind_rows(tibble(label=names(ours_neo), source="Ours", compartment="neoplastic"),
    tibble(label=names(ours_nn), source="Ours", compartment="non-neoplastic"),
    tibble(label=names(neftel), source="Neftel 2019", compartment="neoplastic"),
    tibble(label=names(richards), source="Richards 2021", compartment="neoplastic"),
    tibble(label=names(nom_sel), source="Nomura 2025", compartment="neoplastic"))
  keep <- vapply(sigs, function(g) length(g) >= 5, logical(1)); list(sigs=sigs[keep], meta=meta[keep,])
}
