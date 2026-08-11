# 09-1 — score & correlate one external dataset (self-contained; v2)
#
# Re-scores a dataset FROM RAW (so v2 is portable / standalone) and writes the signature x
# signature Spearman correlation matrices that the joint figure (13-3) assembles:
#   results/31-external-validation/31-corr-<ds>-neoplastic.csv   full neoplastic panel (26 sigs) on MALIGNANT cells  [MAIN]
#   results/31-external-validation/31-corr-<ds>-allcells.csv     full 40-sig panel on a per-cell-type-capped ALL-cells set [SUPP]
# + results/31-external-validation/31-scores-<ds>.rds            per-cell score matrices + compartment meta (hand-back)
#
# Usage:  Rscript 09-1-external-validation-score.R <Soni|GBmap|Sussman|Suter>
# Nomura is scored HPC-side (data is HPC-only); its correlation is derived from the pulled
# C18c-Nomura-scores.rds by 09-1b (documented external input). Determinism: set.seed(42).

suppressMessages({library(Seurat); library(dplyr); library(readr); library(tibble); library(Matrix)})
set.seed(42)
R  <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")  # publish repo root (data/external, user-supplied)
V  <- R
out_dir <- file.path(R, "results/31-external-validation")   # heavy/regenerable (git-ignored)
res_dir <- file.path(R, "results/31-external-validation")  # small portable corr CSVs (tracked)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
SIG_DIR <- file.path(R, "references/signatures")

neftel_keep   <- c("Neftel_Cell_2019_AC","Neftel_Cell_2019_MES1","Neftel_Cell_2019_MES2",
                   "Neftel_Cell_2019_NPC1","Neftel_Cell_2019_NPC2","Neftel_Cell_2019_OPC")
richards_keep <- c("Developmental_Richards","Injury_Response_Richards")
# Nomura references use the COALESCED 9 states (+ cell cycle), built from mp_state_map in
# build_panel_full — NOT the 15 raw MPs (see there).

read_sig <- function(stem, suff) { d <- read_csv(file.path(SIG_DIR, sprintf("%s_original_%s.csv", stem, suff)), show_col_types = FALSE)
  if ("original" %in% names(d)) d$gene[d$original] else d$gene }
build_panel_full <- function(species) {
  suff <- if (species=="mouse") "mm" else "hs"
  gs  <- readRDS(file.path(V, sprintf("references/genesets/Richards_NatCancer_2021_GeneSets%s.rds", if (species=="mouse") "_Mouse" else "")))
  nom <- readRDS(file.path(V, sprintf("references/genesets/Nomura_NatGenet_2025_GeneSets%s.rds", if (species=="mouse") "_Mouse" else "")))
  ns  <- attr(nom, "mp_state_map")
  ours_neo <- list(
    "Neopl-ACR (ours)"=read_sig("neopl_ACR",suff), "Neopl-ECM (ours)"=read_sig("neopl_ECM",suff),
    "Neopl-COP (ours)"=read_sig("neopl_COP",suff), "Neopl-NC (ours)"=read_sig("neopl_NC",suff),
    "Neopl-OPC (ours)"=read_sig("neopl_OPC",suff),
    "Neopl-CC/Bulk (ours)"=unique(c(read_sig("neopl_CC-I",suff),read_sig("neopl_CC-II",suff),read_sig("neopl_CC-III",suff),read_sig("neopl_Bulk",suff))))
  nn_stems <- c("Astrocyte","Microglia","Macrophage","OPC-COP-OLG","Neural","Endothelial","Mural","Fibroblast","Choroid","Ependymal","Dendritic","Neutrophil","NKTB")
  ours_nn <- c(list("Astrocyte R (ours)"=read_sig("astrocyte_R",suff)),
               setNames(lapply(nn_stems, function(s) read_sig(paste0("nonneopl_",s),suff)), paste0("NN-",nn_stems," (ours)")))
  neftel <- gs[neftel_keep];   names(neftel)   <- gsub("Neftel_Cell_2019_","Neftel ",names(neftel))
  richards <- gs[richards_keep]; names(richards) <- c("Richards Developmental","Richards Injury-Response")
  # Adopt Nomura's FINAL MP->state decisions: coalesce NEU-like (MP9+MP14) and Stress (MP10+MP15)
  # by GENE UNION; single-MP states as-is; DROP unassigned (MP1 RP / MP11 Doublet / MP12 LQ).
  nom_states <- c("OPC-like","AC-like","Hypoxia","MES-like","NPC-like","GPC-like",
                  "NEU-like","Stress","Cilia-like","Cell cycle")
  nom_sel <- lapply(nom_states, function(st)
    unique(unlist(nom[paste0(attr(nom,"prefix"), ns$mp[ns$state == st])], use.names = FALSE)))
  names(nom_sel) <- paste0("Nomura ", nom_states)
  sigs <- c(ours_neo, ours_nn, neftel, richards, nom_sel)
  meta <- bind_rows(
    tibble(label=names(ours_neo), source="Ours",          compartment="neoplastic"),
    tibble(label=names(ours_nn),  source="Ours",          compartment="non-neoplastic"),
    tibble(label=names(neftel),   source="Neftel 2019",   compartment="neoplastic"),
    tibble(label=names(richards), source="Richards 2021", compartment="neoplastic"),
    tibble(label=names(nom_sel),  source="Nomura 2025",   compartment="neoplastic"))
  keep <- vapply(sigs, function(g) length(g)>=5, logical(1)); list(sigs=sigs[keep], meta=meta[keep,])
}
score_panel <- function(obj, panel, map_fn=function(g) g) {
  feats <- rownames(obj); present <- lapply(panel$sigs, function(g) intersect(map_fn(g), feats))
  ok <- vapply(present, length, integer(1))>=5; present <- present[ok]; meta <- panel$meta[ok,]
  cat(sprintf("  scoring %d sigs on %d cells\n", length(present), ncol(obj)))
  for (i in seq_along(present)) { obj <- AddModuleScore(obj, features=list(present[[i]]), name=paste0("S_",i), seed=42)
    obj[[names(present)[i]]] <- obj[[paste0("S_",i,"1")]] }
  list(score_mat=as.matrix(obj@meta.data[, names(present), drop=FALSE]), meta=meta)
}
ds_by <- function(obj, col, n) { g <- obj@meta.data[[col]]
  idx <- unlist(lapply(split(seq_len(ncol(obj)), g), function(ii) if (length(ii)>n) sample(ii,n) else ii)); obj[, sort(idx)] }
write_corr <- function(sm, meta, stub) {
  cm <- cor(sm, method="spearman", use="pairwise.complete.obs")
  df <- as.data.frame(round(cm,4)) |> rownames_to_column("signature")
  # attach source/compartment for the row signature so 31-2 can facet/annotate
  df <- left_join(df, meta[,c("label","source","compartment")], by=c("signature"="label"))
  write_csv(df, file.path(res_dir, paste0(stub,".csv"))); cat("  wrote", stub, "\n")
}

CONFIGS <- list(
  Soni    = list(path=file.path(R,"data/external/Soni_GBM_020925.rds"), species="mouse", ensg=FALSE,
                 malignant=function(md) md$cell_type=="neoplastic cell", celltype_col="cell_type"),
  GBmap   = list(path=file.path(R,"data/external/Core-GBmap.rds"), species="human", ensg=TRUE,
                 prefilter=function(o) o[, o$assay=="10x 3' v2"], malignant=function(md) md$annotation_level_1=="Neoplastic", celltype_col="annotation_level_3"),
  Sussman = list(path=file.path(R,"data/external/Sussman_pediatric_GBM_2026/Neoplastic_snRNAseq_020726.rds"),
                 species="human", ensg=FALSE, malignant=function(md) rep(TRUE,nrow(md)), celltype_col=NULL),
  Suter   = list(path=file.path(R,"data/external/Suter_human_GBM_2026/GSE229779_seuratObj.RDS"),
                 species="human", ensg=FALSE, malignant=function(md) md$CellType=="Neoplastic", celltype_col="CellType"),
  # Nomura / GBM-CARE: rebuilt LOCALLY from GSE274546 + GBM-CARE-WT annotations (see
  # data/external/Nomura_human_GBM_2026), so it now scores identically to the others -
  # replaces the 09-1b interim coalesced proxy for 31-corr-Nomura-neoplastic.csv.
  Nomura  = list(path=Sys.getenv("NOMURA_VALIDATION_OBJ",
                   unset=file.path(R, "data/external/nomura-gbm-care/nomura_validation_seurat.rds")),
                 species="human", ensg=FALSE, malignant=function(md) md$Malignant %in% TRUE, celltype_col="CellType")
)

run <- function(name) {
  cfg <- CONFIGS[[name]]; stopifnot(!is.null(cfg)); cat("\n===== 09-1:", name, "=====\n")
  so <- readRDS(cfg$path); if (!is.null(cfg$prefilter)) so <- cfg$prefilter(so); cat("cells:", ncol(so), "\n")
  map_fn <- function(g) g
  if (isTRUE(cfg$ensg)) { fd <- readRDS(file.path(R,"data/external/fdata_full.rds"))
    s2e <- setNames(fd$ensembl_gene_id, fd$external_gene_name_ensembl); map_fn <- function(g) unname(na.omit(s2e[g])) }
  panel <- build_panel_full(cfg$species)
  # MAIN: neoplastic panel on malignant cells
  malig <- so[, cfg$malignant(so@meta.data)]; cat("malignant:", ncol(malig), "\n")
  pN <- list(sigs=panel$sigs[panel$meta$label[panel$meta$compartment=="neoplastic"]],
             meta=panel$meta[panel$meta$compartment=="neoplastic",])
  scN <- score_panel(malig, pN, map_fn)
  write_corr(scN$score_mat, scN$meta, paste0("31-corr-",name,"-neoplastic"))
  # By-driver-model breakdown (Soni spans PDGFB / Nf1 / EGFRvIII): mean neoplastic-panel score
  # per model over malignant cells -> portable table for the Fig S4 by-model panel.
  if (name=="Soni" && "groupid" %in% colnames(malig@meta.data)) {
    grp  <- malig@meta.data$groupid
    bym  <- sapply(split(seq_along(grp), grp), function(ii) colMeans(scN$score_mat[ii,,drop=FALSE]))
    bymdf <- as.data.frame(round(bym,4)) |> rownames_to_column("signature")
    bymdf <- left_join(bymdf, scN$meta[,c("label","source","compartment")], by=c("signature"="label"))
    write_csv(bymdf, file.path(res_dir, "31-Soni-mean-by-model.csv"))
    cat("  wrote 31-Soni-mean-by-model (models:", paste(names(table(grp)),collapse=", "), ")\n")
  }
  # SUPP: full panel on capped all-cells (if the object has non-malignant cells)
  scA <- NULL
  if (!is.null(cfg$celltype_col)) {
    aobj <- ds_by(so, cfg$celltype_col, 3000); scA <- score_panel(aobj, panel, map_fn)
    write_corr(scA$score_mat, scA$meta, paste0("31-corr-",name,"-allcells"))
  }
  saveRDS(list(neo=scN, all=scA, dataset=name), file.path(out_dir, paste0("31-scores-",name,".rds")))
  cat("DONE", name, "\n")
}
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>=1) run(args[1]) else cat("usage: Rscript 31-1-score-and-correlate.R <", paste(names(CONFIGS),collapse="|"), ">\n")
