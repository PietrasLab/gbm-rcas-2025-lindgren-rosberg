# 09-4 — Nomura (GBM-CARE): marker-anchored recovery of the malignant Neopl-ACR state
#        GFAP-high-gated NEOPLASTIC cells, characterised INDEPENDENTLY of our signature.
#
# INPUT: nomura_validation_seurat.rds is a STRATIFIED, STATE-BALANCED DOWNSAMPLE of the
#   GBM-CARE cohort (Nomura 2025 / Spitzer 2025; 121 samples, ~429k annotated nuclei), built
#   to ~120k nuclei: ALL astrocytes retained, malignant nuclei balanced by sample x State,
#   non-malignant by sample x CellType (see build_nomura_validation_object.R; STAR "Public
#   external ... GBM datasets"). The 72,158 malignant nuclei scored here are that subset
#   (59 patients contribute to the per-patient GFAP correlation). NB: state-balancing is
#   applied to the malignant compartment BEFORE this GFAP-vs-program correlation.
#
# PROMOTED PANEL: the per-patient GFAP-correlation panel (gCor) is manuscript Figure S6E
# (Fig-7 supplement). The other panels here remain reviewers-only / provisional.
#
# WHY: external "GFAP+ tumour-cell" / reactive-astrocyte-like evidence.
#   Rather than scoring OUR signature and asking whether it rises (self-referential),
#   we anchor on a single canonical marker (GFAP) to DEFINE candidate GFAP-high neoplastic cells,
#   then ask independently whether those cells carry the broader reactive/AC-MES program
#   (our Astrocyte-R & Neopl-ACR + Neftel/Nomura/Richards). Parallels mouse GFAP+RFP+ (Fig 7B).
#
# GATE: THREE GFAP bins for neoplastic cells, keyed to the cohort's OWN non-malignant astrocytes
#   (the true GFAP reference): GFAP-high (>=90th pct of astrocyte GFAP; "very high"), GFAP-mid
#   (50-90th pct; elevated but not highest), GFAP-low (<50th pct). The high bin is a high-confidence
#   tail (~5% of tumour cells); the mid bin shows the intermediate/dose-response step. The
#   recurrence-fraction test and the sensitivity sweep (75/90/95th) key on the HIGH bin.
#
# INDEPENDENCE: GFAP is a member of both Astrocyte-R and Neopl-ACR, so GFAP is REMOVED from EVERY
#   signature before scoring (leave-one-out on the gate gene) -> enrichment reflects the OTHER genes.
#
# TERMINOLOGY: Astrocyte-R / ACR = NON-neoplastic reactive astrocyte; Neopl-ACR (NACR) = the
#   MALIGNANT look-alike. This script characterises the MALIGNANT compartment (sibling of 09-3,
#   which tests the NON-neoplastic compartment). Richards dev/injury = the FULL (long) lists
#   (InHouse_BulkRNAseq_2019_*GSC, ~4k genes; = Richards_*_2021 full in Fig 4).
#
# HONEST EXPECTATION: tests the IDENTITY of the state. Recurrence-fraction test is expected NULL,
#   consistent with the malignant Neopl-ACR score being null Primary->Recurrent in this cohort.
#
# Determinism: set.seed(42); AddModuleScore(seed = 42).

suppressMessages({library(Seurat); library(dplyr); library(readr); library(tidyr)
                  library(tibble); library(ggplot2)})
set.seed(42)

R       <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")   # publish repo root
NOMURA  <- Sys.getenv("NOMURA_VALIDATION_DIR", unset = file.path(R, "data/external/nomura-gbm-care"))
OBJ     <- file.path(NOMURA, "nomura_validation_seurat.rds")
SIG_DIR <- file.path(R, "references/signatures")
GS_DIR  <- file.path(R, "references/genesets")
out_dir <- file.path(R, "results/09-4-nomura-gfap-neoplastic")
fig_dir <- file.path(R, "manuscript-figures/figure-5") # gCor panel = manuscript Fig S6E; other panels reviewers-only.
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GATE_GENE  <- "GFAP"     # human object -> HUGO uppercase
GATE_PCT   <- 0.90       # HIGH cutoff = 90th pct of GFAP across non-malig astrocytes (very high)
GATE_MID_PCT <- 0.50     # MID cutoff = 50th pct (median) -> mid bin = 50-90th pct (elevated, not highest)
SWEEP_PCTS <- c(0.75, 0.90, 0.95)
MIN_GRP    <- 10         # min GFAP-high and GFAP-low neoplastic cells per patient for the paired test
MIN_MALIG  <- 50         # min malignant cells per patient x timepoint for the fraction test
HI <- "GFAP-high neopl"; MID <- "GFAP-mid neopl"; LO <- "GFAP-low neopl"
PR_COL     <- c("Primary" = "#6BAED6", "Recurrent" = "#D9A55A")

## ---- signature panel (mirrors 09-1 build_panel_full, human/hs) --------------------------------
read_sig <- function(stem, suff = "hs") {
  d <- read_csv(file.path(SIG_DIR, sprintf("%s_original_%s.csv", stem, suff)), show_col_types = FALSE)
  if ("original" %in% names(d)) d$gene[d$original] else d$gene
}
build_panel <- function() {
  gs  <- readRDS(file.path(GS_DIR, "Richards_NatCancer_2021_GeneSets.rds"))
  nom <- readRDS(file.path(GS_DIR, "Nomura_NatGenet_2025_GeneSets.rds"))
  ns  <- attr(nom, "mp_state_map")
  ours_neo <- list(
    "Neopl-ACR (ours)"    = read_sig("neopl_ACR"), "Neopl-ECM (ours)" = read_sig("neopl_ECM"),
    "Neopl-COP (ours)"    = read_sig("neopl_COP"), "Neopl-NC (ours)"  = read_sig("neopl_NC"),
    "Neopl-OPC (ours)"    = read_sig("neopl_OPC"),
    "Neopl-CC/Bulk (ours)"= unique(c(read_sig("neopl_CC-I"), read_sig("neopl_CC-II"),
                                     read_sig("neopl_CC-III"), read_sig("neopl_Bulk"))))
  ours_ar <- list("Astrocyte R (ours)" = read_sig("astrocyte_R"))
  neftel_keep <- c("Neftel_Cell_2019_AC","Neftel_Cell_2019_MES1","Neftel_Cell_2019_MES2",
                   "Neftel_Cell_2019_NPC1","Neftel_Cell_2019_NPC2","Neftel_Cell_2019_OPC")
  neftel <- gs[neftel_keep]; names(neftel) <- gsub("Neftel_Cell_2019_","Neftel ",names(neftel))
  # FULL (long) Richards lists = InHouse_BulkRNAseq_2019_*GSC (~4k genes; = Richards_*_2021 in Fig 4).
  richards <- gs[c("InHouse_BulkRNAseq_2019_DevelopmentalGSC","InHouse_BulkRNAseq_2019_InjuryResponseGSC")]
  names(richards) <- c("Richards Developmental (full)","Richards Injury-Response (full)")
  nom_states <- c("OPC-like","AC-like","Hypoxia","MES-like","NPC-like","GPC-like",
                  "NEU-like","Stress","Cilia-like","Cell cycle")
  nom_sel <- lapply(nom_states, function(st)
    unique(unlist(nom[paste0(attr(nom,"prefix"), ns$mp[ns$state == st])], use.names = FALSE)))
  names(nom_sel) <- paste0("Nomura ", nom_states)
  sigs <- c(ours_ar, ours_neo, neftel, richards, nom_sel)
  meta <- bind_rows(
    tibble(label = names(ours_ar),  source = "Ours (non-neopl ref)"),
    tibble(label = names(ours_neo), source = "Ours (neoplastic)"),
    tibble(label = names(neftel),   source = "Neftel 2019"),
    tibble(label = names(richards), source = "Richards 2021"),
    tibble(label = names(nom_sel),  source = "Nomura 2025"))
  sigs <- lapply(sigs, function(g) setdiff(g, GATE_GENE))    # leave-one-out on the gate gene
  keep <- vapply(sigs, function(g) length(g) >= 5, logical(1))
  list(sigs = sigs[keep], meta = meta[keep, ])
}

## ---- load object + labels ---------------------------------------------------------------------
so <- readRDS(OBJ)
if ("RNA" %in% Assays(so)) DefaultAssay(so) <- "RNA"
stopifnot(GATE_GENE %in% rownames(so))
so$Timepoint <- as.character(so$Timepoint)
so$PvR       <- ifelse(so$Timepoint == "Primary", "Primary", "Recurrent")
so$is_astro  <- grepl("astro", tolower(as.character(so$CellType)))
so$is_malig  <- so$Malignant %in% TRUE
TOTAL <- ncol(so)   # denominator for "% of all cells" labels on the violins
cat(sprintf("cells %d | astrocytes %d | malignant %d\n", TOTAL, sum(so$is_astro), sum(so$is_malig)))

dat <- tryCatch(GetAssayData(so, layer = "data"), error = function(e) NULL)
if (is.null(dat) || nrow(dat) == 0 || length(dat@x) == 0) {
  cat("normalising (data layer empty)\n"); so <- NormalizeData(so, verbose = FALSE)
}
gfap_vec <- as.numeric(GetAssayData(so, layer = "data")[GATE_GENE, ])
so$GFAP_expr <- gfap_vec

## ---- gate: THREE GFAP bins keyed to NON-MALIGNANT astrocytes ---------------------------------
ref_gfap  <- gfap_vec[so$is_astro & !so$is_malig]
cutoff_of <- function(pct) as.numeric(quantile(ref_gfap, pct))
CUT     <- cutoff_of(GATE_PCT)       # HIGH cutoff (also the reference-plot dashed line)
CUT_MID <- cutoff_of(GATE_MID_PCT)   # MID cutoff
cat(sprintf("non-malig astrocytes: %d (%.1f%% GFAP+) | cutoffs: mid(%.0fth)=%.3f  high(%.0fth)=%.3f\n",
            length(ref_gfap), 100 * mean(ref_gfap > 0),
            100 * GATE_MID_PCT, CUT_MID, 100 * GATE_PCT, CUT))

so$gfap_bin <- ifelse(so$GFAP_expr >= CUT, HI, ifelse(so$GFAP_expr >= CUT_MID, MID, LO))
malig <- so[, so$is_malig]
malig$gfap_grp <- malig$gfap_bin
cat(sprintf("malignant cells: %d | GFAP-high %d (%.1f%%) | GFAP-mid %d (%.1f%%) | GFAP-low %d\n",
            ncol(malig), sum(malig$gfap_grp == HI), 100 * mean(malig$gfap_grp == HI),
            sum(malig$gfap_grp == MID), 100 * mean(malig$gfap_grp == MID), sum(malig$gfap_grp == LO)))

## ---- score the (GFAP-excluded) panel on ALL malignant cells ----------------------------------
panel <- build_panel()
feats <- rownames(malig)
present <- lapply(panel$sigs, function(g) intersect(g, feats))
ok <- vapply(present, length, integer(1)) >= 5
present <- present[ok]; pmeta <- panel$meta[ok, ]
cat(sprintf("scoring %d signatures on %d malignant cells\n", length(present), ncol(malig)))
for (i in seq_along(present)) {
  malig <- AddModuleScore(malig, features = list(present[[i]]), name = paste0("S_", i), seed = 42)
  malig[[names(present)[i]]] <- malig[[paste0("S_", i, "1")]]
}
## also score the REFINED (Astrocyte-R-overlap removed) Neopl-ACR for the GFAP correlation panel
acr_ref <- setdiff(readr::read_csv(file.path(SIG_DIR, "neopl_ACR_refined_hs.csv"), show_col_types = FALSE)$gene, GATE_GENE)
malig <- AddModuleScore(malig, features = list(intersect(acr_ref, rownames(malig))), name = "ACRref_", seed = 42)
malig[["Neopl-ACR (refined)"]] <- malig[["ACRref_1"]]

## ---- GFAP-high vs GFAP-low neoplastic: PATIENT-LEVEL paired test (honest unit = patient) -------
md <- malig@meta.data
sig_names <- names(present)
per_patient <- lapply(sig_names, function(sn) {
  d <- md %>% transmute(Patient, gfap_grp, val = .data[[sn]]) %>%
    group_by(Patient, gfap_grp) %>% summarise(n = n(), m = mean(val), .groups = "drop")
  whi <- d %>% filter(gfap_grp == HI); wlo <- d %>% filter(gfap_grp == LO)
  w <- inner_join(whi, wlo, by = "Patient", suffix = c("_hi", "_lo")) %>%
    filter(n_hi >= MIN_GRP, n_lo >= MIN_GRP)
  if (nrow(w) < 3) return(NULL)
  w$delta <- w$m_hi - w$m_lo
  tibble(signature = sn, n_paired = nrow(w),
         mean_hi = mean(w$m_hi), mean_lo = mean(w$m_lo),
         mean_delta = mean(w$delta), se_delta = sd(w$delta) / sqrt(nrow(w)), n_up = sum(w$delta > 0),
         p = wilcox.test(w$m_hi, w$m_lo, paired = TRUE)$p.value)
})
sig_cmp <- bind_rows(per_patient) %>%
  left_join(pmeta, by = c("signature" = "label")) %>%
  mutate(padj = p.adjust(p, "BH")) %>% arrange(desc(mean_delta))
cat("\n=== GFAP-high vs GFAP-low neoplastic (patient-paired mean delta; GFAP excluded from sigs) ===\n")
print(as.data.frame(sig_cmp %>% select(signature, source, n_paired, mean_delta, p, padj)))

## ---- recurrence: GFAP-high fraction of malignant cells per patient x PvR (expected NULL) -------
frac <- md %>% group_by(Patient, PvR) %>%
  summarise(n_malig = n(), gfap_hi_frac = mean(gfap_grp == HI), .groups = "drop") %>%
  filter(n_malig >= MIN_MALIG)
paired_frac <- {
  w <- frac %>% select(Patient, PvR, gfap_hi_frac) %>%
    pivot_wider(names_from = PvR, values_from = gfap_hi_frac) %>%
    filter(!is.na(Primary) & !is.na(Recurrent))
  if (nrow(w) >= 3) {
    w$delta <- w$Recurrent - w$Primary
    tibble(metric = "gfap_hi_frac", n_paired = nrow(w), mean_delta = mean(w$delta),
           n_up = sum(w$delta > 0), n_down = sum(w$delta < 0),
           p = wilcox.test(w$Recurrent, w$Primary, paired = TRUE)$p.value)
  } else tibble(metric = "gfap_hi_frac", n_paired = nrow(w), mean_delta = NA_real_,
                n_up = NA_integer_, n_down = NA_integer_, p = NA_real_)
}
cat("\n=== GFAP-high neoplastic fraction Primary->Recurrent (paired) ===\n")
print(as.data.frame(paired_frac))

## ---- gate sensitivity sweep -------------------------------------------------------------------
sweep <- bind_rows(lapply(SWEEP_PCTS, function(p) {
  cut_p <- cutoff_of(p)
  hi <- md[["Astrocyte R (ours)"]][malig$GFAP_expr >= cut_p]
  lo <- md[["Astrocyte R (ours)"]][malig$GFAP_expr <  cut_p]
  tibble(pct = p, cutoff = cut_p, n_hi = sum(malig$GFAP_expr >= cut_p),
         n_lo = sum(malig$GFAP_expr < cut_p), astroR_delta = mean(hi) - mean(lo))
}))
cat("\n=== gate sensitivity (Astrocyte-R delta GFAP-high vs GFAP-low, per cutoff) ===\n")
print(as.data.frame(sweep))

## ---- SUMMARY tables for the violins (compact; NOT per-cell, to keep repo small) ---------------
focus <- intersect(c("Astrocyte R (ours)","Neopl-ACR (ours)","Neftel AC","Nomura AC-like",
                     "Neftel MES1","Richards Injury-Response (full)"), names(present))
vln_summary <- bind_rows(lapply(focus, function(sn) {
  md %>% select(gfap_grp, val = all_of(sn)) %>% group_by(gfap_grp) %>%
    summarise(signature = sn, n = n(), mean = mean(val), sd = sd(val),
              q25 = quantile(val, .25), median = median(val), q75 = quantile(val, .75),
              .groups = "drop")
})) %>% select(signature, gfap_grp, n, mean, sd, q25, median, q75)

## ---- write analysis outputs (results/) --------------------------------------------------------
write_csv(sig_cmp,     file.path(out_dir, "09-4-Nomura-gfap-neopl-signature-comparison.csv"))
write_csv(frac,        file.path(out_dir, "09-4-Nomura-gfap-neopl-fraction-perpatient.csv"))
write_csv(paired_frac, file.path(out_dir, "09-4-Nomura-gfap-neopl-fraction-stats.csv"))
write_csv(sweep,       file.path(out_dir, "09-4-Nomura-gfap-gate-sensitivity.csv"))
write_csv(vln_summary, file.path(out_dir, "09-4-Nomura-gfap-focus-violin-summary.csv"))
saveRDS(list(sig_cmp = sig_cmp, frac = frac, paired_frac = paired_frac, sweep = sweep,
             vln_summary = vln_summary, cutoff = CUT, cutoff_mid = CUT_MID, gate_pct = GATE_PCT,
             n_malig = ncol(malig), gate_gene = GATE_GENE),
        file.path(out_dir, "09-4-Nomura-gfap-neopl.rds"))

## ---- provisional supp panels (manuscript-figures/) + source-data companions -------------------
# append per-group cell count + % of all cells to a factor's level labels
add_lab <- function(f, total) {
  f <- as.factor(f); n <- as.integer(table(f)[levels(f)])
  levels(f) <- sprintf("%s\n(n=%s, %.1f%%)", levels(f), format(n, big.mark = ","), 100 * n / total)
  f
}

# Panel A: per-signature GFAP-high minus GFAP-low delta across the panel (patient-paired means)
cmpA <- sig_cmp %>% mutate(sig = factor(signature, levels = rev(signature)),
                           star = as.character(cut(padj, c(-Inf,1e-3,1e-2,5e-2,Inf), c("***","**","*",""))))
gA <- ggplot(cmpA, aes(mean_delta, sig, fill = source)) +
  geom_col(width = .72) +
  geom_errorbarh(aes(xmin = mean_delta - se_delta, xmax = mean_delta + se_delta),
                 height = .3, linewidth = .3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  geom_text(aes(x = mean_delta + se_delta, label = star), hjust = -0.2, size = 3) +
  labs(title = "Nomura: GFAP-high vs GFAP-low neoplastic cells (patient-paired mean delta +/- SE)",
       subtitle = "GFAP excluded from every signature (leave-one-out). Stars = BH-adjusted paired Wilcoxon. PROVISIONAL / reviewers-only.",
       x = "mean(GFAP-high neopl) - mean(GFAP-low neopl)", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"), legend.position = "bottom")
# [reviewer-only panel removed from publication repo]

# Panel B (GFAP-vs-Neopl-ACR scatter) DROPPED: the promoted manuscript panel is the
# correlation-bar panel only (Figure S6E). It now carries BOTH Neopl-ACR (original)
# and (refined) bars, so the original-vs-refined robustness is already shown there and the
# standalone scatter is redundant. Refined Neopl-ACR is still scored above (feeds the correlation panel).

# Panel A' (CONTINUOUS): GFAP as a continuous vector correlated with each program per malignant cell,
# per-patient Spearman aggregated -> mean +/- SE (replaces the high/low binning of Panel A).
# includes BOTH Neopl-ACR (original) and Neopl-ACR (refined) so this panel is self-consistent with
# the scatter (Panel B); "Neopl-ACR (ours)" is relabelled "(original)" so the two versions read cleanly.
sig_cols <- c(names(present), "Neopl-ACR (refined)")
pmeta2   <- bind_rows(pmeta, tibble(label = "Neopl-ACR (refined)", source = "Ours (neoplastic)"))
cor_pp <- bind_rows(lapply(sig_cols, function(sn) {
  md %>% group_by(Patient) %>%
    summarise(rho = suppressWarnings(cor(GFAP_expr, .data[[sn]], method = "spearman")), n = n(), .groups = "drop") %>%
    filter(n >= MIN_GRP, !is.na(rho)) %>%
    summarise(signature = sn, n_pat = n(), mean_rho = mean(rho), se_rho = sd(rho) / sqrt(n()))
})) %>% left_join(pmeta2, by = c("signature" = "label")) %>%
  mutate(signature = ifelse(signature == "Neopl-ACR (ours)", "Neopl-ACR (original)", signature)) %>%
  arrange(desc(mean_rho)) %>%
  mutate(sig = factor(signature, levels = rev(signature)))
gCor <- ggplot(cor_pp, aes(mean_rho, sig, fill = source)) +
  geom_col(width = .72) +
  geom_errorbarh(aes(xmin = mean_rho - se_rho, xmax = mean_rho + se_rho), height = .3, linewidth = .3, colour = "grey25") +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey40") +
  labs(title = "GBM-CARE: per-cell GFAP correlation with each program (continuous)",
       subtitle = "per-patient Spearman rho(GFAP, module score) across malignant cells (state-balanced downsample), mean +/- SE. GFAP excluded from every signature.",
       x = "Spearman rho (GFAP vs module score)", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"), legend.position = "bottom")
ggsave(file.path(fig_dir, "Figure-S5-nomura-gfap-correlation-per-signature.pdf"), gCor, width = 7.5, height = 5)
write_csv(cor_pp %>% select(signature, source, n_pat, mean_rho, se_rho),
          file.path(fig_dir, "Figure-S5-nomura-gfap-correlation-per-signature-sourcedata.csv"))

# Panel B: three-bin focus violins for the key AC/reactive signatures (dose-response)
vln <- md %>% select(gfap_grp, all_of(focus)) %>%
  pivot_longer(-gfap_grp, names_to = "signature", values_to = "score")
vln$signature <- factor(vln$signature, levels = focus)
vln$gfap_grp  <- add_lab(factor(vln$gfap_grp, levels = c(LO, MID, HI)), TOTAL)
gB <- ggplot(vln, aes(gfap_grp, score, fill = gfap_grp)) +
  geom_violin(scale = "width", linewidth = .2) +
  geom_boxplot(width = .12, outlier.shape = NA, linewidth = .25, fill = "white") +
  facet_wrap(~signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("grey70", "#E8A33D", "#C0504D"), guide = "none") +  # low, mid, high
  labs(title = "Reactive / AC-MES programs across GFAP bins in neoplastic cells (Nomura)",
       subtitle = "Per malignant cell; GFAP-low/mid/high = <50th / 50-90th / >=90th pct of astrocyte GFAP. GFAP excluded from sigs. Labels: n, % of all cells. PROVISIONAL / reviewers-only.",
       x = NULL, y = "module score") +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"),
        strip.text = element_text(size = 7.5, face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))
# [reviewer-only panel removed from publication repo]

# Panel C: paired GFAP-high neoplastic fraction Primary->Recurrent (expected null)
fr <- frac; fr$PvR <- factor(fr$PvR, levels = c("Primary","Recurrent"))
gC <- ggplot(fr, aes(PvR, gfap_hi_frac, group = Patient)) +
  geom_line(alpha = .4, colour = "grey55") +
  geom_point(aes(colour = PvR), size = 1.3, alpha = .7) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "#F86814", linewidth = 1.1) +
  scale_colour_manual(values = PR_COL, guide = "none") +
  labs(title = "GFAP-high neoplastic fraction Primary->Recurrent (Nomura, paired)",
       subtitle = sprintf("paired n=%s, p=%s. PROVISIONAL / reviewers-only.",
                          paired_frac$n_paired, signif(paired_frac$p, 2)),
       x = NULL, y = "GFAP-high fraction of malignant cells") +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"))
# [reviewer-only panel removed from publication repo]

## ---- REFERENCE panels: GFAP expression across the WHOLE dataset, by cell group ----------------
allmd <- so@meta.data
allmd$GFAP <- so$GFAP_expr
neopl_lab <- c("GFAP-high neopl" = "Neoplastic GFAP-high", "GFAP-mid neopl" = "Neoplastic GFAP-mid",
               "GFAP-low neopl" = "Neoplastic GFAP-low")
allmd$group  <- ifelse(allmd$is_malig, neopl_lab[allmd$gfap_bin], as.character(allmd$CellType))
allmd$group5 <- ifelse(allmd$is_malig, neopl_lab[allmd$gfap_bin],
                       ifelse(allmd$is_astro, "Non-neopl astrocyte", "Non-neopl non-astrocyte"))

gfap_expr_summary <- allmd %>% group_by(group) %>%
  summarise(n = n(), pct_of_all = 100 * n() / TOTAL, pct_gfap_high = 100 * mean(GFAP >= CUT),
            mean = mean(GFAP), sd = sd(GFAP), q25 = quantile(GFAP, .25),
            median = median(GFAP), q75 = quantile(GFAP, .75), .groups = "drop") %>%
  arrange(desc(median))
write_csv(gfap_expr_summary, file.path(out_dir, "09-4-Nomura-gfap-expression-by-group-summary.csv"))
bin_summary <- allmd %>% group_by(group5) %>%
  summarise(n = n(), pct_of_all = 100 * n() / TOTAL, pct_gfap_high = 100 * mean(GFAP >= CUT),
            mean = mean(GFAP), q25 = quantile(GFAP, .25), median = median(GFAP),
            q75 = quantile(GFAP, .75), .groups = "drop")
write_csv(bin_summary, file.path(out_dir, "09-4-Nomura-gfap-expression-by-bin-summary.csv"))

allmd$fillcat <- ifelse(allmd$group == "Neoplastic GFAP-high", "Neoplastic GFAP-high",
                 ifelse(allmd$group == "Neoplastic GFAP-mid",  "Neoplastic GFAP-mid",
                 ifelse(allmd$group == "Neoplastic GFAP-low",  "Neoplastic GFAP-low",
                 ifelse(grepl("astro", tolower(allmd$group)), "Astrocyte (non-neopl)", "Other non-neopl"))))
fill_cols <- c("Neoplastic GFAP-high" = "#C0504D", "Neoplastic GFAP-mid" = "#E8A33D",
               "Neoplastic GFAP-low" = "grey60", "Astrocyte (non-neopl)" = "#2E8B57",
               "Other non-neopl" = "grey82")

# Panel D: per-cell-type GFAP violin
neo_grps <- c("Neoplastic GFAP-high", "Neoplastic GFAP-mid", "Neoplastic GFAP-low")
nn_order <- gfap_expr_summary %>% filter(!group %in% neo_grps) %>% arrange(desc(median)) %>% pull(group)
allmd$fillcat <- factor(allmd$fillcat, levels = names(fill_cols))
allmd$group   <- add_lab(factor(allmd$group, levels = c(neo_grps, nn_order)), TOTAL)
gD <- ggplot(allmd, aes(group, GFAP, fill = fillcat)) +
  geom_violin(scale = "width", linewidth = .15) +
  geom_boxplot(width = .12, outlier.shape = NA, linewidth = .2, fill = "white") +
  geom_hline(yintercept = c(CUT_MID, CUT), linetype = "dashed", colour = "#C0504D", linewidth = .35) +
  scale_fill_manual(values = fill_cols, name = NULL) +
  labs(title = "GFAP expression by cell group (Nomura, whole dataset)",
       subtitle = sprintf("dashed = mid / high cutoffs (%.2f / %.2f = 50th / 90th pct of astrocyte GFAP). Labels: n, %% of all cells. PROVISIONAL / reviewers-only.", CUT_MID, CUT),
       x = NULL, y = "GFAP (log-normalised)") +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"),
        axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
# [reviewer-only panel removed from publication repo]

# Panel E: collapsed GFAP violin — three neoplastic bins vs non-neoplastic reference
allmd$group5 <- add_lab(factor(allmd$group5, levels = c("Neoplastic GFAP-high", "Neoplastic GFAP-mid",
                       "Neoplastic GFAP-low", "Non-neopl astrocyte", "Non-neopl non-astrocyte")), TOTAL)
gE <- ggplot(allmd, aes(group5, GFAP, fill = group5)) +
  geom_violin(scale = "width", linewidth = .2) +
  geom_boxplot(width = .12, outlier.shape = NA, linewidth = .2, fill = "white") +
  geom_hline(yintercept = c(CUT_MID, CUT), linetype = "dashed", colour = "#C0504D", linewidth = .35) +
  scale_fill_manual(values = c("#C0504D", "#E8A33D", "grey60", "#2E8B57", "grey82"), guide = "none") +
  labs(title = "GFAP expression: three neoplastic GFAP bins vs non-neoplastic reference (Nomura)",
       subtitle = sprintf("dashed = mid / high cutoffs (%.2f / %.2f). Labels: n, %% of all cells. PROVISIONAL / reviewers-only.", CUT_MID, CUT),
       x = NULL, y = "GFAP (log-normalised)") +
  theme_bw(base_size = 9) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"),
        axis.text.x = element_text(angle = 20, hjust = 1))
# [reviewer-only panel removed from publication repo]

cat("\nDONE 09-4 ->\n  results:", out_dir, "\n  panels: ", fig_dir, "\n")
