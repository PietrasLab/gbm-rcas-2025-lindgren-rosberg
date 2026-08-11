# 16-3 — Figure 7: named composition niches (spatial map + enrichment + frequency)
# ============================================================================
# The composition niches (KNN-mean of RCTD weight vectors -> joint k-means, 12 niches)
# with data-driven descriptive names, rendered in the canonical spatial style. Three panels:
#   (1) named niche spatial map (spatial_bin_grid; normal niches paled, tumour/reactive saturated)
#   (2) cell-type enrichment heatmap per named niche (row z-scored mean RCTD weight)
#   (3) named niche frequency by condition (equal-width dodged bars incl. zero-count niches)
# Reads the precomputed niche assignments + per-bin RCTD weights from
# results/08-visium-hd-derived (env VISIUM_HD_DERIVED_DIR); no re-clustering here.
# ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(cli) })
proj_root <- path.expand("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
fig_dir   <- file.path(proj_root, "scripts/R/03-manuscript-figures")
source(file.path(proj_root, "scripts/R/00-0-source-functions.R"))    # .color_pal
source(file.path(fig_dir, "_fig1d-spatial-style.R")); source(file.path(fig_dir, "_mask-outlines.R"))
plot_dir  <- file.path(proj_root, "manuscript-figures/figure-7"); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
DERIVED   <- Sys.getenv("VISIUM_HD_DERIVED_DIR",
  unset = file.path(proj_root, "results/08-visium-hd-derived"))
SAMP_ORD  <- c("HB_hm01","TP_1083","TP_987","TR_03d_990")
.cond <- function(s) ifelse(grepl("^HB",s),"Healthy",ifelse(grepl("^TP",s),"Primary","Recurrent"))
META <- c("sample_id","condition","barcode","neopl_frac","reactive_astro_frac","in_tumor","x","y")

# niche -> descriptive name. Named by DOMINANT (raw mean-weight) component, NOT rare z-enriched
# programs: Neopl-ACR/ECM are rare in EVERY niche (max raw ~0.04/0.06), so a niche is never named
# after them. Class = summed malignant (Neopl-*) fraction: Tumor >0.5, Interface 0.15-0.5, Normal <0.15.
NICHE_NAME <- c(
  "1" ="Parenchyma",
  "2" ="Peritumoral parenchyma",
  "3" ="Tumor margin (neuronal / reactive-astro)",  # neuron-rich infiltrative edge (Neopl 0.25)
  "4" ="Tumor: Proneural (OPC)",                     # Neopl-OPC 0.62
  "5" ="Tumor: Bulk/OPC (reactive)",                 # Neopl-Bulk 0.26 + OPC 0.22 + AstroR 0.11
  "6" ="Tumor: OPC/Bulk",                            # Neopl-OPC 0.43 + Bulk 0.31
  "7" ="White matter",
  "8" ="Reactive-astro / immune",                    # AstroR 0.24 + myeloid/immune (neuron-poor)
  "9" ="Oligodendrocyte",                            # OPC-COP-OLG 0.72
  "10"="Homeostatic astrocyte (TE-NT)",              # Astrocyte TE-NT 0.28 (not vascular)
  "11"="Tumor: Bulk",                                # Neopl-Bulk 0.73
  "12"="Choroid/ventricle")
# legend / x-axis order: group by class -> Tumor (4), Interface (2), Normal (6)
NAME_LV <- unname(NICHE_NAME[as.character(c(4,6,5,11, 3,8, 2,1,7,9,10,12))])

nb <- read_csv(file.path(DERIVED, "bin-niche-16um.csv.gz"), show_col_types = FALSE)
W  <- bind_rows(lapply(SAMP_ORD, function(s) read_csv(file.path(DERIVED, sprintf("rctd-weights-%s-16um.csv.gz", s)), show_col_types = FALSE)))
CT <- setdiff(names(W), META)
d  <- left_join(nb, W[, c("sample_id","barcode",CT)], by = c("sample_id","barcode"))
d  <- d[!is.na(d[[CT[1]]]), ]
d$niche_name <- factor(NICHE_NAME[as.character(d$niche)], NAME_LV)
d$cond <- factor(.cond(d$sample_id), c("Healthy","Primary","Recurrent")); d$sample_id <- factor(d$sample_id, SAMP_ORD)
write_csv(tibble(niche = 1:12, name = NAME_LV), file.path(plot_dir, "Figure-7-niches-names.csv"))

# named 12-colour palette, keyed by name (order-independent). TUMOR niches = saturated warm
# (purple/magenta/red/orange); INTERFACE = distinct saturated (coral, teal) so the recurrent
# risers pop; NORMAL = pale/desaturated (spatial-plot-styling). No blue on a niche (the Tg
# tumor outline is dark navy).
NPAL <- c(
  "Tumor: Proneural (OPC)"                   = "#6A3D9A",  # purple   \
  "Tumor: OPC/Bulk"                          = "#E7298A",  # magenta   > TUMOR (saturated warm)
  "Tumor: Bulk/OPC (reactive)"               = "#E31A1C",  # red      /
  "Tumor: Bulk"                              = "#FF7F00",  # orange   /
  "Tumor margin (neuronal / reactive-astro)" = "#FB6A4A",  # coral     > INTERFACE (distinct sat.)
  "Reactive-astro / immune"                  = "#35978F",  # teal     /
  "Peritumoral parenchyma"                   = "#CBC4D9",  # pale purple \
  "Parenchyma"                               = "#DBDBDB",  # grey        \
  "White matter"                             = "#DDEED6",  # pale green   > NORMAL (paled)
  "Oligodendrocyte"                          = "#CFE8AC",  # pale olive  /
  "Homeostatic astrocyte (TE-NT)"            = "#C6DBEF",  # pale blue   /
  "Choroid/ventricle"                        = "#FED976")  # pale yellow/

# (1) named niche spatial map
outlines <- bind_rows(lapply(SAMP_ORD, mask_outlines)); outlines$sample_id <- factor(outlines$sample_id, SAMP_ORD)
p_map <- spatial_bin_grid(d, "niche_name", SAMP_ORD, outlines,
  title = "Composition niches (named; dark-blue outline = Tg tumor, grey = tissue, dashed = excluded)",
  # absent-level legend keys (tumor niches missing from the healthy first panel) are
  # coloured by the seed-layer fix inside spatial_bin_grid(); drop=FALSE keeps all 12.
  fill_scale = scale_fill_manual(values = NPAL, name = "niche", drop = FALSE))
ggsave(file.path(plot_dir, "Figure-7-niches-named-map.pdf"), p_map, width = 15, height = 5.0)

# (2) named enrichment heatmap
enr <- d |> group_by(niche_name) |> summarise(across(all_of(CT), mean), .groups="drop") |>
  pivot_longer(all_of(CT), names_to="cell_type", values_to="mw") |> group_by(cell_type) |> mutate(z=as.numeric(scale(mw))) |> ungroup()
enr$cell_type <- factor(enr$cell_type, rev(CT))
ggsave(file.path(plot_dir, "Figure-7-niches-celltype-heatmap.pdf"),
  ggplot(enr, aes(niche_name, cell_type, fill=z)) + geom_tile() +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", name="z") +
    labs(title="Cell-type enrichment per (named) niche", x=NULL, y=NULL) +
    theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=40, hjust=1)), width=9, height=6)
write_csv(enr, file.path(plot_dir, "Figure-7-niches-celltype-heatmap-sourcedata.csv"))

# (3) named niche frequency by condition
# tweak: complete all condition x niche combinations so every bar has equal width and
# zero-count niches are drawn at zero (position_dodge2 preserve="single").
fc <- d |> count(cond, niche_name) |> group_by(cond) |> mutate(frac=n/sum(n)) |> ungroup() |>
  complete(cond, niche_name, fill = list(n = 0, frac = 0))
ggsave(file.path(plot_dir, "Figure-7-niches-frequency-by-condition.pdf"),
  ggplot(fc, aes(niche_name, frac, fill=cond)) +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
    scale_fill_manual(values=c(Healthy=.color_pal[["sample_type"]][["Healthy"]], Primary=.color_pal[["sample_type"]][["Primary"]], Recurrent=.color_pal[["sample_type"]][["Recurrent"]]), name=NULL) +
    labs(title="Niche frequency by condition (named)", x=NULL, y="fraction of bins") +
    theme_bw(base_size=9) + theme(axis.text.x=element_text(angle=40, hjust=1)), width=11, height=5)
write_csv(fc, file.path(plot_dir, "Figure-7-niches-frequency-by-condition-sourcedata.csv"))

# (3b) — TUMOUR-restricted niche frequency: tg_expanded (Tg mask, matches the 16-6 cell-type
# distribution panel), Primary vs Recurrent. Restricting to tumour bins removes the normal-tissue dilution
# that masks the tumour-internal shift: the two reactive niches rise (Tumor margin ~0.07 -> 0.34;
# Reactive-astro/immune ~0.04 -> 0.28) while the malignant niches (Proneural/OPC-Bulk/Bulk) drop to
# ~0, i.e. the compact malignant niches give way to reactive/immune tissue at recurrence (the ACR
# signal redistributes into the rising reactive niches; reconciles with the cell-type ACR rise).
# Healthy has no tumour bins so it drops out. NB: Recurrent = n=1 (TR_03d_990, 3d post-RT, acute);
# niches at exactly 0 are an n=1 artefact -> descriptive only.
fc_tum <- d |> filter(tg_expanded, cond %in% c("Primary","Recurrent")) |> mutate(cond = droplevels(cond)) |>
  count(cond, niche_name) |> group_by(cond) |> mutate(frac=n/sum(n)) |> ungroup() |>
  complete(cond, niche_name, fill = list(n = 0, frac = 0))
ggsave(file.path(plot_dir, "Figure-7-niches-frequency-TUMOR-PTvTR.pdf"),
  ggplot(fc_tum, aes(niche_name, frac, fill=cond)) +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
    scale_fill_manual(values=c(Primary=.color_pal[["sample_type"]][["Primary"]], Recurrent=.color_pal[["sample_type"]][["Recurrent"]]), name=NULL) +
    labs(title="Niche frequency in tumour (Tg mask): Primary vs Recurrent", x=NULL, y="fraction of tumour bins") +
    theme_bw(base_size=9) + theme(axis.text.x=element_text(angle=40, hjust=1)), width=6.5, height=5)
write_csv(fc_tum, file.path(plot_dir, "Figure-7-niches-frequency-TUMOR-PTvTR-sourcedata.csv"))
cli::cli_alert_success("16-3 named niches (map + heatmap + frequency + tumour-restricted) -> {plot_dir}")
