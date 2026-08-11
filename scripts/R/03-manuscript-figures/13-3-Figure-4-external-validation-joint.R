# 13-3 — External validation: our neoplastic programs vs reference GBM states (reviewers-only)
# Home = Figure S4.
#
# FIGURE LAYER (runs now, no raw data needed): reads the portable correlation CSVs in
# results/31-external-validation/ (tracked) and renders the manuscript PDFs + source data.
# Re-scoring from raw is 09-1-external-validation-score.R (needs data/external/, HPC for Nomura).
#
# NOTE ON NAMING: the reference panel uses the SHORT Richards dev/injury (top-250) and the
# COALESCED Nomura 9 states (not the 15 raw MPs) — methodologically the manuscript choice.
# Display labels ("Richards Developmental", "Nomura AC-like") are kept concise here
# (not harmonised to the Richards_*_2021 display scheme used in Fig 4).

suppressMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2)
                  library(ComplexHeatmap); library(circlize); library(grid)})
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")
RES <- "results/31-external-validation"          # portable corr CSVs + sourcedata (tracked)
O   <- "manuscript-figures/figure-4"             # figure PDFs (tracked); provisional S4
dir.create(O, showWarnings = FALSE, recursive = TRUE)

# ---- MAIN: 6 neoplastic programs x 20 reference states, Spearman across malignant cells ----
datasets  <- c("Soni","GBmap","Sussman","Suter","Nomura")
our_progs <- c("Neopl-ACR (ours)","Neopl-ECM (ours)","Neopl-COP (ours)",
               "Neopl-NC (ours)","Neopl-OPC (ours)","Neopl-CC/Bulk (ours)")
ref_order <- c("Neftel AC","Neftel MES1","Neftel MES2","Neftel NPC1","Neftel NPC2","Neftel OPC",
  "Nomura AC-like","Nomura MES-like","Nomura Hypoxia","Nomura GPC-like","Nomura NPC-like",
  "Nomura OPC-like","Nomura NEU-like","Nomura Stress","Nomura Cilia-like",
  "Richards Developmental","Richards Injury-Response")

rows <- list()
for (d in datasets) {
  f <- file.path(RES, paste0("31-corr-", d, "-neoplastic.csv"))
  if (!file.exists(f)) { message("missing ", f, " - skipping"); next }
  cm <- read_csv(f, show_col_types = FALSE)
  for (op in intersect(our_progs, cm$signature)) {
    r <- cm[cm$signature == op, , drop = FALSE]
    for (rf in intersect(ref_order, colnames(cm)))
      rows[[length(rows)+1]] <- tibble(dataset = d, our = op, ref = rf, corr = r[[rf]][1])
  }
}
df <- bind_rows(rows)
df$ref     <- factor(df$ref, levels = rev(ref_order[ref_order %in% df$ref]))
df$dataset <- factor(df$dataset, levels = datasets[datasets %in% df$dataset])
df$our     <- factor(df$our, levels = our_progs)

# display labels (Nomura data key -> GBM-CARE) + a per-dataset colour key for the highlight strip
ds_show <- c(Soni = "Soni", GBmap = "GBmap", Sussman = "Sussman", Suter = "Suter", Nomura = "GBM-CARE")
ds_cols <- c(Soni = "#4E79A7", GBmap = "#59A14F", Sussman = "#B07AA1", Suter = "#9C755F", Nomura = "#E15759")
dl   <- levels(df$dataset)          # datasets present, in column order (data keys)
xcol <- unname(ds_cols[dl])         # colour per column, matched to the axis breaks
# colour bar just below the tiles: one scalar-fill block per dataset (kept off the corr scale)
strip_layers <- lapply(seq_along(dl), function(k)
  annotate("rect", xmin = k - 0.5, xmax = k + 0.5, ymin = -0.05, ymax = 0.42, fill = xcol[k]))

g <- ggplot(df, aes(dataset, ref, fill = corr)) +
  geom_tile(colour = "grey92") +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 2.3) +
  strip_layers +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1,1)) +
  scale_x_discrete(labels = ds_show) +
  scale_y_discrete(expand = expansion(add = c(1.1, 0.6))) +
  facet_wrap(~our, ncol = 3) +
  labs(title = "External validation - our neoplastic programs vs reference states",
       subtitle = "Spearman correlation across malignant cells, per dataset.",
       x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = xcol, face = "bold"),
        panel.grid = element_blank(), plot.margin = margin(6, 6, 12, 6))
ggsave(file.path(O, "Figure-S4-extval-neoplastic-programs-vs-states.pdf"), g, width = 12, height = 9)
write_csv(df, file.path(O, "Figure-S4-extval-neoplastic-programs-vs-states-sourcedata.csv"))
cat("saved main external-validation figure (", nlevels(df$dataset), "datasets )\n")

# [reviewer-only supplements (all-cells correlation heatmaps, Soni-by-model) removed from publication repo]
cat("done\n")
