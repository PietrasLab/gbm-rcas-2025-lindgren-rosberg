# 15-3 — Supp Fig 6 (Figure 6 supplement): TNF / COMPLEMENT / IL1 sender->receiver
# ============================================================================
# Alignment with the Bhat/Schmitt/Richards + Liddelow myeloid->astrocyte
# injury-response axis: do TNFa / C1q / IL1A pathways appear, and who signals
# to whom? Also supports the myeloid->astrocyte link and the astrocyte- vs
# myeloid-derived signalling comparison.
# Single panel = the sender->receiver bubble.
#
# Finding:
#   TNF, COMPLEMENT and IL1 ARE present in the Level_3ACM SecretedSignaling
#   network. IL1A/IL1B from myeloid cells (Microglia, Macrophage, DC, Neutrophil)
#   reach several stromal/reactive targets incl. Astrocyte R, consistent with the
#   Liddelow/Richards axis, extended by additional TWEAK/MIF/GRN myeloid signalling.
#
# Panel: sender (source) -> receiver (target) bubble, faceted pathway x condition;
#        dot size = summed CellChat communication probability.
# ============================================================================

suppressMessages({
  library(CellChat); library(dplyr); library(ggplot2); library(readr)
})
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

cc.dir   <- "./results/07-3-interactions"
plot.dir <- "./manuscript-figures/figure-6"
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)

FOCUS <- c("TNF", "COMPLEMENT", "IL1")

# ---- pull TNF/COMPLEMENT/IL1 interactions from each condition object ---------
extract_pw <- function(rds, condition) {
  cc <- readRDS(file.path(cc.dir, rds))
  present <- unique(c(FOCUS[FOCUS %in% cc@netP$pathways],
                      cc@netP$pathways[grepl("^TNF$|^IL1$|^COMPLEMENT$", cc@netP$pathways)]))
  if (length(present) == 0) return(NULL)
  df <- subsetCommunication(cc, signaling = present)
  df$condition <- condition
  df
}

message("Loading Primary/Recurrent Level_3ACM SecretedSignaling objects...")
df <- bind_rows(
  extract_pw("07-3-CellChat_object_Level_3ACM_PrimaryTumors_SecretedSignaling.rds",   "Primary"),
  extract_pw("07-3-CellChat_object_Level_3ACM_RecurrentTumors_SecretedSignaling.rds", "Recurrent")
) %>%
  filter(pathway_name %in% FOCUS) %>%
  mutate(pathway_name = factor(pathway_name, levels = FOCUS),
         condition    = factor(condition, levels = c("Primary", "Recurrent")))

# ---- aggregate to summed sender->receiver strength (the plotted table) -------
agg <- df %>%
  group_by(pathway_name, source, target, condition) %>%
  summarise(strength = sum(prob), n = n(), .groups = "drop")

# ---- panel: sender -> receiver bubble ---------------------------------------
p <- ggplot(agg, aes(source, target, size = strength, colour = pathway_name)) +
  geom_point(alpha = 0.85) +
  facet_grid(pathway_name ~ condition) +
  scale_size_area(max_size = 8, name = "summed\nstrength") +
  scale_colour_manual(values = c(TNF = "#D55E00", COMPLEMENT = "#0072B2", IL1 = "#009E73"),
                      guide = "none") +
  labs(title = "TNF / COMPLEMENT / IL1: sender -> receiver (Level_3ACM SecretedSignaling)",
       subtitle = "Dot = summed CellChat probability. IL1A/IL1B from myeloid cells reach Astrocyte R (Liddelow/Richards axis) among other stromal targets.",
       x = "sender (source)", y = "receiver (target)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot.dir, "Figure-S6-cellchat-TNF-COMP-IL1-sender-receiver.pdf"),
       p, width = 11, height = 9)

# ---- source data (plotted points) -------------------------------------------
write_csv(agg, file.path(plot.dir,
          "Figure-S6-cellchat-TNF-COMP-IL1-sender-receiver-sourcedata.csv"))

message("Done. Figure + source data written to ", plot.dir)
