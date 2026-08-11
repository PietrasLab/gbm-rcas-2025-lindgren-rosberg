
## Look at raw counts values of the different Tg-s.
## rcas_both_3bin: A binned proxy, defined on raw count values of the Tg-s
## Objective - classify cells into bins wether positive / ambiguous / negative for transgene expression

# ---- Env ----
require(conflicted)
require(tidyverse)
require(qs2)
require(cli)
require(glue)
suppressMessages(require(Seurat))
suppressMessages(require(Matrix))
suppressMessages(require(gridExtra))
suppressMessages(require(ggplot2))
# library(scCustomize)
# require(peRcebe)
require(CellMetaVerse)
require(kableExtra)
require(scCustomize)
require(scDblFinder)
require(BiocParallel)

set.seed(169)

# ---- Set path ----
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

# ---- Source functions ----
source("./scripts/R/00-0-source-functions.R")

# ---- Load data -----
cli::cli_alert_info("Loading Seurat oject (non downsampled) ")
sd.merged <- qs2::qs_read("./data/processed/seurat/02-1-flex-raw-s16-c92737.qs") %>%
  SeuratObject::JoinLayers()
cell_barcodes <- paste0(sd.merged@meta.data$sample_id, "_", colnames(sd.merged))
# removing trailing _digit AT the END
cell_barcodes <- gsub("_\\d+$", "", cell_barcodes)
sd.merged@meta.data$cell_barcode <- cell_barcodes
colnames(sd.merged) <- cell_barcodes
cli::cli_alert_info("Seurat object contains {ncol(sd.merged)} cells and {nrow(sd.merged)} features.")

# ---- Load sample and probe feature metadata ----
pdata <- readRDS("./metadata/00-sample-metadata-flex-s16.Rds")
cli::cli_alert_info("Loading FlexProbes feature metadata.")
fdata <- read.delim("./metadata/00-probeset-flex.csv", header = T, sep=",")
str(fdata)
head(fdata)



# ---- Call Raw Counts for Transgene custom probes ----
## call different variants of rcas cell calling (non normalized)
# Bins and Individual Probe Calls
custom.probes <- fdata[grep("Tg-", fdata$external_gene_name),]$external_gene_name

custom.probes
custom.probes.bin <- gsub("Tg-","",paste0(custom.probes,"_bin"))
custom.probes.bin <- gsub("-","_",custom.probes.bin)
custom.probes.call <- gsub("Tg-","",paste0(custom.probes,"_call"))
custom.probes.call <- gsub("-","_",custom.probes.call)

probeCalls <- data.frame(
  cell_name=names(sd.merged$orig.ident),
  orig_ident = as.factor(sd.merged$orig.ident)
)
rownames(probeCalls) <- NULL
str(probeCalls)


for( i in 1:length(custom.probes)){
  i.name <- custom.probes[i]
  b.name <- custom.probes.bin[i]
  c.name <- custom.probes.call[i]

  my.breaks <-  c(0, 1, 2, 6, Inf)
  my.labels <- paste(
    "count",
    as.character(my.breaks)[-length(my.breaks)],
    as.character(my.breaks)[-1],
    sep="-")

  probeCalls[,b.name] <- cut(
    # sd.merged@assays$RNA@counts[i.name, ],
    LayerData(sd.merged, assay = "RNA", layer = "counts")[i.name,],
    breaks = c(0, 1, 2, 6, Inf),
    labels = c("count-00","count-01","count-02-05", "count-06+"),
    include.lowest = T, right = F)


  p <- peRcebe::barplot_stacked(
    probeCalls,
    group.var = "orig_ident",
    color.var = b.name,
    scaled.y = F)
  plot(p)
}
colnames(probeCalls)
table(probeCalls[probeCalls$orig_ident=="TP_01b","hPDGFB_2_bin"], probeCalls[probeCalls$orig_ident=="TP_01b","hPDGFB_nHA_bin"])


## ---- Create Joint Bins for the hPDGFB probes (fixed) ----
# Start from the nHA bin
table(probeCalls$hPDGFB_nHA_bin)
table(probeCalls$hPDGFB_2_call)
table(probeCalls$hPDGFB_2_bin)

probeCalls$hPDGFB_bin <- probeCalls$hPDGFB_nHA_bin
probeCalls$hPDGFB_bin[
  probeCalls$hPDGFB_nHA=="count-00" &
    probeCalls$hPDGFB_2_call!="count-00"
] <- probeCalls$hPDGFB_2_call[
  probeCalls$hPDGFB_nHA=="count-00" &
    probeCalls$hPDGFB_2_call!="count-00"
]
table(probeCalls$hPDGFB_bin)


probeCalls$RFP_bin <- cut((
  LayerData(sd.merged, assay = "RNA", layer = "counts")["Tg-RFP-1", ] +
    LayerData(sd.merged, assay = "RNA", layer = "counts")["Tg-RFP-2", ]
),
  breaks = c(0, 1, 2, 6, Inf),
  labels = c("count-00","count-01","count-02-05", "count-06+"),
  include.lowest = T, right = F)


table(probeCalls$RFP_bin)
table(probeCalls$RFP_1_bin, probeCalls$RFP_bin)
table(probeCalls$RFP_2_bin, probeCalls$RFP_bin)

table(probeCalls$RFP_bin, probeCalls$hPDGFB_bin)
# rcas_any = rcas_call
# rcas_both   1+ in both
# rcas_strict 2+ in both
probeCalls$rcas_bin <- probeCalls$hPDGFB_bin
probeCalls$rcas_bin[probeCalls$RFP_bin!="count-00" & probeCalls$hPDGFB_bin=="count-00"] <-
  probeCalls$RFP_bin[probeCalls$RFP_bin!="count-00" & probeCalls$hPDGFB_bin=="count-00"]
table(probeCalls$rcas_bin)

probeCalls$rcas_call <- "rcas_neg"
probeCalls$rcas_call[probeCalls$rcas_bin!="count-00" ] <- "rcas_pos"
table(probeCalls$rcas_call)

# both call: require only 1+ in expression
probeCalls$rcas_both <- "rcas_neg"
probeCalls$rcas_both[probeCalls$hPDGFB_bin!="count-00"] <- "hpdgfb_pos"
probeCalls$rcas_both[probeCalls$RFP_bin!="count-00"] <- "rfp_pos"
probeCalls$rcas_both[probeCalls$hPDGFB_bin!="count-00" & probeCalls$RFP_bin!="count-00"] <- "rcas_pos"
table(probeCalls$rcas_both)

# a slightly stricter one (require 2+ in expression)
probeCalls$rcas_both_stricter <- "rcas_neg"
probeCalls$rcas_both_stricter[ !probeCalls$hPDGFB_bin %in% c("count-00","count-01")] <- "hpdgfb_pos"
probeCalls$rcas_both_stricter[ !probeCalls$RFP_bin  %in% c("count-00","count-01") ] <- "rfp_pos"
probeCalls$rcas_both_stricter[
  !probeCalls$hPDGFB_bin %in% c("count-00","count-01") &
    !probeCalls$RFP_bin %in% c("count-00","count-01")
] <- "rcas_pos"
table(probeCalls$rcas_both_stricter)
table(probeCalls$rcas_both_stricter, probeCalls$rcas_both)
table(probeCalls$rcas_both_stricter, probeCalls$rcas_call)

probeCalls$rcas_strict <- "rcas_neg"
probeCalls$rcas_strict[probeCalls$hPDGFB_bin %in% c("count-02-05","count-06+") & probeCalls$RFP_bin %in% c("count-02-05","count-06+") ] <- "rcas_strict"
table(probeCalls$rcas_strict)

##. rcas_both_3bin - separate into 3 bins based
probeCalls <- probeCalls %>%
  mutate(rcas_both_3bin = "rcas_ambiguous") %>%
  mutate(rcas_both_3bin = if_else(rcas_both_stricter=="rcas_pos", "rcas_pos", rcas_both_3bin)) %>%
  mutate(rcas_both_3bin = if_else(rcas_both=="rcas_neg", "rcas_neg", rcas_both_3bin))
table(probeCalls$rcas_both_3bin)

table(probeCalls$rcas_both_3bin, probeCalls$rcas_both_stricter)
table(probeCalls$rcas_both_3bin, probeCalls$rcas_both)
table(probeCalls$rcas_both_stricter, probeCalls$rcas_both)
str(probeCalls$rcas_both_stricter)
colnames(probeCalls)


# > table(probeCalls$rcas_both_3bin)
#
# rcas_ambiguous       rcas_neg       rcas_pos
# 17700          37964          37073

# save dataframe with transgene probeCalls rcas_both_3bin
df <- probeCalls %>% dplyr::select(cell_name, rcas_both_3bin)
readr::write_csv(df, "./results/03-qc/03-5-rcas-transGenes-probeCalls_rcas_both_3bin.csv")
readr::write_csv(probeCalls, "./results/03-qc/03-4-rcas-transGenes-probeCalls_full_table.csv")





