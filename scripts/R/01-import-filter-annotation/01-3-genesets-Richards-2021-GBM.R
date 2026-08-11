
# Signatures collected and presented in Richards et al 2021,
# A number of Glioblastoma gene signuatures (including definition of a stem cells GSC signature)
# https://www.nature.com/articles/s43018-020-00154-9
# Gradient of Developmental and Injury Response transcriptional states defines functional vulnerabilities underpinning glioblastoma heterogeneity
# Nat Cancer. 2021 Feb;2(2):157-173. doi: 10.1038/s43018-020-00154-9. Epub 2021 Jan 4.
# Richards et al

# From supplementary. Load, and curate the gene lists. Translate to Mouse homologues
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

richards.gsc.sigs <- "./data/genesets/RichardsWhitley_SuppTable_7.csv"
richards.gbm.sigs <- "./data/genesets/RichardsWhitley_SuppTable_6.csv"


# Read CSV file
library(readr)
library(tidyr)
library(dplyr)

# Read CSV file
gene_sets <- readr::read_csv(file.path(richards.dir, richards.gbm.sigs))

# Convert wide format to long format (if needed)
gene_sets_long <- gene_sets %>%
  pivot_longer(cols = everything(), names_to = "geneset", values_to = "gene") %>%
  filter(!is.na(gene))  # Remove empty entries

# convert to list
gene_sets_list <- gene_sets %>%
  pivot_longer(cols = -1, names_to = "geneset", values_to = "gene") %>%  # Convert wide to long format
  filter(!is.na(gene)) %>%  # Remove NA values
  group_by(geneset) %>%
  summarise(genes = list(gene), .groups = "drop") %>%  # Convert to list format
  tibble::deframe()  # Convert to a named list

# Check the output
str(gene_sets_list)

# gsc signatures
gsc.df <- readr::read_csv(file.path(richards.dir, richards.gsc.sigs))
head(gsc.df)
table(gsc.df$Group)
# group on Group, filter on logg2FoldCahnge > 1
my.df <- gsc.df %>%
  dplyr::filter(abs(log2FoldChange) > 2) %>%
  dplyr::filter(padj < 0.001) %>%
  dplyr::arrange(desc(abs(log2FoldChange))) %>%
  group_by(Group)
# get top N for each group
topN <- 250
gsc.list <- my.df %>%
  slice_head(n = topN) %>%
  summarise(genes = list(Gene), .groups = "drop") %>%
  tibble::deframe()
str(gsc.list)



# append to the gene_sets_list
names(gsc.list) <- paste0(names(gsc.list),"_Richards")
gene_sets_list <- c(gene_sets_list, gsc.list)

# save as rds
saveRDS(gene_sets_list, file.path("./references/genesets/Richards_NatCancer_2021_GeneSets.rds"))

# Convert Entrez IDs to Mouse Gene Symbols
library(org.Hs.eg.db)
library(AnnotationDbi)
library(purrr)
library(homologene)

mouse_symbols <- map(gene_sets_list, ~ {
  homologs <- human2mouse(.x)  # 9606 = Human, 10090 = Mouse
  return(homologs$mouseGene)  # Extract mouse gene symbols
})

saveRDS(mouse_symbols, file.path("./data/genesets/Richards_NatCancer_2021_GeneSets_Mouse.rds"))
