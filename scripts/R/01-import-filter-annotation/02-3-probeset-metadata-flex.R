# ---- get feature/probe data from 10x CSV ----
# Get feature metadata

# Format the probe csv file to read into R

# Get gene metadata from embl datbase.Use the probe csv to get
# correct ENSG identifiers. The gene symbols from 10x gene tables do not entirely
# match ensembl (possibly different versions). The ENSG ids should be more stabile
# * Get gene name : ensg id probe match from csv table
# * get ensembl feature data using biomart (chromosome etc)
# * use probe table ensembl id to get fe
# * keep the ensemble gene symbol in fdata, but continue use the 10x gene symbol
# * match fdata table to seurat gene symbols and add full metadata to object
#
# Note that the seurat object carries 19073 rows (not 19475) as seen in raw probe csv file
#
# Note that when importing the seurat object, gene symbols are used to generate seurat feature ids. Some gene symbols are duplicated (i.e. a few probes with separate ensembl ids are mapped to the same gene symbol). When these are imported into seurat the gene symbol for duplicates are appended with a .1 .2
# etc

# >**_NOTE_**
# Mind that a few duplicated gene symbols exits

require(dplyr)
require(tidyr)
require(stringr)
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

probe.table <- read.delim(
  file.path("references/cellranger/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A-flex-rcas-v3-split.csv"),
  header = F, sep='|')
str(probe.table)
tail(probe.table)
colnames(probe.table) <- c(
  "sequence","ensembl_gene_id","external_gene_name","probe_hash"
)
str(probe.table) # 55552
str(unique(probe.table$ensembl_gene_id)) # 19478

# Access feature-level meta-data
# initialize connection to mart, may take some time if the sites are
# unresponsive.
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
# fetch chromosome info plus some other annotations
fdata.ensembl <- try(biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id", "external_gene_name",
    "description", "gene_biotype", "chromosome_name", "start_position"
  ),
  mart = mart, useCache = F)
)

# rename gene symbol from ensembl (keep it but dont use it)
fdata.ensembl <- fdata.ensembl %>%
  dplyr::rename(external_gene_name_ensembl = external_gene_name)

# Merge the two tables
probe.table <- probe.table %>%
  dplyr::select(2:4,1) %>%
  left_join(y = fdata.ensembl, by = "ensembl_gene_id")

str(probe.table) # 55546
str(unique(probe.table$ensembl_gene_id)) # 19478, OK! all probes kept

# save probe table with proper annotations
probe.table.file <-  "./metadata/00-probeset-flex.csv"
write.table(
  probe.table,
  probe.table.file,
  col.names = T, row.names = F, sep=","
)

# read proble table
probe.table <- read.delim(
  probe.table.file,
  header = T, sep=","
)
# str(probe.table[grepl("DEPRECATED",probe.table$ensembl_gene_id),]) # 59 deprecated
# str(probe.table[is.na(probe.table$external_gene_name_ensembl),]) # 9 custom
# str(probe.table[is.na(probe.table$external_gene_name),]) # 0 ok!
str(probe.table) # 55552
probe.table.filtered <- probe.table %>%
  dplyr::filter(!str_detect(ensembl_gene_id, "^DEPRECATED")) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)
str(probe.table.filtered) # 19419

# a few gene symbols map to multiple ENMUSG ids, e.g. Aldoa ENSMUSG00000030695
# i.e. rows have duplicated external_gene_name but unique ensembl id
probe.table.filtered %>%
  dplyr::filter(
    duplicated(external_gene_name) |
      duplicated(external_gene_name, fromLast = TRUE))

probe.table.filtered[
  probe.table.filtered$external_gene_name=="Aldoa",]
fdata.ensembl[fdata.ensembl$external_gene_name_ensembl=="Aldoa", ]
# WHERE ARE THESE in the h5 or sdata objects (since this uses gene symbol)


# create an empty feature metadata frame for seurat object
# based on its gene symbol rownames. match from proble table
fdata.sds <- sd.list[[i]][["RNA"]][[]]
fdata.sds$external_gene_name <- rownames(sd.list[[i]])
str(fdata.sds) # 19073
# fdata.sds[grepl(".", rownames(fdata.sds)),]
# probe.table[probe.table$external_gene_name=="Jakmip1",]
# 17 entries in the Seurat object stem from duplicated gene symbols
#
dup.symbols <- fdata.sds %>%
  dplyr::filter(str_detect(external_gene_name, "\\."))
dup.symbols <- gsub("\\.[0-9]$", "", dup.symbols$external_gene_name)
fdata.sds[grepl("H2-M10",fdata.sds$external_gene_name), ]
# unfortunately, the order of entries are not neccesarily sequential

# for these duplicated gene symbols we cannot assign proper ensembl ids
# for now, set all these to NA with regards to metadata
# performed by removing these genes from the ensemble/probe table
# also remove symbols that are duplicated in the ensembl table
str(probe.table.filtered) # 19416
str(fdata.sds) # 19067
probe.table.filtered <- probe.table.filtered %>%
  dplyr::filter( ! ( external_gene_name %in% dup.symbols) ) %>%
  dplyr::filter(
    ! ( duplicated(external_gene_name) )
  )
str(probe.table.filtered) # 19411
fdata <- fdata.sds %>%
  left_join(y=probe.table.filtered, by="external_gene_name")

str(fdata) # 19073 OK!
fdata.file <-  "./metadata00-probeset-flex.csv"
write.table(
  fdata, file = fdata.file,
  sep=',', row.names = F, col.names = T
)

