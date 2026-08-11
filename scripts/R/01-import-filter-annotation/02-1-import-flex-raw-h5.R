# BiocManager::install("hdf5r")
setwd("~/Developer/git-pietras-lab/gbm-rcas-2025-lindgren-rosberg")

require(conflicted)
require(tidyverse)
require(Seurat)
require(Matrix)
require(qs2)
require(cli)
require(kableExtra)
require(readr)


# ---- define parameters ----
sample.metadata.raw.file = "./metadata/00-sample-metadata-raw-s16.csv"
sample.metadata.flex.file = "./metadata/00-sample-metadata-flex-s16.Rds"

# change this to your raw data path (h5 data available at ArrayExpress)
data.raw.path = "./data/processed/cellranger"


cli(cli_alert_info("importing sample metadata, {.var pdata} table, from file: {.file {sample.metadata.raw.file}}"))
pdata.raw <- readr::read_csv(file = sample.metadata.raw.file)
cli(cli_alert_info("Dimensions: {dim(pdata.raw)[1]} samples, {dim(pdata.raw)[2]} metadata variables"))
pdata.raw # A tibble: 16 × 9
kableExtra::kbl(table(pdata.raw$disease_state))
kableExtra::kbl(table(pdata.raw$sample_type))
kableExtra::kbl(table(pdata.raw$mouse_id))
# row.names(pdata) <- pdata$id

# read preproc flex metadata file
pdata <- readRDS(sample.metadata.flex.file)

# ---- Import Cellranger h5 data as Seurat Object----
# Read all cellranger h5 files for 16 smaples included in final data set
# - Load Cellranger h5's
# - Convert to Seurat object & Save
# - add sample metadata (pdata)


sd.list <- list() # list for storing each individual seurat object
h5.files <- list.files(
  path = data.raw.path,
  pattern="sample_filtered_feature_bc_matrix.h5",
  full.names = T, recursive = T, include.dirs = T
)


for (i in 1:nrow(pdata.raw)){ # i <- 1

  sample.i <- pdata.raw$sample_id[i]
  cellranger.i <- pdata.raw$cellranger_id[i]
  # h5.file <- pdata$h5.file[i]
  h5.file.i <- h5.files[grep(cellranger.i, h5.files)] # get this sample

  stopifnot(length(h5.file.i)==1)
  stopifnot(grep(cellranger.i, h5.file.i)==1)
  stopifnot(file.exists(h5.file.i))

  cli::cli_alert_info(" ... Reading h5, sample: {sample.i}, cellranger id: {cellranger.i},
    {.file {h5.file.i}}")

  h5 <- Seurat::Read10X_h5(
    filename=h5.file.i, use.names = TRUE,
    unique.features = TRUE
  )

  cli::cli_alert_info(" converting to Seurat object and appending list ")
  sd.list[[i]] <- suppressWarnings(
    Seurat::CreateSeuratObject(h5, project=sample.i)
  )

  cli::cli_alert_info(" adding metadata ")
  for(m in 1:ncol(pdata)){
    sd.list[[i]][[colnames(pdata)[m]]] <- pdata[i,m]
  }

  rm(h5, h5.file.i, sample.i)
  suppressMessages(gc())
}
names(sd.list) <- pdata$id
cli::cli_alert_success("Ok, loaded all sample h5 {names(sd.list)}")
gc()
# check that all feature dims are same
lapply(sd.list, dim)
stopifnot(length(unique(unlist(lapply(sd.list, function(x) dim(x)[1]))))==1)




# ---- Merge Samples -----
# Merge samples into one Seurat object `sd.layered`
cli::cli_alert_info("Merging {length(sd.list)} seurat objects in {.var sd.list} into one oject")

x <- unlist(lapply( 1:length(sd.list), function(x) {paste0("sd.list[[",x,"]]")}))
y <- paste0(x[-1], collapse = ", ")
merge.string <- paste0("merge(", x[1], ", c(" ,y , "), add.cell.ids = pdata$id)", collapse=', ')
sd.layered <- eval(parse(text=merge.string))
cli::cli_alert_info(paste0(
  "  Dim: ", cli::col_br_yellow(dim(sd.layered))
))
cli::cli_alert_success("Done! {.var sd.layered} generated.")

cli::cli_alert("saving data to {.file {params$seurat.save.file }}")

seurat.save.file =  "./data/processed/seurat/02-1-flex-raw-s16-c92737.qs"
system.time(
  qs2::qs_save(sd.layered, seurat.save.file)
)




