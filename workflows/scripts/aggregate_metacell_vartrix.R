# Redirect stdout/stderr to log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

library(Seurat)
library(magrittr)
library(Matrix)

# get variables
path_mc_label <- snakemake@input[["mc_label"]]
path_alt_mtx <- snakemake@input[["alt_mtx"]]
path_ref_mtx <- snakemake@input[["ref_mtx"]]
path_row_idx <- snakemake@input[["row_idx"]]
path_col_idx <- snakemake@input[["col_idx"]]
path_out_alt <- snakemake@output[["alt_mtx"]]
path_out_ref <- snakemake@output[["ref_mtx"]]
path_out_row <- snakemake@output[["row_idx"]]
path_out_col <- snakemake@output[["col_idx"]]

# # for debugging
# path_mc_label <- "outputs/mcRigor/{sample}/metacell.csv" %>% glue::glue(sample="ptr_tenx_tsv2")
# path_alt_mtx <- "outputs/VariantCalling/vartrix/{sample}/alt.mtx" %>% glue::glue(sample="ptr_tenx_tsv2")
# path_ref_mtx <- "outputs/VariantCalling/vartrix/{sample}/ref.mtx" %>% glue::glue(sample="ptr_tenx_tsv2")
# path_row_idx <- "outputs/VariantCalling/vawk/{sample}.snv.loci.txt" %>% glue::glue(sample="ptr_tenx_tsv2")
# path_col_idx <- "outputs/Remapping/renamer/{sample}/barcodes.tsv" %>% glue::glue(sample="ptr_tenx_tsv2")

# read files
mc_label <- read.csv(path_mc_label, row.names = 1)
mc_label[["Metacell"]] <- paste0("mc-", mc_label[["Metacell"]])
alt_mtx <- ReadMtx(
  mtx = path_alt_mtx,
  cells = path_col_idx,
  cell.column = 1,
  features = path_row_idx,
  feature.column = 1
) %>% CreateSeuratObject # underscores '_' in Scaffold names were replaced by dashes '-'
ref_mtx <- ReadMtx(
  mtx = path_ref_mtx,
  cells = path_col_idx,
  cell.column = 1,
  features = path_row_idx,
  feature.column = 1
) %>% CreateSeuratObject # underscores '_' in Scaffold names were replaced by dashes '-'

# aggregate alt/ref counts per MetaCell
col_mc <- colnames(mc_label)
alt_mtx@meta.data[col_mc] <- mc_label[Cells(alt_mtx), col_mc]
ref_mtx@meta.data[col_mc] <- mc_label[Cells(ref_mtx), col_mc]
alt_mtx <- subset(alt_mtx, mcRigor_sc == "trustworthy")
ref_mtx <- subset(ref_mtx, mcRigor_sc == "trustworthy")
mc_alt_mtx <- AggregateExpression(alt_mtx, assays="RNA", group.by = "Metacell")$RNA
mc_ref_mtx <- AggregateExpression(ref_mtx, assays="RNA", group.by = "Metacell")$RNA

# write aggregated alt/ref mtx to file
writeMM(mc_alt_mtx, file = path_out_alt)
writeMM(mc_ref_mtx, file = path_out_ref)
write.table(
  data.frame(dimname_1 = mc_alt_mtx@Dimnames[[1]]),
  path_out_row,
  col.names = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(dimname_1 = mc_alt_mtx@Dimnames[[2]]),
  path_out_col,
  col.names = FALSE,
  row.names = FALSE
)

# Cleanup connections
sink(type = "message")
sink()
