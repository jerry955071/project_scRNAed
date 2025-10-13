library(Seurat)
library(magrittr)
library(R.utils)
library(future)
source("scripts/run-seurat/utils.R", chdir = TRUE)

# # parse arguments
parse_arguments()

# # set up future
# maxSize_gb <- 50
# options(future.globals.maxSize= maxSize_gb*1000*1024^2)
# plan("multicore", workers = as.integer(threads))
# plan()


# create output directory
if (!dir.exists(dirname(output_pca_projection_csv))) {
    dir.create(dirname(output_pca_projection_csv), recursive = TRUE)
}

# load data
object.list <- sapply(
    samples2integrate %>% strsplit(",") %>% unlist(),
    function(x) {
        csv2seurat(
            csv_file = x,
            project_name = basename(x) %>% tools::file_path_sans_ext()
        )
    }
)
sample.tree <- build_seurat_sample_tree(integration_order)
print(sample.tree)

# qc data
qc_results <- sapply(
    object.list,
    function(x) {
        qc_seurat_object(
            seurat_object = x,
            n_feature_cutoff_on_sample = as.integer(n_feature_cutoff_on_sample),
            n_barcode_cutoff_on_sample = as.integer(n_barcode_cutoff_on_sample)
        )
    }
)
if (any(!qc_results)) {
    stop("QC failed")
}

# run "custom integration pipeline 2"
combined <- custom_integration_pipeline_2(
    object.list = object.list,
    sample.tree = sample.tree,
    n_pcs = as.integer(n_pcs),
    nfeatures = as.integer(nfeatures)
)
# message("Seurat result: ", str(combined))


# write PC
write.csv(
    x = Embeddings(combined, reduction = "pca"),
    file = output_pca_projection_csv,
    row.names = TRUE,
    quote = FALSE
)

# write to RDS
saveRDS(combined, output_seurat_combined_object_rds)
