library(Seurat)
library(magrittr)

# source the utility functions
source("scripts/run-seurat/utils.R", chdir = TRUE)

# import data
ptr <- csv2seurat(
    "outputs/seurat_integration/ortho_umi_count/tung_batch1_rs17.csv",
    "tung_batch1_rs17"
)
egr <- csv2seurat(
    "outputs/seurat_integration/ortho_umi_count/egr_nextseq_batch1_rs17.csv",
    "egr_nextseq_batch1_rs17"
)

# sctransform
ptr %<>% SCTransform(vst.flavor = "v2", verbose = TRUE) %>%
    RunPCA(npcs = 30, verbose = TRUE)

egr %<>% SCTransform(vst.flavor = "v2", verbose = TRUE) %>%
    RunPCA(npcs = 30, verbose = TRUE)

# integrate
object.list <- list(ptr = ptr, egr = egr)
features <- SelectIntegrationFeatures(
    object.list,
    nfeatures = 2000,
    verbose = TRUE
)
object.list <- PrepSCTIntegration(
    object.list = object.list,
    anchor.features = features,
    verbose = TRUE
)
anchors <- FindIntegrationAnchors(
    object.list = object.list,
    normalization.method = "SCT",
    anchor.features = features
)
combined.sct <- IntegrateData(
    anchorset = anchors,
    sample.tree = integration_order
    normalization.method = "SCT"
)
combined.sct <- RunPCA(
    combined.sct,
    verbose = TRUE
)

# write PC
write.csv(
    x = Embeddings(combined.sct, reduction = "pca")[, 1:30],
    file = "test/new_seurat/ptr-egr-30pc.csv",
    row.names = TRUE
    quote = FALSE
)

# write corrected counts
write.csv(
    x = GetAssayData(combined.sct, "SCT", "data") %>% t,
    file = "test/new_seurat/ptr-egr-corrected-data.csv",
    row.names = TRUE,
    quote = FALSE
)