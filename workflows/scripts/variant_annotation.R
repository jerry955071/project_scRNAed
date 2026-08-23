# Redirect stdout/stderr to Snakemake log file
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

# 1. Setup
# Load required packages
library(VariantAnnotation)
library(GenomicFeatures)
library(BSgenome)
library(Biostrings)

# get snakemake variables
path_ref <- snakemake@input[["path_ref"]]
path_gff <- snakemake@input[["path_gff"]]
path_hom <- snakemake@input[["path_hom"]]
path_rna_variants <- snakemake@input[["path_rna_variants"]]

path_rna_editing_vcf <- snakemake@params[["path_rna_editing_vcf"]] # This is read from param
path_variant_location <- snakemake@output[["path_variant_location"]]
path_variant_annotation <- snakemake@output[["path_variant_annotation"]]

param_genome <- snakemake@params[["genome"]] # Ptr


# 2. Load annotation and genome fasta
ref <- FaFile(path_ref)
txdb <- makeTxDbFromGFF(path_gff)


# 3. Load variants and predict coding changes
# Homozygous loci
cat("Loading homozygous loci from: ", path_hom, "\n")
hom_ref <- readVcf(path_hom)

# RNA variants
cat("Loading RNA variant loci from: ", path_rna_variants, "\n")
vcf <- readVcf(path_rna_variants, genome = param_genome)

# Subset RNA variants to homozgous loci (putatively RNA editing)
vcf_hom <- subsetByOverlaps(vcf, hom_ref)
writeVcf(
    vcf_hom,
    path_rna_editing_vcf,
    index = TRUE
)

# predict consequence of RNA editing
variants <- locateVariants(vcf_hom, txdb, AllVariants())
coding <- predictCoding(vcf_hom, txdb, ref)

# write var_loc to file
var_loc <- data.frame(
    "CHROM" = as.character(variants@seqnames),
    "POS" = variants@ranges@start
)
var_loc$LOCATION <- mcols(variants)$LOCATION
var_loc$LOCSTART <- mcols(variants)$LOCSTART
var_loc$GENEID <- mcols(variants)$GENEID
write.table(
    var_loc,
    path_variant_location,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# write consequence to file 
out_data <- data.frame(
    "CHROM" = as.character(coding@seqnames),
    "POS" = coding@ranges@start
)
out_data$REF <- as.character(mcols(coding)$REF); dim(out_data)
ALT <- sapply(mcols(coding)$ALT, paste0, collapse = ",")
out_data$ALT <- as.character(ALT); dim(out_data)
out_data$GENEID <- as.character(mcols(coding)$GENEID); dim(out_data)
PLOC <- sapply(mcols(coding)$PROTEINLOC, paste0, collapse = ",")
PLOC <- as.character(PLOC)
out_data$PROTEINLOC <- PLOC; dim(out_data)
out_data$CONSEQUENCE <- as.character(mcols(coding)$CONSEQUENCE); dim(out_data)
out_data$REFCODON <- as.character(mcols(coding)$REFCODON); dim(out_data)
out_data$VARCODON <- as.character(mcols(coding)$VARCODON); dim(out_data)
out_data$REFAA <- as.character(mcols(coding)$REFAA); dim(out_data)
out_data$VARAA <- as.character(mcols(coding)$VARAA); dim(out_data)

write.table(
    out_data,
    path_variant_annotation,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# Cleanup Snakemake log file connections
sink(type = "message")
sink()

# # 4. Subset transcripts overlap with variants
# # Get CDS or exons for transcripts
# # cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
# # exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
# # write.csv(
# #     as.data.frame(exons),
# #     "tests/test-VariantAnnotation2/exons_gtf.csv",
# #     row.names = FALSE,
# #     quote = FALSE
# # )

# # Get transcripts with variations on exon
# # cds_nonsynonymous_mod <- subsetByOverlaps(
# #     cds,
# #     coding_effective
# # )
# # cds_synonymous_mod <- subsetByOverlaps(
# #     cds,
# #     coding[(mcols(coding)$CONSEQUENCE) == "synonymous"]
# # )
# exons_mod <- subsetByOverlaps(
#     exons,
#     rowRanges(vcf_hom)
# )
# write.csv(
#     as.data.frame(exons_mod),
#     "tests/test-VariantAnnotation2/exons_mod_gtf.csv",
#     row.names = FALSE,
#     quote = FALSE
# )

# # Get CDS or exon sequences
# # cds_seqs_nonsynonymous_mod <- extractTranscriptSeqs(ref_mod, cds_nonsynonymous_mod)
# # cds_seqs_synonymous_mod <- extractTranscriptSeqs(ref_mod, cds_synonymous_mod)
# exons_seqs_mod <- extractTranscriptSeqs(ref_mod, exons_mod)
# exons_seqs <- extractTranscriptSeqs(ref, exons_mod)
# ```

# 5. Write to FASTA
# ```{r}
# # writeXStringSet(cds_seqs_nonsynonymous_mod, "tests/test-VariantAnnotation/modified_cds_nonsynonymous.fa")
# # writeXStringSet(cds_seqs_synonymous_mod, "tests/test-VariantAnnotation/modified_cds_synonymous.fa")
# writeXStringSet(exons_seqs_mod, "tests/test-VariantAnnotation2/edited_exons.fa")
# writeXStringSet(exons_seqs, "tests/test-VariantAnnotation2/original_exons.fa")
# ```

# 6. Use `seqkit translate` to get protein sequences (in shell)
# ```{bash}
# # performs 6-frame translation on modified exons
# docker run -v $(pwd):/data -w /data -u $(id -u) --rm staphb/seqkit:2.10.0 \
#     seqkit translate \
#         -f 6 \
#         -F \
#         --min-len 30 \
#         tests/test-VariantAnnotation2/edited_exons.fa \
#     > tests/test-VariantAnnotation2/edited_protein.fa

# docker run -v $(pwd):/data -w /data -u $(id -u) --rm staphb/seqkit:2.10.0 \
#     seqkit translate \
#         -f 6 \
#         -F \
#         --out-subseqs \
#         --min-len 1 \
#         --trim \
#         tests/test-VariantAnnotation2/edited_exons.fa \
#     > tests/test-VariantAnnotation2/edited_protein_6frame_all_ORF.fa

# # performs 6-frame translation on original exons
# docker run -v $(pwd):/data -w /data -u $(id -u) --rm staphb/seqkit:2.10.0 \
#     seqkit translate \
#         -f 6 \
#         -F \
#         --min-len 30 \
#         tests/test-VariantAnnotation2/original_exons.fa \
#     > tests/test-VariantAnnotation2/original_protein.fa
# ```

# 7. Keep either the "first" of the "longest" open reading frame
# ```{r}
# find_longest_orf <- function(seq) {
#   matches <- gregexpr("M[^\\*]*\\*", seq)[[1]]
#   if (matches[1] == -1) return(NA)

#   match_lengths <- attr(matches, "match.length")
#   substrings <- substring(seq, matches, matches + match_lengths - 2)

#   substrings[which.max(nchar(substrings))]
# }
# ```
# ```{r}
# # Read modified protein sequences
# original_proteins <- readAAStringSet("tests/test-VariantAnnotation2/original_protein.fa")
# modified_proteins <- readAAStringSet("tests/test-VariantAnnotation2/edited_protein.fa")

# # fast filter protein without M
# modified_proteins_with_M <- modified_proteins[grepl("M", modified_proteins)]; length(modified_proteins_with_M)

# # extract protein from 1st ORF (probably uORF)
# ltrimmed_proteins <- sub("[^M]*M", "M", modified_proteins_with_M)
# ltrimmed_proteins_with_stop <- ltrimmed_proteins[grepl("[*]", ltrimmed_proteins)]; length(ltrimmed_proteins_with_stop)
# rtrimmed_proteins <- sub("[*].*", "", ltrimmed_proteins_with_stop)
# final_proteins <- rtrimmed_proteins[nchar(rtrimmed_proteins) >= 30]; length(final_proteins)
# final_proteins <- AAStringSet(final_proteins)

# # Write first ORF proteins to FASTA
# writeXStringSet(final_proteins, "tests/test-VariantAnnotation2/edited_protein_6frame_first_ORF.fa")

# # extract longest ORF
# longest_protein <- sapply(modified_proteins_with_M, function(x){find_longest_orf(as.character(x))})
# longest_AAStringSet <- AAStringSet(longest_protein[!is.na(longest_protein)])

# # write longest ORF protein to FASTA
# writeXStringSet(longest_AAStringSet, "tests/test-VariantAnnotation2/edited_protein_6frame_longest_ORF.fa")
# ```

# ```{r}
# # Read original protein sequences
# original_proteins <- readAAStringSet("tests/test-VariantAnnotation/original_exons_protein.fa")
# original_proteins_with_M <- original_proteins[grepl("M", original_proteins)]
# ltrimmed_proteins <- sub("[^M]*M", "M", original_proteins_with_M)
# ltrimmed_proteins_with_stop <- ltrimmed_proteins[grepl("[*]", ltrimmed_proteins)]
# rtrimmed_proteins <- sub("[*].*", "", ltrimmed_proteins_with_stop)
# final_proteins <- rtrimmed_proteins[nchar(rtrimmed_proteins) >= 30]
# final_proteins <- AAStringSet(final_proteins)

# # Write final proteins to FASTA
# writeXStringSet(final_proteins, "tests/test-VariantAnnotation/original_protein_6frames.fa")
# ```
