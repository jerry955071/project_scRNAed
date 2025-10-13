# rule samtools_index_any:
#     input:
#         "{any}.bam"
#     output:
#         "{any}.bam.bai"
#     log:
#         "logs/Snakefile/samtools_index/{any}.bam.log"
#     shell:
#         "samtools index {input} 1> {log} 2> {log}"

# NOTE: This rule is confusing Snakemake, and I forget when it was used. Thus commented until needed.
# rule bgzip:
#     input:
#         "{file}"
#     output:
#         "{file}.gz"
#     log:
#         "logs/gzip/{file}.log"
#     shell:
#         """
#         docker run \
#             {docker_mount} \
#             -u $(id -u) \
#             --rm \
#             staphb/htslib:1.21 \
#                 bgzip -c {input} \
#             1> {output} \
#             2> {log}
#         """
