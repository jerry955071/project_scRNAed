# Snakemake implementation of the souporcell workflow
# Reference: https://github.com/wheaton5/souporcell

# Note: Modularizing Snakemake workflow
# Reference: https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html
# Modularization in Snakemake comes at four different levels.
# level 1: wrapper
# level 2: small Snakefiles
# level 3: module statements

# load config file
configfile: "configs/config.json"

# include sub-workflows
include: "CellRanger.smk"
include: "Remapping.smk"
include: "VariantCalling.smk"
# include: "VariantCalling-DNA.smk"
# include: "VariantCalling-DNA-Freebayes-by-sample.smk"
include: "VariantCalling-DNA-GATK.smk"
include: "LCM-RNA-editing.smk"
include: "Bulk-scRNA-editing.smk"

# call for output files
rule VariantCalling:
    input:
        alt_mtx=expand(
            "outputs/VariantCalling/vartrix/{sample}/alt.mtx",
            sample=[i["name"] for i in config["samples"]]
        ),
        ref_mtx=expand(
            "outputs/VariantCalling/vartrix/{sample}/ref.mtx",
            sample=[i["name"] for i in config["samples"]]
        ),
        snv_loci=expand(
            "outputs/VariantCalling/vawk/{sample}.snv.loci.txt",
            sample=[i["name"] for i in config["samples"]]
        ),
        barcodes=expand(
            "outputs/Remapping/renamer/{sample}/barcodes.tsv",
            sample=[i["name"] for i in config["samples"]]
        )

rule VariantCalling_DNA_GATK_egr:
    input:
        hom_alt="outputs/VariantCalling-DNA/gatk_joint/egr/hom_alt.vcf",
        hom_ref="outputs/VariantCalling-DNA/gatk_joint/egr/hom_ref.vcf",
        het="outputs/VariantCalling-DNA/gatk_joint/egr/het.vcf",
        snp="outputs/VariantCalling-DNA/gatk_joint/egr/snp.vcf",


rule ptr_tenx_batch1:
    input:
        alt_mtx="outputs/VariantCalling/vartrix/ptr_tenx_batch1/alt.mtx",
        ref_mtx="outputs/VariantCalling/vartrix/ptr_tenx_batch1/ref.mtx",
        snv_loci="outputs/VariantCalling/vawk/ptr_tenx_batch1.snv.loci.txt",
        barcodes="outputs/Remapping/renamer/ptr_tenx_batch1/barcodes.tsv",
        gatk_output="outputs/VariantCalling-DNA/gatk_joint/Ptr.filtered.vcf.gz"

rule lcm_rna_editing:
    input:
        expand(
            "outputs/LCM-RNA-editing/subset-reditable/{sample}-subset.tsv",
            sample=[
                "fiber-rep1",
                "fiber-rep2",
                "fiber-rep3",
                "vessel-rep1",
                "vessel-rep2",
                "vessel-rep3",
                "ray-rep1",
                "ray-rep2",
                "ray-rep3"
            ]
        )
        
rule categorized_sites:
    input:
        hom_alt="outputs/VariantCalling-DNA/gatk_joint/hom_alt.vcf",
        hom_ref="outputs/VariantCalling-DNA/gatk_joint/hom_ref.vcf",
        het="outputs/VariantCalling-DNA/gatk_joint/het.vcf",
        snp="outputs/VariantCalling-DNA/gatk_joint/snp.vcf"



# all-purpose rules
# =========================================================
rule samtools_index_any:
    input:
        "{any}.bam"
    output:
        "{any}.bam.bai"
    log:
        "logs/Snakefile/samtools_index/{any}.bam.log"
    shell:
        "samtools index {input} 1> {log} 2> {log}"


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
