# Snakemake setups
configfile: "configs/config.json"
wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]]),
    species="|".join([i["species"] for i in config["references"]]),
    random_seed="[0-9]+",
    gff_ext="gff|gff3"

# set docker mount points
docker_mount = ""
for volume in config["volumes"]:
    docker_mount += "-v %s:%s:%s " % (
        volume["host"],
        volume["container"],
        volume["mode"]
    )
    if volume["is_workspace"]:
        docker_mount += "-w %s " % volume["container"]


# rule "all"
rule bulk_scRNA_editing:
    input:
        using_mimimap_bam="outputs/bulk-scRNA-editing/reditools1-minimap2/ptr_tenx_batch1.txt",
        using_star_bam="outputs/bulk-scRNA-editing/reditools1-star/ptr_tenx_batch1.txt"

# 1. extract bam records from Cell Ranger output that contribute to UMI
# Reads contribute to UMI counting:
# (1) The read is confidently mapped in the sense orientation to a single gene
# (2) This read is representative for the molecule and can be treated as a UMI count
# (3) Confidently assigned feature barcode
# Reference: https://www.10xgenomics.com/cn/analysis-guides/tutorial-navigating-10x-barcoded-bam-files
rule extract_umi_reads:
    input:
        cell_ranger_bam="outputs/CellRanger/count/{sample}/outs/possorted_genome_bam.bam"
    output:
        umi_reads="outputs/bulk-scRNA-editing/extract_umi_reads/{sample}.bam"
    log:
        "logs/bulk-scRNA-editing/extract_umi_reads/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm biocontainers/samtools:v1.9-4-deb_cv1 \
            scripts/extract-umi-reads.sh \
                {input} {output} \
                > {log} 2>&1
        """


# 2. Call RNA editing using REDItools, where strand are inferred by gene annotation
# -X	Sort annotation files. It requires grep and sort unix executables.
# -G	Infer strand by GFF annotation (must be GFF and sorted, otherwise use -X). 
#       Sorting requires grep and sort unix executables.
rule reditools1_star:
    threads: 32
    input:
        bam="outputs/bulk-scRNA-editing/extract_umi_reads/{sample}.bam",
        assembly="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa",
        annotation="references/Ptr/annotation/Ptrichocarpa_533_v4.1.gene_exons.gtf.gz"
    output:
        out_dir=directory("outputs/bulk-scRNA-editing/reditools1-star/{sample}/"),
        reditable="outputs/bulk-scRNA-editing/reditools1-star/{sample}.txt"
    params:
        mapq=255
    log:
        "logs/bulk-scRNA-editing/reditools1-star/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm chiaenu/reditools1:1.3 \
            REDItoolDnaRna.py \
                -i {input.bam} \
                -f {input.assembly} \
                -t {threads} \
                -o {output.out_dir} \
                -m 0,{params.mapq} \
                -G {input.annotation} \
                > {log} 2>&1
        cp {output.out_dir}/DnaRna_*/outTable_* {output.reditable}
        """

rule reditools1_star_all:
    threads: 32
    input:
        bam="outputs/CellRanger/count/{sample}/outs/possorted_genome_bam.bam",
        assembly="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa",
        annotation="references/Ptr/annotation/Ptrichocarpa_533_v4.1.gene_exons.gtf.gz"
    output:
        out_dir=directory("outputs/bulk-scRNA-editing/reditools1-star-all/{sample}/"),
        reditable="outputs/bulk-scRNA-editing/reditools1-star-all/{sample}.txt"
    params:
        mapq=255
    log:
        "logs/bulk-scRNA-editing/reditools1-star-all/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm chiaenu/reditools1:1.3 \
            REDItoolDnaRna.py \
                -i {input.bam} \
                -f {input.assembly} \
                -t {threads} \
                -o {output.out_dir} \
                -m 0,{params.mapq} \
                -G {input.annotation} \
                > {log} 2>&1
        cp {output.out_dir}/DnaRna_*/outTable_* {output.reditable}
        """

rule reditools1_minimap:
    threads: 32
    input:
        bam="outputs/Remapping/samtools/{sample}/minitagged_sorted.bam",
        assembly="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa",
        annotation="references/Ptr/annotation/Ptrichocarpa_533_v4.1.gene_exons.gtf.gz"
    output:
        out_dir=directory("outputs/bulk-scRNA-editing/reditools1-minimap2/{sample}/"),
        reditable="outputs/bulk-scRNA-editing/reditools1-minimap2/{sample}.txt"
    params:
        mapq=30
    log:
        "logs/bulk-scRNA-editing/reditools1-minimap2/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm chiaenu/reditools1:1.3 \
            REDItoolDnaRna.py \
                -i {input.bam} \
                -f {input.assembly} \
                -t {threads} \
                -o {output.out_dir} \
                -m 0,{params.mapq} \
                -G {input.annotation} \
                > {log} 2>&1
        cp {output.out_dir}/DnaRna_*/outTable_* {output.reditable}
        """