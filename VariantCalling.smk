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

# ================= Custom functions =================
# query: query from list of dict by key-value pair
from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    return [x for x in d if x[k] == v][0]

# get_ref_by_sample: get reference genome for sample (used in minimap2 rule)
def get_ref_by_sample(wildcards):
    species = query(config["samples"], "name", wildcards.sample)["species"]
    return query(config["references"], "species", species)["assembly"]

# call_chromosomes_vcf_by_sample: call for vcf files per chromosome (used in bcftools_merge rule)
def call_chromosomes_vcf_by_sample(wildcards):
    species = query(config["samples"], "name", wildcards.sample)["species"]
    chromosomes = query(config["references"], "species", species)["chromosomes"]
    return [
        f"outputs/VariantCalling/freebayes/{wildcards.sample}/{chrom}.vcf"
        for chrom in chromosomes
    ]

# ==== Calling candidate variants using freebayes ====
# 1. Calling candidates per chromosome
rule freebayes_by_chromosome:
    input:
        minitagged_sorted_bam="outputs/Remapping/samtools/{sample}/minitagged_sorted.bam",
        reference_fasta=get_ref_by_sample
    output:
        vcf="outputs/VariantCalling/freebayes/{sample}/{chromosome}.vcf"
    log:
        "logs/VariantCalling/freebayes/{sample}/{chromosome}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            --rm \
            chiaenu/freebayes:1.3.6-1 \
                freebayes \
                    -f {input.reference_fasta} \
                    -iXu \
                    -C 2 \
                    -q 20 \
                    -n 3 \
                    -E 1 \
                    -m 30 \
                    --min-coverage 6 \
                    --limit-coverage 100000 \
                    {input.minitagged_sorted_bam} \
                    1> {output.vcf} \
                    2> {log}
        """

# 2. Merging vcf files
rule bcftools_merge:
    input:
        call_chromosomes_vcf_by_sample
    output:
        vcf="outputs/VariantCalling/bcftools_merge/{sample}.vcf"
    log:
        "logs/VariantCalling/bcftools_merge/{sample}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            --rm \
            biocontainers/bcftools:v1.9-1-deb_cv1 \
                bcftools merge \
                    -O v \
                    -o {output.vcf} \
                    {input}
        """

# ==== Calling single-cell variants using vartrix ====
rule get_vartrix:
    output:
        vartrix="src/vartrix/vartrix_linux"
    log:
        "logs/VariantCalling/get_vartrix.log"
    shell:
        """
        wget https://github.com/10XGenomics/vartrix/releases/download/v1.1.22/vartrix_linux
        chmod +x vartrix_linux
        mv vartrix_linux src/vartrix/vartrix_linux
        """

rule vartrix:
    threads: 24
    input:
        bam="outputs/Remapping/samtools/{sample}/minitagged_sorted.bam",
        barcodes="outputs/Remapping/renamer/{sample}/barcodes.tsv",
        vcf="outputs/VariantCalling/bcftools_merge/{sample}.vcf",
        fasta="outputs/Remapping/renamer/{sample}/fq.fq"
    output:
        ref_matrix="outputs/VariantCalling/vartrix/{sample}/ref.mtx",
        alt_matrix="outputs/VariantCalling/vartrix/{sample}/alt.mtx"
    shell:
        """
        src/vartrix/vartrix_linux \
            --umi \
            --mapq 30 \
            -b {input.bam} \
            -c {input.barcodes} \
            --scoring-method coverage \
            --threads {threads} \
            --ref-matrix {output.ref_matrix} \
            --out-matrix {output.alt_matrix} \
            -v {input.vcf} \
            --fasta {input.fasta} \
            1> {log} \
            2> {log}
        """