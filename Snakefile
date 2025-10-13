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
include: "Utilities.smk"
include: "CellRanger.smk"
include: "Remapping.smk"
include: "VariantCalling.smk"
include: "VariantCalling-DNA-GATK.smk"
include: "LCM-RNA-editing.smk"
include: "Bulk-scRNA-editing.smk"
include: "VariantAnnotation.smk"

# Custom functions used by all workflows
from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    """Return the first dictionary in a list of dictionaries where the value of key k matches v."""
    return [x for x in d if x[k] == v][0]

def query_all(d:List[dict], k:str, v:str, k_out:str) -> List[str]:
    """Return a list of values from key k_out in a list of dictionaries where the value of key k matches v."""
    return [x[k_out] for x in d if x[k] == v]

def _genome_assembly(wildcards):
    """Get genome assembly file path by sample name"""
    species = query(
        config["samples"] + \
        config["samples-dna"] + \
        config["samples-lcm"], 
        "name", 
        wildcards.sample
    )["species"]
    return query(config["references"], "species", species)["assembly"]

# wildcard constraints
wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]]),
    species="|".join([i["species"] for i in config["references"]]),


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


# rules
rule all:
    input:
        "outputs/Snakefile/VariantCalling/done.txt",
        "outputs/Snakefile/VariantCalling-DNA-GATK/done.txt",
        "outputs/Snakefile/LCM-RNA-editing/done.txt"


# call scRNA polymorphism per single-cell samples
SC_SAMPLES = [i["name"] for i in config["samples"]]
rule VariantCalling:
    input:
        alt_mtx=expand(
            "outputs/VariantCalling/vartrix/{sample}/alt.mtx",
            sample=SC_SAMPLES
        ),
        ref_mtx=expand(
            "outputs/VariantCalling/vartrix/{sample}/ref.mtx",
            sample=SC_SAMPLES
        ),
        snv_loci=expand(
            "outputs/VariantCalling/vawk/{sample}.snv.loci.txt",
            sample=SC_SAMPLES
        ),
        barcodes=expand(
            "outputs/Remapping/renamer/{sample}/barcodes.tsv",
            sample=SC_SAMPLES
        )
    output:
        "outputs/Snakefile/VariantCalling/done.txt"
    shell:
        """
        echo date +%Y-%m-%d_%H:%M:%S > {output}
        echo Created files: >> {output}
        echo {input} >> {output}
        """

# perform genotyping (mostly for identifying homozygous sites) at scRNA polymorphic sites per species
rule VariantCalling_DNA_GATK:
    input:
        "outputs/VariantCalling-DNA/gatk_joint/ptr/hom_ref.vcf",
        "outputs/VariantCalling-DNA/gatk_joint/egr/hom_ref.vcf"
    output:
        "outputs/Snakefile/VariantCalling-DNA-GATK/done.txt"
    shell:
        """
        echo date +%Y-%m-%d_%H:%M:%S > {output}
        echo Created files: >> {output}
        echo {input} >> {output}
        """


# characterize RNA editing profile of fiber, vessel, and ray cells
rule LCM_RNA_editing:
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
    output:
        "outputs/Snakefile/LCM-RNA-editing/done.txt"
    shell:
        """
        echo date +%Y-%m-%d_%H:%M:%S > {output}
        echo Created files: >> {output}
        echo {input} >> {output}
        """
        



