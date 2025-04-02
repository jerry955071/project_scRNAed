# First time trying modularizing Snakemake workflow
# Reference: https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html
# Modularization in Snakemake comes at four different levels.
# level 1: wrapper
# level 2: small Snakefiles
# level 3: module statements

include: "CellRanger.smk"

rule all:
    input:
        "outputs/CellRanger/count/ptr_tenx_batch1/outs/possorted_genome_bam.bam",
