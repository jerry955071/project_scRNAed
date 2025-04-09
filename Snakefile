# Snakemake implementation of the souporcell workflow
# Reference: https://github.com/wheaton5/souporcell

# Note: Modularizing Snakemake workflow
# Reference: https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html
# Modularization in Snakemake comes at four different levels.
# level 1: wrapper
# level 2: small Snakefiles
# level 3: module statements

# include sub-workflows
include: "CellRanger.smk"
include: "Remapping.smk"
include: "VariantCalling.smk"

# call for output files
rule all:
    input:
        "outputs/VariantCalling/vartrix/ptr_tenx_batch1/alt.mtx"


