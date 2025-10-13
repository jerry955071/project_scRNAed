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

# TODO: 
# Put tests/test-VariantAnnotation/Variant-annotation.Rmd &
# make-modified-ref.sh here

# convert concatenated vcf to bcf files
rule vcf_to_bcf:
    input:
        "outputs/VariantCalling/bcftools_concat/{sample}.vcf"
    output:
        "outputs/VariantCalling/bcftools_concat/{sample}.bcf"
    log:
        "logs/VariantAnnotation/vcf_to_bcf/{sample}.log"
    shell:
        """
        # convert to bcf
        docker run {docker_mount} -u $(id -u) --rm staphb/bcftools:1.21 \
            bcftools view -Ob {input} -o {output} \
        1> {log} \
        2> {log}
        
        # index
        docker run {docker_mount} -u $(id -u) --rm staphb/bcftools:1.21 \
            bcftools index {output} \
        1>> {log} \
        2>> {log}
        """


# # bgzip concatenated bcf files
# rule bgzip_concatenated:
#     input:
#         "outputs/VariantCalling/bcftools_concat/{sample}.bcf"
#     output:
#         "outputs/VariantCalling/bcftools_concat/{sample}.bcf.gz"
#     log:
#         "logs/VariantAnnotation/bgzip_concatenated/{sample}.log"
#     shell:
#         """
#         docker run {docker_mount} -u $(id -u) --rm staphb/htslib:1.21 \
#                 bgzip -c {input} \
#             1> {output} \
#             2> {log}
#         """

# # tabix concatenated bcf.gz files
# rule tabix_concatenated:
#     input:
#         "outputs/VariantCalling/bcftools_concat/{sample}.bcf.gz"
#     output:
#         "outputs/VariantCalling/bcftools_concat/{sample}.bcf.gz.tbi"
#     log:
#         "logs/VariantAnnotation/tabix_concatenated/{sample}.log"
#     shell:
#         """
#         docker run {docker_mount} -u $(id -u) --rm staphb/htslib:1.21 \
#                 tabix -s 1 -b 2 -e 2 -c 'R' {input} \
#             1> {log} \
#             2> {log}
#         """


# bcftools merge
rule bcftools_merge:
    input:
        lambda wildcards: expand(
            "outputs/VariantCalling/bcftools_concat/{sample}.bcf",
            sample=query_all(
                d=config["samples"],
                k="species",
                v=wildcards.species,
                k_out="name"
            )
        )
    output:
        "outputs/VariantAnnotation/bcftools_merge/{species}.vcf"
    log:
        "logs/VariantAnnotation/bcftools_merge/{species}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm staphb/bcftools:1.21 \
            bcftools merge \
                -Ob \
                -o {output} \
                --force-samples \
                {input} \
            1> {log} \
            2> {log}
        """

# bcftools csq
rule bcftools_csq:
    input:
        bcf="outputs/VariantAnnotation/bcftools_merge/{species}.bcf",
        reference=lambda wildcards: query(config["references"], "species", wildcards.species)["assembly"],
        gff3=lambda wildcards: query(config["references"], "species", wildcards.species)["annotation-gff"]
    output:
        gff4csq="outputs/VariantAnnotation/bcftools_csq/{species}.gff.gz",
        bcf="outputs/VariantAnnotation/bcftools_csq/{species}.bcf"
    log:
        "logs/VariantAnnotation/bcftools_csq/{species}.log"
    shell:
        """
        # covert gff to format recognizable by bcftools-csq
        cat {input.gff3} | src/bcftools-1.22/misc/gff2gff | gzip -c > {output.gff4csq}

        # Basic usage
        docker run {docker_mount} -u $(id -u) --rm staphb/bcftools:1.21 \
            bcftools csq \
                -f {input.reference} \
                -g {output.gff4csq} \
                -Ob \
                -o {output.bcf} \
                -v 2 \
                {input.bcf} \
            1> {log} \
            2> {log}

        # bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' out.bcf
        """

