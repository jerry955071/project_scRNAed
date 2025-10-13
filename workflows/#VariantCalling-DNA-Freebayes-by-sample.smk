# Snakemake setups
configfile: "configs/config.json"
wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples"]] + [i["name"] for i in config["samples-dna"]]),
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



# ====== add read group ======
rule samtools_addreplacerg:
    input:
        "/home/b05b01002/HDD/project_RE/outputs/sorted_bam/Ptr/DNA/{sample}.sorted.bam"
    output:
        "outputs/VariantCalling-DNA/samtools_addreplacerg/{sample}.sorted.addrg.bam"
    log:
        "logs/VariantCalling-DNA/samtools_addreplacerg/{sample}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            biocontainers/samtools:v1.9-4-deb_cv1 \
                samtools addreplacerg \
                    -r 'ID:{wildcards.sample}' \
                    -r 'SM:{wildcards.sample}' \
                    -o {output} \
                    {input} \
            > {log}
        """


# ====== mark duplicates (PCR duplicates) ======
rule sambamba_markdup_by_sample:
    threads: 8
    input:
        "outputs/VariantCalling-DNA/samtools_addreplacerg/{sample}.sorted.addrg.bam"
    output:
        marked_bam="outputs/VariantCalling-DNA/sambamba_markdup/{sample}.sorted.addrg.markdup.bam",
        temp_outdir=directory("outputs/VariantCalling-DNA/sambamba_markdup/{sample}.tmp")
    log:
        "logs/VariantCalling-DNA/sambamba_markdup/{sample}.log"
    shell:
        """
        src/sambamba/sambamba-1.0.1-linux-amd64-static \
            markdup \
            -t {threads} \
            -p \
            --tmpdir {output.temp_outdir} \
            {input} {output.marked_bam} \
            2> {log}
            1> {log}
        """

# ==== Calling candidate variants using freebayes ====
# 1. Calling candidates per chromosome
rule freebayes_by_chromosome_DNA_by_sample:
    input:
        minitagged_sorted_bam=expand(
            "outputs/VariantCalling-DNA/sambamba_markdup/{sample}.sorted.addrg.markdup.bam",
            sample=[i["name"] for i in config["samples-dna"]]
        ),
        reference_fasta="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa"
    output:
        vcf="outputs/VariantCalling-DNA/freebayes-by-sample/{chromosome}.vcf"
    log:
        "logs/VariantCalling-DNA/freebayes-by-sample/{chromosome}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            --rm \
            chiaenu/freebayes:1.3.6-1 \
                freebayes \
                    -f {input.reference_fasta} \
                    -r {wildcards.chromosome} \
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

# ==== This rule is defined in VariantCalling.smk ====
# # 2. Compressing vcf files with bgzip (will create temporary .gz files)
# rule bgzip:
#     input:
#         "{file}"
#     output:
#         temp("{file}.gz")
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
# ====================================================

# 3. Merging vcf files (using gzipped vcf files)
rule bcftools_concat_DNA_by_sample:
    input:
        expand(
            "outputs/VariantCalling-DNA/freebayes-by-sample/{chromosome}.vcf",
            chromosome=config["references"][0]["chromosomes"]
        )
    output:
        vcf="outputs/VariantCalling-DNA/bcftools_concat/Ptr-by-sample.vcf"
    log:
        "logs/VariantCalling-DNA/bcftools_concat/Ptr-by-sample.log"
    shell:
        """
        # merge vcf files
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            staphb/bcftools:1.21 \
                bcftools concat \
                    -O v \
                    -o {output.vcf} \
                    {input} \
            1> {log} \
            2> {log}
        """

# 4. Extract variants from the merged vcf file 
rule vawk_filter_extract_by_sample:
    input:
        "outputs/VariantCalling-DNA/bcftools_concat/Ptr-by-sample.vcf"
    output:
        "outputs/VariantCalling-DNA/vawk_filter_extract/Ptr-by-sample.snv.loci.txt"
    log:
        "logs/VariantCalling-DNA/vawk_filter_extract/Ptr-by-sample.log"
    shell:
        """
        # write header to the output file
        echo -e "CHROM:POS\tREF\tALT\tQUAL\tDP\tRO\tAO\tQR\tQA\tAF" > {output}

        # extract variants using vawk
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
                chiaenu/vawk:0.0.2 \
                    vawk \
                        '{{print $1":"$2,$4,$5,$6,I$DP,I$RO,I$AO,I$QR,I$QA,I$AF}}' \
                        {input} \
            1>> {output} \
            2> {log}
        """
