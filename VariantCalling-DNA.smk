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


# ====== mark duplicates (PCR duplicates) ======
rule sambamba_markdup:
    threads: 8
    input:
        "/home/b05b01002/HDD/project_RE/outputs/sorted_bam/Ptr/DNA/Ptr.sorted.bam"
    output:
        marked_bam="outputs/VariantCalling-DNA/sambamba_markdup/Ptr.sorted.markdup.bam",
        temp_outdir=directory("outputs/VariantCalling-DNA/sambamba_markdup/Ptr.tmp")
    log:
        "logs/VariantCalling-DNA/sambamba_markdup/Ptr.log"
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
rule freebayes_by_chromosome_DNA:
    input:
        minitagged_sorted_bam="outputs/VariantCalling-DNA/sambamba_markdup/Ptr.sorted.markdup.bam",
        reference_fasta="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa"
    output:
        vcf="outputs/VariantCalling-DNA/freebayes/{chromosome}.vcf"
    log:
        "logs/VariantCalling-DNA/freebayes/{chromosome}.log"
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
rule bcftools_merge_DNA:
    input:
        expand(
            "outputs/VariantCalling-DNA/freebayes/{chromosome}.vcf",
            chromosome=config["references"][0]["chromosomes"]
        )
    output:
        vcf="outputs/VariantCalling-DNA/bcftools_merge/Ptr.vcf"
    log:
        "logs/VariantCalling-DNA/bcftools_merge/Ptr.log"
    shell:
        """
        # merge vcf files
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            staphb/bcftools:1.21 \
                bcftools merge \
                    -O v \
                    -o {output.vcf} \
                    --no-index \
                    --force-samples \
                    {input} \
            1> {log} \
            2> {log}
        """

# 4. Extract variants from the merged vcf file 
# Only whose quality is above a certain threshold
# Output fields:
# (1) CHROM:POS (chromosome:position)
# (2) DP (depth of coverage)
# (5) AF (alternate allele frequency)
rule vawk_filter_extract:
    input:
        "outputs/VariantCalling-DNA/bcftools_merge/Ptr.vcf"
    output:
        "outputs/VariantCalling-DNA/vawk_filter_extract/Ptr.snv.loci.txt"
    log:
        "logs/VariantCalling-DNA/vawk_filter_extract/Ptr.log"
    params:
        qual_threshold=20
    shell:
        """
        # write header to the output file
        echo -e "CHROM:POS\tDP\tAF" > {output}

        # extract variants using vawk
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
                chiaenu/vawk:0.0.2 \
                    vawk \
                        '{if ($6>{params.qual_threshold}) print $1:$2,I$DP,I$AF}' \
                        {input} \
            1>> {output} \
            2> {log}
        """
