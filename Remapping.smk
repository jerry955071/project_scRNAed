# wildcard_constraints:
#     sample="|".join([i["name"] for i in config["samples"]]),
#     species="|".join([i["species"] for i in config["references"]])


# ================== Rename BAM to FASTQ ==================
rule renamer:
    threads: 1
    input:
        bam="outputs/CellRanger/count/{sample}/outs/possorted_genome_bam.bam",
        barcodes_gz="outputs/CellRanger/count/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        barcodes="outputs/Remapping/renamer/{sample}/barcodes.tsv",
        fq="outputs/Remapping/renamer/{sample}/fq.fq"
    log:
        "logs/Remapping/renamer/{sample}.log"
    shell:
        """
        # unzip barcodes.tsv.gz
        gunzip -c {input.barcodes_gz} 1> {output.barcodes} 2> {log}

        # rename bam to fastq
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/pysam:0.23.0 \
                python src/souporcell/renamer.py \
                --bam {input.bam} \
                --barcodes {output.barcodes} \
                --out {output.fq} \
            1>> {log} \
            2>> {log}
        """

# ================== Remapping with minimap2 ==================
rule minimap2:
    threads: 24
    input:
        ref=lambda wildcards: _genome_assembly(wildcards),
        fq="outputs/Remapping/renamer/{sample}/fq.fq"
    output:
        sam="outputs/Remapping/minimap2/{sample}/minimap.sam"
    log:
        "logs/Remapping/minimap2/{sample}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            staphb/minimap2:2.28 \
                minimap2 \
                    -ax splice -t {threads} -G50k -k 21 -w 11 --sr \
                    -A2 -B8 -O12,32 -E2,1 -r200 -p.5 -N20 -f1000,5000 \
                    -n2 -m20 -s40 -g2000 -2K50m --secondary=no \
                    -o {output.sam} \
                    {input.ref} {input.fq} \
                1> {log} \
                2> {log}
        """

# ================== retag SAM to BAM ==================
rule retag:
    input:
        sam="outputs/Remapping/minimap2/{sample}/minimap.sam"
    output:
        tagged_bam="outputs/Remapping/retag/{sample}/minitagged.bam"
    log:
        "logs/Remapping/retag/{sample}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            chiaenu/pysam:0.23.0 \
                python src/souporcell/retag.py \
                    --sam {input.sam} \
                    --out {output.tagged_bam} \
            1> {log} \
            2> {log}
        """

# ================== sort BAM with samtools ==================
rule samtools_sort:
    threads: 24
    input:
        tagged_bam="outputs/Remapping/retag/{sample}/minitagged.bam"
    output:
        sorted_bam="outputs/Remapping/samtools/{sample}/minitagged_sorted.bam"
    log:
        "logs/Remapping/samtools-sort/{sample}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools sort \
                -@ {threads} \
                -o {output.sorted_bam} \
                {input.tagged_bam} \
                1> {log} \
                2> {log}
        """

# ================== index BAM with samtools ==================
# NOTE: replaced by the all purpose rule `samtools_index_any` in Utilities.smk
rule samtools_index:
    threads: 24
    input:
        sorted_bam="outputs/Remapping/samtools/{sample}/minitagged_sorted.bam"
    output:
        index="outputs/Remapping/samtools/{sample}/minitagged_sorted.bam.bai"
    log:
        "logs/Remapping/samtools-index/{sample}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools index \
                -@ {threads} \
                {input.sorted_bam} \
                1> {log} \
                2> {log}
        """