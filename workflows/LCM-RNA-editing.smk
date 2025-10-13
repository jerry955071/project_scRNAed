# Reference: https://github.com/BioinfoUNIBA/REDItools2?tab=readme-ov-file#7-reditools-20-options
# Snakemake setups
configfile: "configs/config.json"
wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples-lcm"]]),
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

# function: get fastq by sample name (used in rule: fastp)
def get_fastq_by_sample(wildcards):
    """
    Get the fastq file names by sample name.
    """
    for sample in config["samples-lcm"]:
        if sample["name"] == wildcards.sample:
            return {
                "r1": sample["fastq_r1"],
                "r2": sample["fastq_r2"]
            }
    raise ValueError("Sample not found: %s" % wildcards.sample)

# The standard REDItools workflow
# 1. Run REDItools on RNA generating RNA.tsv
# 2. Run REDItools on DNA generating DNA.tsv
# 3. Annotate RNA.tsv with DNA.tsv
# 
# We only run step 1 because we've done variant calling with DNA data

# ==== 0. Pre-process FASTQ files ====
rule fastp:
    threads: 16
    input:
        r1=lambda wildcards: get_fastq_by_sample(wildcards)["r1"],
        r2=lambda wildcards: get_fastq_by_sample(wildcards)["r2"]
    output:
        r1="outputs/LCM-RNA-editing/fastp/{sample}_1.fq.gz",
        r2="outputs/LCM-RNA-editing/fastp/{sample}_2.fq.gz"
    log:
        "logs/LCM-RNA-editing/fastp/{sample}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) staphb/fastp:0.24.1 \
            fastp \
                --in1 {input.r1} \
                --in2 {input.r2} \
                --out1 {output.r1} \
                --out2 {output.r2} \
                --compression 4 \
                --verbose \
                --detect_adapter_for_pe \
                --cut_front \
                --cut_tail \
                --correction \
                --cut_window_size 4 \
                --cut_mean_quality 20 \
                --length_required 50 \
                --thread {threads} \
                --json /dev/null \
                --html /dev/null \
            > {log} 2>&1
        """

# ==== 1. Align RNA reads to genome ====
rule star_index:
    threads: 1
    input:
        assembly="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa",
        gtf="references/Ptr/annotation/Ptrichocarpa_533_v4.1.gene_exons.gtf"
    output:
        directory("references/Ptr/star_index")
    log:
        "logs/LCM-RNA-editing/star_index.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) alexdobin/star:2.7.10a_alpha_220506 \
            STAR \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir {output} \
                --genomeFastaFiles {input.assembly} \
                --sjdbGTFfile {input.gtf} \
                --sjdbOverhang 149 \
            > {log} 2>&1
        """

rule star_align:
    threads: 16
    input:
        r1="outputs/LCM-RNA-editing/fastp/{sample}_1.fq.gz",
        r2="outputs/LCM-RNA-editing/fastp/{sample}_2.fq.gz",
        index="references/Ptr/star_index"
    output:
        outdir=directory("outputs/LCM-RNA-editing/star/{sample}"),
        bam="outputs/LCM-RNA-editing/star/{sample}/Aligned.sortedByCoord.out.bam"
    log:
        "logs/LCM-RNA-editing/star-align/{sample}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) alexdobin/star:2.7.10a_alpha_220506 \
            STAR \
                --runThreadN {threads} \
                --genomeDir {input.index} \
                --readFilesIn {input.r1} {input.r2} \
                --readFilesCommand gunzip -c \
                --outFileNamePrefix {output.outdir}/ \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterType BySJout \
                --outFilterMultimapNmax 20 \
                --alignIntronMax 1000000 \
                --alignMatesGapMax 1000000 \
                --outSAMattributes NH HI AS NM MD \
                --outSAMmultNmax 1 \
            > {log} 2>&1
        """

rule star_align_to_fastq:
    """
    Convert STAR output BAM to FASTQ files for remapping with minimap2.
    This is necessary because minimap2 does not support BAM input directly.
    """
    threads: 1
    input:
        bam="outputs/LCM-RNA-editing/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        qsort_bam="outputs/LCM-RNA-editing/minimap_remap/{sample}/star-mapped.qsort.bam",
        r1="outputs/LCM-RNA-editing/minimap_remap/{sample}/star-mapped.r1.fq",
        r2="outputs/LCM-RNA-editing/minimap_remap/{sample}/star-mapped.r2.fq"
    log:
        "logs/LCM-RNA-editing/star-align-to-fastq/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools sort \
            -n \
            -o {output.qsort_bam} \
            {input.bam} \
            > {log} 2>&1

        docker run {docker_mount} -u $(id -u) --rm staphb/bedtools:2.31.1 \
            bedtools bamtofastq \
                -i {output.qsort_bam} \
                -fq {output.r1} \
                -fq2 {output.r2} \
            >> {log} 2>&1
        """

rule star_align_to_1_fastq:
    """
    Convert STAR output BAM to FASTQ files for remapping with minimap2.
    This is necessary because minimap2 does not support BAM input directly.
    """
    threads: 1
    input:
        bam="outputs/LCM-RNA-editing/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        fastq="outputs/LCM-RNA-editing/minimap_remap/{sample}/star-mapped.fq",
    log:
        "logs/LCM-RNA-editing/star-align-to-1-fastq/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools fastq \
            {input.bam} \
            1> {output.fastq} \
            2> {log}
        """

# References:
# https://lh3.github.io/2025/04/18/short-rna-seq-read-alignment-with-minimap2
# https://lh3.github.io/minimap2/minimap2.html
rule minimap_remap:
    threads: 24
    input:
        assembly="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa",
        fq="outputs/LCM-RNA-editing/minimap_remap/{sample}/star-mapped.fq",
    output:
        remapped_sam="outputs/LCM-RNA-editing/minimap_remap/{sample}/remapped.sam",
        remapped_sorted_bam="outputs/LCM-RNA-editing/minimap_remap/{sample}/remapped.sorted.bam",
        remaaped_sorted_bam_bai="outputs/LCM-RNA-editing/minimap_remap/{sample}/remapped.sorted.bam.bai"
    log:
        minimap2="logs/LCM-RNA-editing/minimap_remap/{sample}-minimap2.log",
        samtools_sort="logs/LCM-RNA-editing/minimap_remap/{sample}-samtools-sort.log",
        samtools_index="logs/LCM-RNA-editing/minimap_remap/{sample}-samtools-index.log"
    shell:
        """
        # triggers rerun
        docker run {docker_mount} -u $(id -u) --rm staphb/minimap2:2.29 \
            minimap2 \
                -ax splice:sr \
                -t {threads} \
                {input.assembly} \
                {input.fq} \
                1> {output.remapped_sam} \
                2> {log.minimap2}
        
        docker run {docker_mount} -u $(id -u) --rm biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools sort \
                -@ {threads} \
                -o {output.remapped_sorted_bam} \
                {output.remapped_sam} \
                1> {log.samtools_sort} \
                2> {log.samtools_sort}
        
        docker run {docker_mount} -u $(id -u) --rm biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools index \
                -@ {threads} \
                {output.remapped_sorted_bam} \
                1> {log.samtools_index} \
                2> {log.samtools_index}
        """

        
rule samtools_index_lcm:
    threads: 8
    input:
        "outputs/LCM-RNA-editing/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "outputs/LCM-RNA-editing/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/LCM-RNA-editing/samtools_index_lcm/{sample}.log"
    shell:
        """
        docker run {docker_mount} -u $(id -u) --rm biocontainers/samtools:v1.9-4-deb_cv1 \
            samtools index \
                -@ {threads} \
                {input} \
            > {log} 2>&1
        """
# ==== 2. Align DNA reads to genome ====
# NOTE: Skipped (see above for explanation)

# ==== 3. REDItools on RNA (on scRNA-variable sites) ====
rule reditools_rna:
    threads: 16
    input:
        bam="outputs/LCM-RNA-editing/minimap_remap/{sample}/remapped.sorted.bam",
        bai="outputs/LCM-RNA-editing/minimap_remap/{sample}/remapped.sorted.bam.bai",
        assembly="references/Ptr/assembly/Ptrichocarpa_533_v4.0.fa"
        # bed_file="outputs/VariantCalling/vawk/ptr_tenx_batch1.snv.loci.bed"
    output:
        "outputs/LCM-RNA-editing/reditools-rna/{sample}.tsv"
    log:
        "logs/LCM-RNA-editing/reditools-rna/{sample}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) chiaenu/reditools3:3.4 \
            python -m reditools analyze \
                    --reference {input.assembly} \
                    --output-file {output} \
                    --strand 1 \
                    --verbose \
                    --variants all \
                    --threads {threads} \
                    {input.bam} \
            > {log} 2>&1
        """

# Subset the RNA editing sites to only those that are variable in the scRNA-seq data
rule bgzip_tabix:
    threads: 1
    input:
        "outputs/LCM-RNA-editing/reditools-rna/{sample}.tsv"
    output:
        gz="outputs/LCM-RNA-editing/bgzip-tabix/{sample}.tsv.gz",
        tbi="outputs/LCM-RNA-editing/bgzip-tabix/{sample}.tsv.gz.tbi"
    log:
        "logs/LCM-RNA-editing/bgzip-tabix/{sample}-subset.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) staphb/htslib:1.21 \
            scripts/bgzip-tabix.sh \
                {input} \
                {output.gz} \
                {output.tbi} \
                > {log} 2>&1
        """

rule subset_reditable:
    threads: 1
    input:
        reditable="outputs/LCM-RNA-editing/bgzip-tabix/{sample}.tsv.gz",
        bed_file="outputs/VariantCalling/vawk/ptr_tenx_batch1.snv.loci.bed"
    output:
        chrom_pos=temp("outputs/LCM-RNA-editing/subset-reditable/{sample}-chrom-pos.tsv"),
        reditable_subset="outputs/LCM-RNA-editing/subset-reditable/{sample}-subset.tsv"
    log:
        "logs/LCM-RNA-editing/subset-reditable/{sample}-subset.log"
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} {{print $1, $2}}' {input.bed_file} 1> {output.chrom_pos} 2> {log}
        docker run {docker_mount} --rm -u $(id -u) staphb/htslib:1.21 \
            tabix -R {output.chrom_pos} {input.reditable} 1> {output.reditable_subset} 2>> {log}
        """
# ==== 4. REDItools on DNA (on scRNA-variable sites) ====
# NOTE: Skipped (see above for explanation)

# ==== 5. REDItools annotate RNA with DNA ====
# NOTE: Skipped (see above for explanation)