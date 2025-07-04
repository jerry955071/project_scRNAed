# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode
# Snakemake setups
configfile: "configs/config.json"
wildcard_constraints:
    sample="|".join([i["name"] for i in config["samples-dna"]]),
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
# get_assembly: get assembly file path by species
def get_assembly_by_sample_name(sample_name:str) -> str:
    species = query(config["samples-dna"], "name", sample_name)["species"]
    assembly = query(config["references"], "species", species)["assembly"]
    return assembly

# === 1. Add read group ===
rule gatk_addreplacerg:
    input:
        lambda wildcards: query(config["samples-dna"], "name", wildcards.sample)["bam_path"]
    output:
        "outputs/VariantCalling-DNA/gatk_addreplacerg/{sample}.sorted.addrg.bam"
    log:
        "logs/VariantCalling-DNA/gatk_addreplacerg/{sample}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk AddOrReplaceReadGroups \
                -I {input} \
                -O {output} \
                -RGID {wildcards.sample} \
                -RGSM {wildcards.sample} \
                -RGLB lib1 \
                -RGPL illumina \
                -RGPU unit1 \
            > {log} 2>&1
        """

# === 2. Mark duplicates ===
rule gatk_markduplicates:
    input:
        "outputs/VariantCalling-DNA/gatk_addreplacerg/{sample}.sorted.addrg.bam"
    output:
        bam="outputs/VariantCalling-DNA/gatk_markduplicates/{sample}.sorted.addrg.markdup.bam",
        metrics="outputs/VariantCalling-DNA/gatk_markduplicates/{sample}.metrics.txt"
    log:
        "logs/VariantCalling-DNA/gatk_markduplicates/{sample}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk --java-options "-Xmx4g" MarkDuplicates \
                -I {input} \
                -O {output.bam} \
                -M {output.metrics} \
                --CREATE_INDEX true \
            > {log} 2>&1
        """

# === 3. Create sequence dictionary ===
rule gatk_create_seq_dict:
    input:
        "{assembly}.fa"
    output:
        "{assembly}.dict"
    log:
        "logs/VariantCalling-DNA/gatk_create_seq_dict/{assembly}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk CreateSequenceDictionary \
                -R {input} \
                -O {output} \
            > {log} 2>&1
        """

# === 3. HaplotypeCaller (GVCF mode) ===
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531532-HaplotypeCaller-Reference-Confidence-Model-GVCF-mode
# TODO: target loci will be merged per species
target_loci = {
    "ptr": "outputs/VariantCalling/ptr/snv.loci.list",
    "egr": "outputs/VariantCalling/egr/snv.loci.list"
}
rule gatk_haplotypecaller:
    input:
        bam="outputs/VariantCalling-DNA/gatk_markduplicates/{sample}.sorted.addrg.markdup.bam",
        reference=lambda wildcards: get_assembly_by_sample_name(wildcards.sample),
        reference_dict=lambda wildcards: get_assembly_by_sample_name(wildcards.sample).replace(".fa", ".dict"),
        intervals=lambda wildcards: target_loci[query(config["samples-dna"], "name", wildcards.sample)["species"]]
    output:
        gvcf="outputs/VariantCalling-DNA/gatk_gvcf/{sample}.g.vcf.gz"
    log:
        "logs/VariantCalling-DNA/gatk_gvcf/{sample}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk --java-options "-Xmx4g" HaplotypeCaller \
                -R {input.reference} \
                -I {input.bam} \
                -O {output.gvcf} \
                -L {input.intervals} \
                -ERC GVCF \
            > {log} 2>&1
        """

# === 4. GenomicsDBImport (combine GVCFs) ===
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/30332006200603-GenomicsDBImport
# TODO: try remove -L; remove rule 'create_interval_list' if is faster
rule gatk_genomicsdbimport:
    threads: 32
    input:
        gvcfs=lambda wildcards: [
            "outputs/VariantCalling-DNA/gatk_gvcf/%s.g.vcf.gz" % query_all(
                d=config["samples-dna"],
                k="species",
                v=wildcards.species,
                k_out="name"
            )
        ],
        intervals="outputs/VariantCalling-DNA/{species}/intervals.list"
    output:
        db_dir=directory("outputs/VariantCalling-DNA/gatk_joint/{species}/db")
    log:
        "logs/VariantCalling-DNA/gatk_joint/{species}/genomicsdb.log"
    params:
        variant_flags=lambda wildcards, input: " ".join(f"--variant {f}" for f in input.gvcfs),
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk --java-options "-Xmx4g" GenomicsDBImport \
                --genomicsdb-workspace-path {output.db_dir} \
                --batch-size 50 \
                -L {input.intervals} \
                --reader-threads {threads} \
                {params.variant_flags} \
            > {log} 2>&1
        """

rule create_intervals_list:
    input:
        lambda wildcards: query(config["references"], "species", wildcards.species)["assembly"] + ".fai"
    output:
        "outputs/VariantCalling-DNA/{species}/intervals.list"
    log:
        "logs/VariantCalling-DNA/gatk_joint/create_intervals_list/{species}.log"
    shell:
        """
        cat {input} | awk '{{print $1}}' > {output}
        """

# === 5. GenotypeGVCFs (joint genotyping) ===
rule gatk_genotypegvcfs:
    input:
        db_dir="outputs/VariantCalling-DNA/gatk_joint/{species}/db",
        reference=lambda wildcards: query(config["references"], "species", wildcards.species)["assembly"],
        intervals=lambda wildcards: {
            "ptr": "outputs/VariantCalling/vawk/ptr_tenx_batch1.snv.loci.list",
            "egr": "outputs/VariantCalling/vawk/egr_tenx_batch1.snv.loci.list"
        }[wildcards.species]
    output:
        vcf="outputs/VariantCalling-DNA/gatk_joint/{species}.raw.vcf.gz"
    log:
        "logs/VariantCalling-DNA/gatk_joint/genotypegvcfs.{species}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk --java-options "-Xmx4g" GenotypeGVCFs \
                -R {input.reference} \
                -V gendb://{input.db_dir} \
                -O {output.vcf} \
                -L {input.intervals} \
                --include-non-variant-sites true \
            > {log} 2>&1
        """

# === 6. Filter raw variants ===
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
rule gatk_variant_filtration:
    input:
        vcf="outputs/VariantCalling-DNA/gatk_joint/{species}.raw.vcf.gz",
        reference=lambda wildcards: query(config["references"], "species", wildcards.species)["assembly"],
    output:
        vcf="outputs/VariantCalling-DNA/gatk_joint/{species}.filtered.vcf.gz"
    log:
        "logs/VariantCalling-DNA/gatk_joint/variant_filter.{species}.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            gatk --java-options "-Xmx4g" VariantFiltration \
                -R {input.reference} \
                -V {input.vcf} \
                -O {output.vcf} \
                --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
                --filter-name "basic_snp_filter" \
            > {log} 2>&1
        """

# Understanding GATK vcf output: https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format
# QUAL and GQ: https://gatk.broadinstitute.org/hc/en-us/articles/360035531392-Difference-between-QUAL-and-GQ-annotations-in-germline-variant-calling
# === 7. Variants to table ===
# Reference genotype confidence model: https://gatk.broadinstitute.org/hc/en-us/articles/360035531532-HaplotypeCaller-Reference-Confidence-Model-GVCF-mode
# rule gatk_varianstotable:
#     input:
#         raw_vcf="outputs/VariantCalling-DNA/gatk_joint/Ptr.raw.vcf.gz",
#     output:
#         tsv="outputs/VariantCalling-DNA/gatk_variantstotable/Ptr.raw.tsv"
#     log:
#         "logs/VariantCalling-DNA/gatk_variantstotable/Ptr.log"
#     shell:
#         """
#         docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
#             gatk VariantsToTable \
#                 -V {input.raw_vcf} \
#                 -F CHROM \
#                 -F POS \
#                 -F TYPE \
#                 -F HET \
#                 -F HOM-REF \
#                 -F HOM-VAR \
#                 -F NSAMPLES \
#                 -F NCALLED \
#                 -GF RGQ \
#                 -GF GQ \
#                 -GF DP \
#                 -O {output.tsv} \
#             > {log} 2>&1
#         """

# ==== 8. Find no-variant sites ====
# TODO: add output snp.vcf
rule categorize_sites:
    input:
        ref=lambda wildcards: query(config["references"], "species", wildcards.species)["assembly"],
        vcf="outputs/VariantCalling-DNA/gatk_joint/{species}.raw.vcf.gz",
    output:
        hom_alt="outputs/VariantCalling-DNA/gatk_joint/{species}/hom_alt.vcf",
        hom_ref="outputs/VariantCalling-DNA/gatk_joint/{species}/hom_ref.vcf",
        het="outputs/VariantCalling-DNA/gatk_joint/{species}/het.vcf",
        snp="outputs/VariantCalling-DNA/gatk_joint/{species}/snp.vcf"
    log:
        "logs/VariantCalling-DNA/gatk_joint/{species}/categorize_sites.log"
    shell:
        """
        docker run {docker_mount} --rm -u $(id -u) broadinstitute/gatk:4.6.1.0 \
            scripts/categorize-sites.sh \
                {input.ref} \
                {input.vcf} \
                {output.hom_alt} \
                {output.hom_ref} \
                {output.het} \
                {output.snp} \
            > {log} 2>&1
        """
