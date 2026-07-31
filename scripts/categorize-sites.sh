#!/bin/bash
# Run with docker image: broadinstitute/gatk:4.6.1.0
# Usage: ./categorize-sites.sh <ref.fa> <raw.vcf> <homozygous-alt.vcf> <homozygous-ref.vcf> <heterozygous.vcf> <snp.vcf>
ref_fa=$1
raw_vcf=$2
hom_alt_vcf=$3
hom_ref_vcf=$4
het_vcf=$5
snp_vcf=$6

# high-confidence homozygous alternate sites
cat $raw_vcf \
    | bcftools view -i 'F_PASS(FMT/GT="RA" & FMT/DP>10 & FMT/GQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="RR" & FMT/DP>10 & FMT/RGQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="AA" & FMT/DP>10 & FMT/GQ>20)>0.2' - \
    > $hom_alt_vcf

# high-confidence homozygous reference sites
cat $raw_vcf \
    | bcftools view -i 'F_PASS(FMT/GT="RA" & FMT/DP>10 & FMT/GQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="AA" & FMT/DP>10 & FMT/GQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="RR" & FMT/DP>10 & FMT/RGQ>20)>0.2' - \
    > $hom_ref_vcf

# low+high-confidence homozygous reference sites
cat $raw_vcf \
    | bcftools view -i 'F_PASS(FMT/GT="RA" & FMT/DP>10 & FMT/GQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="AA" & FMT/DP>10 & FMT/GQ>20)=0' - \
    > outputs/notebooks/UMAP-show-scSNV/genotyping/low_high_confidence_hom_ref.vcf


# high-confidence heterozygous no-variant sites
cat $raw_vcf \
    | bcftools view -i 'F_PASS(FMT/GT="AA" & FMT/DP>10 & FMT/GQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="RR" & FMT/DP>10 & FMT/RGQ>20)=0' - \
    | bcftools view -i 'F_PASS(FMT/GT="RA" & FMT/DP>10 & FMT/GQ>20)>0.2' - \
    > $het_vcf

# hard filtering for high-confidence snps 
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
gatk SelectVariants \
    -R $ref_fa \
    -V $raw_vcf \
    --select-type-to-include SNP \
    -select "QD > 2.0" \
    -select "FS < 60.0" \
    -select "SOR < 3.0" \
    -select "MQ > 40.0" \
    -select "MQRankSum > -12.5" \
    -select "ReadPosRankSum > -8.0" \
    -O $snp_vcf