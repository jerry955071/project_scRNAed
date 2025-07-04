#!/bin/bash
possorted_genome_bam=$1
umi_reads_bam=$2


# Filter BAM file for records with 'xf:i:25' and write that to body_filtered_sam. body_filtered_sam is on the sam format and will not include the BAM header.
samtools view $possorted_genome_bam | LC_ALL=C grep "xf:i:25" > ${umi_reads_bam}.body_filtered_sam

# Extract the BAM header and write to header_filted_sam
samtools view -H $possorted_genome_bam > ${umi_reads_bam}.header_filted_sam

# Combine the header and extracted records
cat ${umi_reads_bam}.header_filted_sam ${umi_reads_bam}.body_filtered_sam > ${umi_reads_bam}_filterd.sam

# Convert SAM to BAM
samtools view -b ${umi_reads_bam}_filterd.sam > $umi_reads_bam

# Clean up temp files
rm ${umi_reads_bam}.body_filtered_sam ${umi_reads_bam}.header_filted_sam ${umi_reads_bam}_filterd.sam
