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

# ================= Custom functions =================
# query: query from list of dict by key-value pair
from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    return [x for x in d if x[k] == v][0]

# set_sample_names_param (used in CellRanger_count)
def set_sample_names_param(wildcards):
    sample = query(config["samples"], "name", wildcards.sample)
    if "sample_names" in sample.keys():
        return "--sample=%s" % sample["sample_names"]
    else:
        return ""

# ================= CellRanger count =================
rule CellRanger_count:
    threads: 8
    resources:
        mem_gb=100
    input:
        CellRanger="src/cellranger-9.0.1/cellranger",
        transcriptome=lambda wildcards: 
            "outputs/CellRanger/mkref/%s" % (
                query(config["samples"], "name", wildcards.sample)["species"],
            ),
        fastq_path=lambda wildcards: 
            query(config["samples"], "name", wildcards.sample)["fastq_path"]
    output:
        outdir=directory("outputs/CellRanger/count/{sample}"),
        out_matrix="outputs/CellRanger/count/{sample}/outs/filtered_feature_bc_matrix.h5",
        out_bam="outputs/CellRanger/count/{sample}/outs/possorted_genome_bam.bam",
        out_barcodes="outputs/CellRanger/count/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        tmpdir=temp(directory("{sample}_count_tmp"))
    params:
        create_bam="true",
        sample_names=lambda wildcards:
            "--sample=%s" % query(config["samples"], "name", wildcards.sample)["sample_names"] \
                if query(config["samples"], "name", wildcards.sample).get("sample_names", False) \
                else "",
        fastq_path=lambda wildcards: 
            ",".join(query(config["samples"], "name", wildcards.sample)["fastq_path"])
    log: 
        "logs/CellRanger/count/{sample}.log"
    shell:
        """
        {input.CellRanger} count \
            --id={output.tmpdir} \
            --transcriptome={input.transcriptome} \
            --fastqs={params.fastq_path} \
            --create-bam={params.create_bam} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            {params.sample_names} \
        2> {log} \
        1> {log} \
        && mkdir -p {output.outdir} \
        && mv {output.tmpdir}/* {output.outdir}
        """ 

# ================= CellRanger mkref =================
rule CellRanger_mkref:
    threads: 16
    resources:
        mem_gb=64
    input:
        CellRanger="src/cellranger-9.0.1/cellranger",
        assembly=lambda wildcards: 
            query(config["references"], "species", wildcards.species)["assembly"],
        # gtf=lambda wildcards: 
        #     query(config["references"], "species", wildcards.species)["annotation-gtf"],
        gtf="outputs/CellRanger/mkgtf/{species}.gtf-filtered"
    output:
        directory("outputs/CellRanger/mkref/{species}")
    log:
        "logs/CellRanger/mkref/{species}.log"
    shell:
        """
        {input.CellRanger} mkref \
            --genome=$(basename {output}) \
            --fasta={input.assembly} \
            --genes={input.gtf} \
            --nthreads={threads} \
            --memgb={resources.mem_gb} \
        2> {log} \
        1> {log} \
        && mv $(basename {output}) $(dirname {output})
        """

# ================= CellRanger mkgtf =================
rule CellRanger_mkgtf:
    input:
        CellRanger="src/cellranger-9.0.1/cellranger",
        gtf=lambda wildcards: 
            query(config["references"], "species", wildcards.species)["annotation-gtf"],
    output:
        "outputs/CellRanger/mkgtf/{species}.gtf-filtered"
    log:
        "logs/CellRanger/mkgtf/{species}.log"
    shell:
        """
        {input.CellRanger} mkgtf \
            {input.gtf} {output} \
            --attribute=gene_biotype:protein_coding \
        2> {log} \
        1> {log}             
        """

# ================= GTF/GFF(GFF3) conversion =================
# convert gff3 to gtf
rule gffread_gff32gtf:
    input:
        gff="{any}.gff3"
    output:
        gtf="{any}.gtf"
    log:
        "{any}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            quay.io/biocontainers/gffread:0.12.1--h8b12597_0 \
                gffread {input.gff} -T -o {output.gtf} \
        2> {log} \
        1> {log}
        """

# convert gff to gtf
rule gffread_gff2gtf:
    input:
        gff="{any}.gff"
    output:
        gtf="{any}.gtf"
    log:
        "{any}.log"
    shell:
        """
        docker run \
            {docker_mount} \
            -u $(id -u) \
            --rm \
            quay.io/biocontainers/gffread:0.12.1--h8b12597_0 \
                gffread {input.gff} -T -o {output.gtf} \
        2> {log} \
        1> {log}
        """

# ================= MARS results to h5  =================
rule mars2hd5:
    conda: "envs/bio-pyps.yaml"
    threads: 1
    input:
        umi_dir=lambda wildcards: 
            query(config["samples"], "name", wildcards.sample)["umi_dir"],
        cds_list=lambda wildcards: 
            query(config["samples"], "name", wildcards.sample)["cds_list"]
    output:
        h5="outputs/mars2hd5/{sample}.h5"
    params:
        umi_criteria=100
    log:
        "logs/mars2hd5/{sample}.log"
    shell:
        """
        python scripts/mars2hd5.py \
            --umi_dir {input.umi_dir} \
            --cds {input.cds_list} \
            --output_hd5 {output.h5} \
            --umi_criteria {params.umi_criteria} \
        2> {log} \
        1> {log}
        """

# ================= get CellRanger =================
rule get_CellRanger:
    output:
        executable="src/cellranger-9.0.1/cellranger"
    log:
        "logs/get_CellRanger.log"
    shell:
        """
        cd src/
        curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1743615920&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=QYxJGDFIH9Uo5YNweuGgg~CB4AsoQv1CzCgq8BXCxSvczdwG0EQ3BJtPlIcxnWYresqxYqGfX1Lj4ei9lGUZjqSvrdiRwg8qRKR-~hIm0ukCEA9r8HmJJBjUQ3EKBePzfsmKJSLyUOL8zO8SuRlg9anUaeEfQcAGHYLR~r6GE2Hjp5mkCjtTUDHxGThekZgeMgkiHdx6uj8breRiiMFspEMSp~3MdrHzldr9DpySjY18TNsJe0JrgiCNg0Ool2O7qnrNCzygcYQChahG1hdDxRb7mBS-qsAHQ9Mo2inAFXrl7YlypmSgtZ2cZfgko5hnsNTV-Bo2a9yeq1~VNzFKKg__" \
            && tar -xzf cellranger-9.0.1.tar.gz \
            && rm cellranger-9.0.1.tar.gz \
        """