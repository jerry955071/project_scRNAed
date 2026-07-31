rule metacell:
    threads: 8
    resources:
        mem_mb=32000,
        runtime=100
    container: "docker://chiaenu/rstudio-metacell:0.3.7"
    input:
        matrix_base_dir="outputs/CellRanger/count/{sample}/outs/filtered_feature_bc_matrix"
    output:
        mc_label="outputs/MetaCell/metacell/{sample}/K={K}/alpha={alpha}/metacell.csv",
        main=directory("outputs/MetaCell/metacell/{sample}/K={K}/alpha={alpha}"),
        mcdb=directory("outputs/MetaCell/metacell/{sample}/K={K}/alpha={alpha}/mcdb"),
        figs=directory("outputs/MetaCell/metacell/{sample}/K={K}/alpha={alpha}/figs")
    params:
        s="{sample}",
        K="{K}",
        alpha="{alpha}"
    log:
        "logs/MetaCell/metacell/{sample}/K={K}/alpha={alpha}/metacell.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            'notebooks-new/MetaCell.Rmd',
            output_dir = '{output.main}',
            output_file = 'metacell.html',
            knit_root_dir = '~',
            params = list(
                sample_name = '{params.s}',
                K = {params.K},
                alpha = {params.alpha},
                path_mcdb = '{output.mcdb}',
                path_figs = '{output.figs}',
                path_10x = '{input.matrix_base_dir}',
                path_output = '{output.main}',
                path_mc_label = '{output.mc_label}'
            )
        )" 2>&1 > {log}
        """

rule mcrigor:
    threads: 8
    resources:
        mem_mb=32000,
        runtime=100
    container: "docker://chiaenu/rstudio-mcrigor:1.0"
    output:
        out_dir=directory("outputs/MetaCell/mcrigor/{sample}"),
        mc_label="outputs/MetaCell/mcrigor/{sample}/metacell.csv",
        mc_res="outputs/MetaCell/mcrigor/{sample}/optimize_res.rds"
    params:
        sample_name="{sample}",
        Ks="30 40 50 60 70",
        alphas="0.6 0.8 1 1.2"
    log:
        "logs/MetaCell/mcrigor/{sample}.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            'notebooks-new/McRigor.Rmd',
            output_dir = '{output.out_dir}',
            output_file = 'mcRigor.html',
            knit_root_dir = '~',
            params = list(
                sample_name = '{params.sample_name}',
                Ks = '{params.Ks}',
                alphas = '{params.alphas}',
                path_output = '{output.out_dir}'
            )
        )" 2>&1 > {log}
        """

rule aggregate_metacell_vartrix:
    threads: 1
    resources:
        mem_mb=20000,
        runtime=100
    container: "docker://chiaenu/rstudio-mcrigor:1.0" # any docker image with package 'Matrix' installed
    input:
        mc_label="outputs/MetaCell/mcrigor/{sample}/metacell.csv",
        alt_mtx="outputs/VariantCalling/vartrix/{sample}/alt.mtx",
        ref_mtx="outputs/VariantCalling/vartrix/{sample}/ref.mtx",
        row_idx="outputs/VariantCalling/vawk/{sample}.snv.loci.txt",
        col_idx="outputs/Remapping/renamer/{sample}/barcodes.tsv"
    output:
        alt_mtx="outputs/MetaCell/aggregate_metacell_vartrix/{sample}/alt.mtx",
        ref_mtx="outputs/MetaCell/aggregate_metacell_vartrix/{sample}/ref.mtx",
        row_idx="outputs/MetaCell/aggregate_metacell_vartrix/{sample}/row_idx.txt",
        col_idx="outputs/MetaCell/aggregate_metacell_vartrix/{sample}/col_idx.txt"
    log:
        "logs/MetaCell/aggregate_metacell_vartrix/{sample}.log"
    script:
        "scripts/aggregate_metacell_vartrix.R"
