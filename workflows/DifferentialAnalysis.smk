# get pseudo-time cell ordering
rule slingshot:
    threads: 1
    resources:
        mem_mb=10000,
        runtime=100
    container: "docker://chiaenu/rmd-slingshot:2.4.0"
    input:
        path_cca_embedding="outputs/notebooks-new/Seurat_integration/cca_embedding.csv",
        path_cell_type_csv="outputs/notebooks-new/Cell_type_annotation/cell_type_annotation.csv"
    output:
        path_pseudo_time="outputs/notebooks-new/Slingshot/pseudotime.csv",
        path_lineage_weight="outputs/notebooks-new/Slingshot/weight.csv",
        path_outdir="outputs/notebooks-new/Slingshot"
    log:
        "logs/DifferentialAnalysis/slingshot.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            'notebooks-new/Slingshot.Rmd',
            output_dir = '{output.path_outdir}',
            output_file = 'slingshot.html',
            knit_root_dir = '~',
            params = list(
                path_cca_embedding = '{input.path_cca_embedding}',
                path_cell_type_csv='{input.path_cell_type_csv}',
                path_pseudo_time='{output.path_pseudo_time}',
                path_lineage_weight='{output.path_lineage_weight}',
                path_outdir='{output.path_outdir}'
            )
        )" 2>&1 > {log}
        """


# perform beta-binomial regression for differential analysis
# samples and corresponding metadata were included in {input.cell_labels}
rule beta_binomial_regression:
    input:
        alt_mtx_template="outputs/MetaCell/aggregate_metacell_vartrix/{{sample}}/alt.mtx",
        ref_mtx_template="outputs/MetaCell/aggregate_metacell_vartrix/{{sample}}/ref.mtx",
        row_idx_template="outputs/MetaCell/aggregate_metacell_vartrix/{{sample}}/row_idx.txt",
        col_idx_template="outputs/MetaCell/aggregate_metacell_vartrix/{{sample}}/col_idx.txt",
        cell_labels="outputs/MetaCell/metacell_label_annotation/metacell_annotation.csv",
        pseudo_time="outputs/notebooks-new/Slingshot/pseudotime.csv",
        lineage_weight="outputs/notebooks-new/Slingshot/weight.csv"