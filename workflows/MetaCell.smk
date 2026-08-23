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
                path_output = '{output.out_dir}',
                path_mc_label_template = 'outputs/MetaCell/metacell/{{sample}}/K={{K}}/alpha={{alpha}}/metacell.csv'
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

rule seurat_integration:
    threads: 1
    resources:
        mem_mb=30000,
        runtime=100
    container: "docker://chiaenu/rstudio-mcrigor:1.0"
    output:
        out_dir=directory("outputs/notebooks-new/Seurat_integration/"),
        path_seurat_rds="outputs/notebooks-new/Seurat_integration/seurat_integrated.rds",
        path_cca_embedding="outputs/notebooks-new/Seurat_integration/cca_embedding.csv"
    log:
        "logs/MetaCell/seurat_integration.log"
    shell:
        """
        # 
        Rscript -e "rmarkdown::render(
            'notebooks-new/Seurat_integration.Rmd',
            output_dir = '{output.out_dir}',
            output_file = 'seurat_integration.html',
            knit_root_dir = '~',
            params = list(
                path_seurat_rds = '{output.path_seurat_rds}',
                path_h5='outputs/CellRanger/count/{{sample}}/outs/filtered_feature_bc_matrix.h5',
                path_cca_embedding='{output.path_cca_embedding}',
                path_outdir='{output.out_dir}'
            )
        )" 2>&1 > {log}
        """

rule setup_reference_cell_label:
    input:
        cell_label_reference="references/ptr_tenx_batch1_rs17_curated.csv"
    output:
        cell_label_out="references/ptr_tenx_batch1_rs17_curated_named.csv"
    log:
        "logs/MetaCell/setup_reference_cell_label.log"
    run:
        import pandas as pd
        cluster_mapping = {
            1: "Vessel",
            2: "FuIP",
            3: "RO",
            4: "RP",
            5: "FuEP",
            6: "FuO",
            7: "Fiber",
            8: "Ray",
            9: "Outlier_9",
            10: "Outlier_10",
            11: "Outlier_11"
        }
        df = pd.read_csv(input.cell_label_reference)
        df["Label"] = df["Cluster"].map(cluster_mapping)
        df.to_csv(output.cell_label_out, index=False)

rule cell_type_annotation:
    threads: 1
    resources:
        mem_mb=30000,
        runtime=100
    container: "docker://chiaenu/rstudio-mcrigor:1.0"
    input:
        integrated_seurat_rds="outputs/notebooks-new/Seurat_integration/seurat_integrated.rds",
        cell_label_reference="references/ptr_tenx_batch1_rs17_curated_named.csv"
    output:
        out_dir=directory("outputs/notebooks-new/Cell_type_annotation/"),
        path_cell_type_csv="outputs/notebooks-new/Cell_type_annotation/cell_type_annotation.csv"
    params:
        ref_name="ptr_tenx_batch1"
    log:
        "logs/MetaCell/cell_type_annotation.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            'notebooks-new/Cell_type_annotation.Rmd',
            output_dir = '{output.out_dir}',
            output_file = 'cell_type_annotation.html',
            knit_root_dir = '~',
            params = list(
                path_cell_type_csv = '{output.path_cell_type_csv}',
                path_seurat_rds='{input.integrated_seurat_rds}',
                path_cell_label_reference='{input.cell_label_reference}',
                ref_name='{params.ref_name}',
                path_outdir='{output.out_dir}'
            )
        )" 2>&1 > {log}
        """

# rule metacell_label_annotation:
#     threads: 1
#     resources:
#         mem_mb=30000,
#         runtime=100
#     container: "docker://chiaenu/rstudio-mcrigor:1.0"
#     input:
#         mc_label_template="outputs/MetaCell/mcrigor/%s/metacell.csv",
#         cell_type_annotation="outputs/notebooks-new/Cell_type_annotation/cell_type_annotation.csv",
#         seurat_obj="outputs/notebooks-new/Seurat_integration/seurat_integrated.rds"
#     output:
#         out_dir=directory("outputs/MetaCell/metacell_label_annotation"),
#         path_metacell_annotated="outputs/MetaCell/metacell_label_annotation/metacell_annotation.csv"
#     log:
#         "logs/MetaCell/metacell_label_annotation.log"
#     shell:
#         """
#         Rscript -e "rmarkdown::render(
#             'notebooks-new/Metacell_label_annotation.Rmd',
#             output_dir = '{output.out_dir}',
#             output_file = 'metacell_label_annotation.html',
#             knit_root_dir = '~',
#             params = list(
#                 path_metacell_annotated = '{output.path_metacell_annotated}',
#                 path_mc_label_template='{input.mc_label_template}',
#                 path_cell_type_annotation='{input.cell_type_annotation}',
#                 path_seurat_obj='{input.seurat_obj}',
#                 path_outdir='{output.out_dir}'
#             )
#         )" 2>&1 > {log}
#         """

# get pseudo-time cell ordering
rule slingshot:
    threads: 1
    resources:
        mem_mb=10000,
        runtime=200
    container: "docker://chiaenu/rmd-slingshot:2.4.0"
    input:
        path_cca_embedding="outputs/notebooks-new/Seurat_integration/cca_embedding.csv",
        path_cell_type_csv="outputs/notebooks-new/Cell_type_annotation/cell_type_annotation.csv"
    output:
        path_pseudo_time="outputs/notebooks-new/Slingshot/pseudotime.csv",
        path_lineage_weight="outputs/notebooks-new/Slingshot/weight.csv",
        path_lineage_prob="outputs/notebooks-new/Slingshot/prob.csv",
        path_branch_id="outputs/notebooks-new/Slingshot/branch_id.csv",
        path_outdir=directory("outputs/notebooks-new/Slingshot")
    params:
        DEBUG="FALSE"
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
                path_lineage_prob='{output.path_lineage_prob}',
                path_branch_id='{output.path_branch_id}',
                path_outdir='{output.path_outdir}',
                DEBUG={params.DEBUG}
            )
        )" 2>&1 > {log}
        """


# summarize metacell pseudotime/cell labels for downstream differential analysis
def call_mc_label_files(wildcards):
    import pandas as pd
    df = pd.read_csv("outputs/notebooks-new/Cell_type_annotation/cell_type_annotation.csv")
    samples = df["orig.ident"].unique()
    return [f"outputs/MetaCell/mcrigor/{s}/metacell.csv" for s in samples]

rule metacell_summary:
    threads: 1
    resources:
        mem_mb=30000,
        runtime=200
    container: "docker://chiaenu/rstudio-rstat:4.4.3"
    input:
        path_mc_label=call_mc_label_files,
        path_cell_type_annotation="outputs/notebooks-new/Cell_type_annotation/cell_type_annotation.csv",
        path_pseudo_time="outputs/notebooks-new/Slingshot/pseudotime.csv",
        path_lineage_weight="outputs/notebooks-new/Slingshot/weight.csv",
        path_lineage_prob="outputs/notebooks-new/Slingshot/prob.csv",
        path_branch_id="outputs/notebooks-new/Slingshot/branch_id.csv",
        path_seurat_obj="outputs/notebooks-new/Seurat_integration/seurat_integrated.rds"
    output:
        path_outdir=directory("outputs/MetaCell/metacell_summary"),
        path_metacell_metadata="outputs/MetaCell/metacell_summary/metacell_metadata.csv"
    params:
        DEBUG="FALSE",
        path_mc_label_template="outputs/MetaCell/mcrigor/%s/metacell.csv"
    log:
        "logs/DifferentialAnalysis/metacell_summary.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            'notebooks-new/Metacell_trajectory_weight.Rmd',
            output_dir = '{output.path_outdir}',
            output_file = 'metacell_summary.html',
            knit_root_dir = '~',
            params = list(
                path_mc_label_template='{params.path_mc_label_template}',
                path_cell_type_annotation='{input.path_cell_type_annotation}',
                path_pseudo_time='{input.path_pseudo_time}',
                path_lineage_weight='{input.path_lineage_weight}',
                path_lineage_prob='{input.path_lineage_prob}',
                path_branch_id='{input.path_branch_id}',
                path_seurat_obj='{input.path_seurat_obj}',
                path_outdir='{output.path_outdir}',
                path_metacell_metadata='{output.path_metacell_metadata}',
                DEBUG='{params.DEBUG}'
            )
        )" 2>&1 > {log}
        """