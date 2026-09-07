# perform binomial regression for differential analysis
# samples and corresponding metadata were included in {input.path_metacell_metadata}
def call_mc_label_files(wildcards):
    import pandas as pd
    df = pd.read_csv("outputs/MetaCell/metacell_summary/metacell_metadata.csv")
    samples = df["sample"].unique()
    return [f"outputs/MetaCell/mcrigor/{s}/metacell.csv" for s in samples]

rule binomial_regression:
    container: "docker://chiaenu/rstudio-rstat:4.4.3"
    threads: 4
    resources:
        mem_mb=100000,
        runtime=40
    input:
        path_mtx_dirs=call_mc_label_files,
        path_metacell_metadata="outputs/MetaCell/metacell_summary/metacell_metadata.csv",
        path_seurat_obj="outputs/notebooks-new/Seurat_integration/seurat_integrated.rds",
        path_hom_loci="outputs/VariantCalling-DNA/gatk_joint/ptr/hom_ref.vcf",
        path_variant_annotation="outputs/VariantAnnotation/variant_annotation/ptr/variant_annotation.txt",
        path_variant_location="outputs/VariantAnnotation/variant_annotation/ptr/variant_location.txt"
    output:    
        path_outdir=directory("outputs/DEA_Regression/binomial_regression/")
    params:
        path_mtx_dirs="outputs/MetaCell/aggregate_metacell_vartrix/%s",
        DEBUG="FALSE"
    log:
        "logs/DifferentialAnalysis/binomial_regression.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            input = 'notebooks-new/DEA_Regression.Rmd',
            output_dir = '{output.path_outdir}',
            output_file = 'binomial_regression.html',
            knit_root_dir = '~',
            params = list(
                path_mtx_dirs='{params.path_mtx_dirs}',
                path_metacell_metadata='{input.path_metacell_metadata}',
                path_seurat_obj='{input.path_seurat_obj}',
                path_hom_loci='{input.path_hom_loci}',
                path_variant_annotation='{input.path_variant_annotation}',
                path_variant_location='{input.path_variant_location}',
                path_outdir='{output.path_outdir}',
                DEBUG='{params.DEBUG}'
            )
        )" 2>&1 > {log}
        """