#/bin/bash
set -euxo pipefail

# iterover combinations of sample & MetaCell parameters
sname=(
  ptr_tenx_tso2 ptr_tenx_tso3 ptr_tenx_tso4 ptr_tenx_tst2
  ptr_tenx_tst3 ptr_tenx_tst4 ptr_tenx_tsv2 ptr_tenx_tsv3 
  ptr_tenx_tsv4 ptr_tenx_tsv5 ptr_tenx_batch2
)
K=(20 50 100)
alpha=(0.5 1 2)

for s in "${sname[@]}"; do
    for k in "${K[@]}"; do
        for a in "${alpha[@]}"; do
            echo "Combination: $s - (K, alpha) = ($k, $a)"
            cmd="rmarkdown::render(
                'notebooks-new/MetaCell.Rmd',
                output_file = '~/local/outputs/notebooks-new/MetaCell/$sname.html',
                knit_root_dir='~/local',
                params = list(
                    sample_name = '$s',
                    K = $k,
                    alpha = $a
                )
              )
            "
            Rscript -e "${cmd}" > ~/local/outputs/notebooks-new/MetaCell/${s}_${k}_${a}.log 2>&1 &
            wait
        done
    done
done



echo "all done"