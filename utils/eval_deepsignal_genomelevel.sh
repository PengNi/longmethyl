#!/bin/bash

set -x
bsbedmethyl=${1}
ds_site_tsv=${2}
result_data=${3}
result_log=${4}

# chroms="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

echo "#### eval deepsignal at genome level"

#python src/correlation_with_any.data.py \
#    --queryfile ${ds_site_tsv} \
#    --targetfile ${bsbedmethyl} \
#    --contig_names ${chroms} --wfile ${result_data} > ${result_log}
python src/correlation_with_any.data.py \
    --queryfile ${ds_site_tsv} \
    --targetfile ${bsbedmethyl} \
    --wfile ${result_data} > ${result_log}

echo "#### eval deepsignal at genome level, DONE"











