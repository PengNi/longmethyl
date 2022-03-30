#!/bin/bash

set -x
dsname=${1}
toolname=${2}
callfn=${3}
sort=${4:-'false'}

# if [[ "${chr_options}" != "" ]] ; then
# 	chr_options="--chrSet ${chr_options}"
# fi

echo "### call freq"
if [ ${sort} == true ] ; then
    PYTHONPATH=src python src/call_modification_frequency.py \
        --input_path ${callfn} \
        --result_file ${dsname}_${toolname}_sitemods_freq.bed \
        --bed --sort > ${dsname}_${toolname}_sitemods_freq.bed.log 2>&1
else
    PYTHONPATH=src python src/call_modification_frequency.py \
        --input_path ${callfn} \
        --result_file ${dsname}_${toolname}_sitemods_freq.bed \
        --bed > ${dsname}_${toolname}_sitemods_freq.bed.log 2>&1
fi
gzip -f ${dsname}_${toolname}_sitemods_freq.bed > ${dsname}_${toolname}_sitemods_freq.bed.gz
rm ${dsname}_${toolname}_sitemods_freq.bed
    
echo "### call freq done"
