#!/bin/bash
# @Author   : Yang Liu
# @FileName : unify_format_for_calls.sh
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome

# Generate unified read-level and/or site-level format of calls
# Usage:
# [prog]  <dsname>   <toolname>  <call-encode> <call-fn> <outd-dir> <num-processors> <step12> <chr-options>
#             1       2               3       4               5            6          7    8
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
