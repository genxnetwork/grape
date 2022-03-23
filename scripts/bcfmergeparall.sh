#!/bin/bash
function vcf_chromosomes() {
    # Fetch contigs from a vcf file.
    bcftools view -h $1 | \
    grep 'contig' | \
    egrep -o "ID=([^,]+)" | \
    sed 's/ID=//g' | \
    tr '\n' ' '
}

function parallel_bcftools_merge() {
    file_set=`echo $@ | egrep -o '(\-l|\-\-file-list)(=|[ ]+)[^ ]+' | tr '=' ' ' | cut -f 2 -d ' '`
    if [ -n "${file_set}" ]
        then
            find_vcf=`cat ${file_set} | head -n 1`
        else
            find_vcf=`echo $@ | tr '\t' '\n' | egrep -o '[^ ]+.vcf.gz' | awk 'NR == 1 { print }' - `
    fi
    contigs=`vcf_chromosomes ${find_vcf}`
    current_dir=$(dirname ${find_vcf})
    hash_merge=`echo "$@" | md5sum | cut -c 1-5`
    output_prefix="${current_dir}/parallel_merge.${hash_merge}"

    parallel --gnu --workdir ${current_dir} \
    --env args -j ${PARALLEL_CORES} \
    'bcftools merge -r {1} -O u ' $@ ' > ' ${output_prefix}'.{1}.bcf' ::: ${contigs}

    order=`echo $contigs | tr ' ' '\n' | awk -v "prefix=${output_prefix}" '{ print prefix "." $0 ".bcf" }'`
    bcftools concat -O v ${order} | grep -v 'parallel_merge' | sed 's/##bcftools_mergeCommand=merge -r I -O u /##bcftools_mergeCommand=merge /g' | bcftools view -O u
    rm ${order}
}