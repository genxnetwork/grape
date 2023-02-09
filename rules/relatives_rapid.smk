import os
import sys


rule index_and_split:
    input: 
        vcf = "preprocessed/data.vcf.gz"
    output: 
        vcf = "vcf/for_rapid_{chrom}.vcf.gz"
    conda:
        "bcftools"
    log:
        "logs/vcf/index_and_split-{chrom}.log"
    benchmark:
        "benchmarks/vcf/index_and_split-{chrom}.txt"
    shell:
        """
        bcftools filter {input} -r {wildcards.chrom} -O z -o vcf/for_rapid_{wildcards.chrom}.vcf.gz | tee {log}
        """


rule estimate_rapid_params:
    input: 
        vcf = rules.index_and_split.output['vcf']
    conda:
        "rapid_params"
    log:
        "logs/rapid/estimate_rapid_params-{chrom}.log"
    output:
        rapid_params = "rapid/params_{chrom}"
    params:
        error_rate = config['rapid_error_rate'],
        min_snps = config['rapid_min_snp'],
        num_runs = config['rapid_num_runs'],
        num_success = config['rapid_num_success']
    script:
        "../scripts/estimate_rapid_params.py"


rule interpolate:
    input:
        vcf = rules.index_and_split.output['vcf'],
        cmmap=CMMAP
    output: "cm/chr{chrom}.cm.g"
    conda:
        "interpolation"
    log:
        "logs/cm/interpolate-{chrom}.log"
    script:
        "../scripts/interpolate.py"

'''
rule erase_dosages:
    input:
        vcf=rules.index_and_split.output[0]
    output:
        vcf='vcf/erased_chr{chrom}.vcf.gz'
    params:
        vcf='vcf/erased_chr{chrom}'
    conda:
        'bcftools'
    shell:
        "bcftools annotate -x 'fmt' {input.vcf} -O z -o {output.vcf}"
'''

rule rapid:
    input:
        vcf = rules.index_and_split.output['vcf'],
        g=rules.interpolate.output
    params:
        num_runs = config['rapid_num_runs'],
        num_success = config['rapid_num_success'],
        min_cm_len=3.0,
        window_size=3,
        output_folder='rapid'
    output:
        seg = "rapid/results_{chrom}.max.gz"
    conda:
        "rapid"
    shell:
        """
            rapidibd -i {input.vcf} -g {input.g} -d {params.min_cm_len} -o {params.output_folder}/chr{wildcards.chrom} \
            -w {params.window_size} -r {params.num_runs} -s {params.num_success}
        """


rule merge_rapid_segments:
    input:
        seg = expand('rapid/chr{chrom}/results.max.gz', chrom=CHROMOSOMES)
    output:
        seg = "rapid/results.max"
    shell:
        """
            for file in {input}
            do
                zcat ${{file}} >> {output}
            done
        """

rule ersa:
    input:
        ibd = rules.merge_rapid_segments.output['seg']
    output:
        "ersa/relatives.tsv"
    conda:
        "ersa"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    params:
        l = config['zero_seg_count'],
        th = config['zero_seg_len'],
        a = config['alpha'],
        t = config['rapid_seg_len'],
        r = '--nomask ' + '-r ' + str(config['ersa_r']) if config.get('weight_mask') else ''
    shell:
        """
        touch {output}
        FILES="{input.ibd}"
        TEMPFILE=ersa/temp_relatives.tsv
        rm -f $TEMPFILE
        rm -f {output}

        for input_file in $FILES; do
            ersa --avuncular-adj -ci -a {params.a} --dmax 14 -t {params.t} -l {params.l} \
                {params.r} -th {params.th} $input_file -o $TEMPFILE |& tee {log}

            if [[ "$input_file" == "${{FILES[0]}}" ]]; then
                cat $TEMPFILE >> {output}
            else
                sed 1d $TEMPFILE >> {output}
            fi
        done
        """


rule postprocess_ersa:
    input:
        ibd = rules.merge_rapid_segments.output['seg'],
        ersa = rules.ersa.output[0]
    output: "results/relatives.tsv"
    conda: "evaluation"
    log: "logs/merge/postprocess-ersa.log"
    script: "../scripts/postprocess_ersa.py"
