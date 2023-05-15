import os
import sys


rule split_for_ibis:
    input:
        bed = "preprocessed/data.bed",
        fam = "preprocessed/data.fam",
        bim = "preprocessed/data_mapped.bim"
    conda:
        "plink"
    output:
        bed = "ibis/{chr}.bed",
        fam = "ibis/{chr}.fam",
        bim = "ibis/{chr}.bim"
    threads: 1
    params:
        bfile = "preprocessed/data",
        chr = lambda wildcards: wildcards.chr
    shell:
        """
            plink --bfile {params.bfile} --chr {params.chr} --make-bed --out ibis/{params.chr} --threads {threads}
        """


rule ibis:
    input:
        bed = "ibis/{chr}.bed",
        fam = "ibis/{chr}.fam",
        bim = "ibis/{chr}.bim"
    conda:
        "ibis"
    output:
        ibd = "ibis/merged_ibis_{chr}.seg"
    log:
        "logs/ibis/run_ibis_{chr}.log"
    benchmark:
        "benchmarks/ibis/run_ibis_{chr}.txt"
    threads: workflow.cores
    params:
        mL = config['ibis_seg_len'],
        mT = config['ibis_min_snp'],
        chr = lambda wildcards: wildcards.chr
    shell:
        """
            ibis {input.bed} {input.bim} {input.fam} -t {threads} -mt {params.mT} -mL {params.mL} -ibd2 -mL2 3 -hbd -f ibis/merged_ibis_{params.chr} |& tee -a {log}
        """


rule merge_ibis_segments:
    input: 
        expand("ibis/merged_ibis_{chr}.seg", chr=CHROMOSOMES)
    output:
        "ibis/merged_ibis.seg"
    shell:
        """
            cat {input} >> {output}
        """


WEIGHTED_IBD_SEGMENTS_FOLDER = 'ibis-weighted'
SNAKEFILE_FOLDER = os.path.dirname(workflow.snakefile)


if config.get('weight_mask'):
    rule ibis_segments_weighing:
        input:
            ibd = rules.ibis.output.ibd,
            script = os.path.join(SNAKEFILE_FOLDER, '../weight/apply_weight_mask.py')
        conda:
            'weight-mask'
        output:
            ibd = os.path.join(WEIGHTED_IBD_SEGMENTS_FOLDER, 'ibis_weighted.seg'),
        params:
            mask = config['weight_mask']
        shell:
            """
            python {input.script} \
                --input-ibd-segments-file {input.ibd} \
                --mask-file {params.mask} \
                --output-ibd-segments-file {output.ibd}
            """
    ibd_segments_file = rules.ibis_segments_weighing.output.ibd
else:
    ibd_segments_file = rules.ibis.output.ibd


checkpoint transform_ibis_segments:
    input:
        ibd = ibd_segments_file,
        fam = "preprocessed/data.fam"
    output:
        bucket_dir = directory("ibd")
    log:
        "logs/ibis/transform_ibis_segments.log"
    conda:
        "evaluation"
    script:
        "../scripts/transform_ibis_segments.py"


def aggregate_input(wildcards):
    checkpoints.transform_ibis_segments.get()
    ids = glob_wildcards(f"ibd/{{id}}.tsv").id
    return expand(f"ibd/{{id}}.tsv", id=ids)


rule ersa:
    input:
        ibd = "ibd/{id}.tsv"
    output:
        relatives="ersa/relatives_{id}.tsv"
    conda:
        "ersa"
    log:
        "logs/ersa/ersa_{id}.log"
    benchmark:
        "benchmarks/ersa/ersa_{id}.txt"
    params:
        l = config['zero_seg_count'],
        th = config['zero_seg_len'],
        a = config['alpha'],
        t = config['ibis_seg_len'],
        r = '--nomask ' + '-r ' + str(config['ersa_r']) if config.get('weight_mask') else ''
    shell:
        """
        ersa --avuncular-adj -ci -a {params.a} --dmax 14 -t {params.t} -l {params.l} \
                {params.r} -th {params.th} {input.ibd} -o {output.relatives} |& tee {log}
        """


rule merge_ersa:
    input:
        expand("ersa/relatives_{id}.tsv", id=glob_wildcards("ibd/{id}.tsv").id)
    output:
        "ersa/relatives.tsv"
    shell:
        """
        FILES="{input}"
        for input_file in $FILES; do
            if [[ "$input_file" == "${{FILES[0]}}" ]]; then
                cat $input_file >> {output}
            else
                sed 1d $input_file >> {output}
            fi
        done
        """


rule postprocess_ersa:
    input:
        ibd=rules.merge_ibis_segments.output[0],
        ersa=rules.merge_ersa.output[0]
    params:
        ibis = True
    output: "results/relatives.tsv"
    conda: "evaluation"
    log: "logs/merge/postprocess-ersa.log"
    script: "../scripts/postprocess_ersa.py"
