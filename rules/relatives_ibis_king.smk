from os import path


rule run_king:
    input:
        bed="preprocessed/data.bed",
        bim="preprocessed/data.bim"
    output:
        king="king/data.seg",
        kinship="king/data.kin",
        kinship0="king/data.kin0",
        segments="king/data.segments.gz"
    params:
        out="king/data",
        kin="king/data"
    threads: workflow.cores
    conda:
        "../envs/king.yaml"
    log:
        "logs/king/run_king.log"
    benchmark:
        "benchmarks/king/run_king.txt"
    shell:
        """
        KING_DEGREE=3

        king -b {input.bed} --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out} |& tee {log}
        king -b {input.bed} --cpus {threads} --kinship --degree 4 --prefix {params.kin} |& tee -a {log}

        # we need at least an empty file for the downstream analysis
        if [ ! -f "{output.king}" ]; then
            touch {output.king}
        fi
        if [ ! -f "{output.kinship}" ]; then
            touch {output.kinship}
        fi
        if [ ! -f "{output.kinship0}" ]; then
            touch {output.kinship0}
        fi
        if [ ! -f "{output.segments}" ]; then
            touch {output.segments}
        fi
        """

rule ibis:
    input:
        bed="preprocessed/data.bed",
        fam="preprocessed/data.fam",
        bim="preprocessed/data_mapped.bim"
    conda:
        "../envs/ibis.yaml"
    output:
        ibd     = "ibis/merged_ibis.seg"
    log:
        "logs/ibis/run_ibis.log"
    benchmark:
        "benchmarks/ibis/run_ibis.txt"
    threads: workflow.cores
    params:
        mL=config['ibis_seg_len'],
        mT=config['ibis_min_snp']
    shell:
        """
        ibis {input.bed} {input.bim} {input.fam} -t {threads} -mt {params.mT} -mL {params.mL} -ibd2 -mL2 3 -hbd -f ibis/merged_ibis |& tee -a {log}
        """


WEIGHTED_IBD_SEGEMNTS_FOLDER = 'ibis-weighted'
SNAKEFILE_FOLDER = os.path.dirname(workflow.snakefile)


if config['weight_mask']:
    rule ibis_segments_weighing:
        input:
            ibd = rules.ibis.output.ibd,
            script = path.join(SNAKEFILE_FOLDER, '../weight/apply_weight_mask.py')
        conda:
            '../envs/weight-mask.yaml'
        output:
            ibd = path.join(WEIGHTED_IBD_SEGEMNTS_FOLDER, 'ibis_weighted.seg'),
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
        "../envs/evaluation.yaml"
    script:
        "../scripts/transform_ibis_segments.py"


def aggregate_input(wildcards):
    checkpoints.transform_ibis_segments.get()
    ids = glob_wildcards(f"ibd/{{id}}.tsv").id
    return expand(f"ibd/{{id}}.tsv", id=ids)

rule ersa:
    input:
        ibd=aggregate_input
    output:
        "ersa/relatives.tsv"
    conda:
        "../envs/ersa.yaml"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    params:
        ersa_l = config['zero_seg_count'],
        ersa_th = config['zero_seg_len'],
        alpha = config['alpha'],
        ersa_t = config['ibis_seg_len']
    shell:
        """
        ERSA_L={params.ersa_l} # the average number of IBD segments in population
        ERSA_TH={params.ersa_th} # the average length of IBD segment in population
        ERSA_T={params.ersa_t} # min length of segment to be considered in segment aggregation

        FILES="{input.ibd}"
        TEMPFILE=ersa/temp_relatives.tsv
        rm -f $TEMPFILE
        rm -f {output}

        for input_file in $FILES; do

            ersa --avuncular-adj -ci --alpha {params.alpha} --dmax 14 -t $ERSA_T -l $ERSA_L -th $ERSA_TH $input_file -o $TEMPFILE  |& tee {log}

            if [[ "$input_file" == "${{FILES[0]}}" ]]; then
                cat $TEMPFILE >> {output}
            else
                sed 1d $TEMPFILE >> {output}
            fi
        done
        """

rule split_map:
    input:
        bim = "preprocessed/data_mapped.bim"
    output: expand("cm/chr{chrom}.cm.map", chrom=CHROMOSOMES)
    params:
        cm_dir='cm'
    conda:
        "../envs/evaluation.yaml"
    script:
        "../scripts/split_map.py"

rule merge_king_ersa:
    input:
        king=rules.run_king.output['king'],
        king_segments=rules.run_king.output['segments'],
        ersa=rules.ersa.output[0],
        kinship=rules.run_king.output['kinship'],
        kinship0=rules.run_king.output['kinship0'],
        cm=expand("cm/chr{chrom}.cm.map", chrom=CHROMOSOMES)
    params:
        cm_dir='cm'
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    log: "logs/merge/merge-king-ersa.log"
    script: "../scripts/merge_king_ersa.py"
