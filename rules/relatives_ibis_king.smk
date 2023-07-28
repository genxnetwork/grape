import os


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


WEIGHTED_IBD_SEGMENTS_FOLDER = 'ibis-weighted'
SNAKEFILE_FOLDER = os.path.dirname(workflow.snakefile)


if config.get('weight_mask'):
    rule ibis_segments_weighing:
        input:
            ibd = rules.ibis.output.ibd,
            script = os.path.join(SNAKEFILE_FOLDER, '../weight/apply_weight_mask.py')
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
    script:
        "../scripts/transform_ibis_segments.py"


def tis_output(wildcards):
    checkpoints.transform_ibis_segments.get(id=wildcards.id)
    return "ibd/{id}.tsv.gz"


rule ersa:
    input:
        ibd = tis_output
    output:
        ibd = temp('temp_ibd/{id}.tsv'),
        relatives="ersa/relatives_{id}.tsv"
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
        zcat {input.ibd} > {output.ibd}
        ersa --avuncular-adj -ci -a {params.a} --dmax 14 -t {params.t} -l {params.l} \
                {params.r} -th {params.th} {output.ibd} -o {output.relatives} |& tee {log}
        """


def aggregate_input(wildcards):
    checkpoints.transform_ibis_segments.get()
    ids = glob_wildcards(f"ibd/{{id}}.tsv.gz").id
    return expand(f"ersa/relatives_{{id}}.tsv", id=ids)


rule merge_ersa:
    input:
        aggregate_input
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


rule split_map:
    input:
        bim = "preprocessed/data_mapped.bim"
    output: expand("cm/chr{chrom}.cm.map", chrom=CHROMOSOMES)
    params:
        cm_dir='cm'
    script:
        "../scripts/split_map.py"

rule merge_king_ersa:
    input:
        king=rules.run_king.output['king'],
        king_segments=rules.run_king.output['segments'],
        ersa=rules.merge_ersa.output[0],
        kinship=rules.run_king.output['kinship'],
        kinship0=rules.run_king.output['kinship0'],
        cm=expand("cm/chr{chrom}.cm.map", chrom=CHROMOSOMES)
    params:
        cm_dir='cm'
    output: "results/relatives.tsv"
    log: "logs/merge/merge-king-ersa.log"
    script: "../scripts/merge_king_ersa.py"
