rule ibis:
    input:
        bed="preprocessed/data.bed",
        fam="preprocessed/data.fam",
        bim="preprocessed/data.bim"
    singularity:
        "docker://genxnetwork/ibis:stable"
    output:
        ibd     = "ibis/merged_ibis.seg"
    log:
        "logs/ibis/run_ibis.log"
    benchmark:
        "benchmarks/ibis/run_ibis.txt"
    threads: workflow.cores
    shell:
        """
        ibis {input.bed} {input.bim} {input.fam} -t {threads} -mt 500 -mL 7 -ibd2 -mL2 3 -hbd -f ibis/merged_ibis |& tee -a {log}
        """

rule transform_ibis_segments:
    input:
        ibd=rules.ibis.output.ibd
    output:
        germline = "ibd/merged_ibd.tsv"
    log:
        "logs/ibis/transform_ibis_segments.log"
    conda:
        "../envs/evaluation.yaml"
    script:
        "../scripts/transform_ibis_segments.py"

rule ersa:
    input:
        ibd=rules.transform_ibis_segments.output['germline']
    output:
        "ersa/relatives.tsv"
    singularity:
        "docker://genxnetwork/ersa:stable"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    params:
        ersa_l = config['zero_seg_count'],
        ersa_th = config['zero_seg_len'],
        alpha = config['alpha']
    shell:
        """
        ERSA_L={params.ersa_l} # the average number of IBD segments in population
        ERSA_TH={params.ersa_th} # the average length of IBD segment in population
        ERSA_T=7.0 # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -ci --alpha {params.alpha} --dmax 14 -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output}  |& tee {log}
        """

rule postprocess_ersa:
    input:
        ibd=rules.ibis.output['ibd'],
        ersa=rules.ersa.output[0]
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    log: "logs/merge/postprocess-ersa.log"
    script: "../scripts/postprocess_ersa.py"
