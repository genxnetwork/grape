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
    singularity:
        "docker://lifebitai/king:latest"
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
    singularity:
        "docker://genxnetwork/ibis:stable"
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
        alpha = config['alpha'],
        ersa_t = config['ibis_seg_len']
    shell:
        """
        ERSA_L={params.ersa_l} # the average number of IBD segments in population
        ERSA_TH={params.ersa_th} # the average length of IBD segment in population
        ERSA_T={params.ersa_t} # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -ci --alpha {params.alpha} --dmax 14 -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output}  |& tee {log}
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
        ibd=rules.transform_ibis_segments.output['germline'],
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