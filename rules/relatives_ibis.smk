rule convert_mapped_to_plink:
    input:
        vcf="preprocessed/data.vcf.gz"
    output:
        plink=expand("plink/{i}.{ext}", i="merged_ibis", ext=PLINK_FORMATS)
    params:
        out = "plink/merged_ibis"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/convert_mapped_to_plink.log"
    benchmark:
        "benchmarks/plink/convert_mapped_to_plink.txt"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out} |& tee {log}
        """

rule ibis_mapping:
    input:
        rules.convert_mapped_to_plink.output
    params:
        input = "plink/merged_ibis"
    singularity:
        "docker://genxnetwork/ibis:stable"
    output:
        "plink/merged_ibis_mapped.bim"
    log:
        "logs/ibis/run_ibis_mapping.log"
    benchmark:
        "benchmarks/ibis/run_ibis_mapping.txt"
    shell:
        """
        (add-map-plink.pl -cm {params.input}.bim {GENETIC_MAP_GRCH37}> plink/merged_ibis_mapped.bim) |& tee -a {log}
        """

rule ibis:
    input:
        rules.ibis_mapping.output
    params:
        input = "plink/merged_ibis"
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
        ibis {params.input}.bed {input} {params.input}.fam -t {threads} -mt 560 -mL 7 -ibd2 -mL2 3 -hbd -f ibis/merged_ibis |& tee -a {log}
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
    shell:
        """
        ERSA_L=0.5 # the average number of IBD segments in population
        ERSA_TH=6.00 # the average length of IBD segment in population
        ERSA_T=7.0 # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -ci --alpha 0.01 --dmax 14 -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output}  |& tee {log}
        """

rule postprocess_ersa:
    input:
        ibd=rules.ibis.output['ibd'],
        ersa=rules.ersa.output[0]
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    log: "logs/merge/postprocess-ersa.log"
    script: "../scripts/postprocess_ersa.py"
