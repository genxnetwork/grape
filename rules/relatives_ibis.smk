rule convert_mapped_to_plink:
    input: rules.prepare_vcf.output
    output: expand("plink/{i}.{ext}", i="merged_ibis", ext=PLINK_FORMATS)
    params:
        out = "plink/merged_ibis"
    conda:
        "../envs/bcf_plink.yaml"
    log:
        "logs/plink/convert_mapped_to_plink.log"
    benchmark:
        "benchmarks/plink/convert_mapped_to_plink.txt"
    shell:
        """
        # leave only chr1..22 because we need to map it later
        bcftools view {input} --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -O z -o vcf/merged_mapped_sorted_22.vcf.gz
        plink --vcf vcf/merged_mapped_sorted_22.vcf.gz --make-bed --out {params.out} |& tee {log}
        """

rule run_king:
    input: rules.convert_mapped_to_plink.output
    output: "king/merged_ibis_king.seg"
    params:
        input = "plink/merged_ibis",
        out = "king/merged_ibis_king"
    threads: workflow.cores
    singularity:
        "docker://lifebitai/king:latest"
    log:
        "logs/king/run_king.log"
    benchmark:
        "benchmarks/king/run_king.txt"
    shell:
        """
        # TODO: add cores
        KING_DEGREE=4

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out} |& tee {log}
        """

rule ibis:
    input:
        rules.convert_mapped_to_plink.output
    params:
        input = "plink/merged_ibis"
    singularity:
        "docker://alexgenx/ibis:stable"
    output:
        ibd     = "ibis/merged_ibis.seg",
        germline= "ibis/merged_ibis.germline"
    shell:
        """
        add-map-plink.pl {params.input}.bim {genetic_map_GRCh37} > plink/merged_ibis_mapped.bim
        # use default params
        ibis {params.input}.bed plink/merged_ibis_mapped.bim {params.input}.fam -f ibis/merged_ibis

        cat ibis/merged_ibis.seg | awk -v OFS='\t' '{{sub(":", "_", $1); sub(":", "_", $2); print $1, $1, $2, $2, $3, $4, $5, 0, 0, $10, $9, "cM", 0, 0, 0}};' > ibis/merged_ibis.germline
        """

rule ersa:
    input:
        ibd=rules.ibis.output['germline']
    output:
        "ersa/relatives.tsv"
    singularity:
        "docker://alexgenx/ersa:stable"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    shell:
        """
        ERSA_L=5.0 # the average number of IBD segments in population
        ERSA_TH=1.5 # the average length of IBD segment
        ERSA_T=1.0 # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output} |& tee {log}
        """

rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        germline=rules.ibis.output['germline'],
        ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/ersa_king.py"
