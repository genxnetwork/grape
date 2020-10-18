rule run_king:
    input: rules.convert_imputed_to_plink.output
    output: "king/merged_imputed_king.seg"
    params:
        input = "plink/merged_imputed",
        out = "king/merged_imputed_king"
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

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out} | tee {log}
        """

rule index_and_split:
    input: rules.merge_imputation_filter.output
    output: "vcf/imputed_chr{chrom}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf/index_and_split-{chrom}.log"
    benchmark:
        "benchmarks/vcf/index_and_split-{chrom}.txt"
    shell:
        """
        bcftools filter {input} -r {wildcards.chrom} -O z -o vcf/imputed_chr{wildcards.chrom}.vcf.gz | tee {log}
        """

rule interpolate:
    input:
        vcf=rules.index_and_split.output[0],
        cmmap=cmmap
    output: "cm/chr{chrom}.cm.g"
    conda:
        "../envs/interpolation.yaml"
    log:
        "logs/cm/interpolate-{chrom}.log"
    script:
        "../scripts/interpolate.py"

rule rapid:
    input:
        vcf=rules.index_and_split.output,
        g=rules.interpolate.output
    singularity:
        "docker://alexgenx/rapid:latest"
    params:
        min_cm_len=2.5,
        window_size=3,
        output_folder='rapid/chr{chrom}'
    output:
        "rapid/chr{chrom}/results.max.gz"
    shell:
        "rapid -i {input.vcf} -g {input.g} -d {params.min_cm_len} -o {params.output_folder} -w {params.window_size} -r 10 -s 2"


rule merge_rapid_segments:
    input:
        expand("rapid/chr{chrom}/results.max.gz", chrom=CHROMOSOMES)
    output:
        ibd='rapid/merged_ibd.tsv'
    conda: "../envs/evaluation.yaml"
    script:
        '../scripts/merge_rapid_ibd.py'

rule ersa_params:
    input:
        king=rules.run_king.output
    output: "ersa/params.txt"
    conda: "../envs/ersa_params.yaml"
    script: "../scripts/estimate_ersa_params.py"

rule ersa:
    input:
        rapid=rules.merge_rapid_segments.output['ibd'],
        estimated=rules.ersa_params.output
    output: "ersa/relatives.tsv"
    singularity:
        "docker://alexgenx/ersa:stable"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    shell:
        """
        #ERSA_L=13.7
        #ERSA_TH=3.197
        #ERSA_T=2.5
        source {input.estimated}
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.rapid} -o {output} | tee {log}
        """

rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        germline=rules.merge_rapid_segments.output['ibd'],
        ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/ersa_king.py"
