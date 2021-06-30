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

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out} |& tee {log}
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
        bcftools filter {input} -r {wildcards.chrom} | bcftools norm --rm-dup none -O z -o vcf/imputed_chr{wildcards.chrom}.vcf.gz |& tee {log}
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

rule erase_dosages:
    input:
        vcf=rules.index_and_split.output[0]
    output:
        vcf='vcf/erased_chr{chrom}.vcf.gz'
    params:
        vcf='vcf/erased_chr{chrom}'
    conda:
        '../envs/bcftools.yaml'
    shell:
        "bcftools annotate -x 'fmt' {input.vcf} -O z -o {output.vcf}"

rule rapid:
    input:
        vcf=rules.erase_dosages.output['vcf'],
        g=rules.interpolate.output
    singularity:
        "docker://genxnetwork/rapid:stable"
    params:
        min_cm_len=1.0,
        window_size=50,
        output_folder='rapid/chr{chrom}'
    output:
        "rapid/chr{chrom}/results.max.gz"
    shell:
        "rapid -i {input.vcf} -g {input.g} -d {params.min_cm_len} -o {params.output_folder} -w {params.window_size} -r 10 -s 6"


rule merge_ibd_segments:
    input:
        expand("rapid/chr{chrom}/results.max.gz", chrom=CHROMOSOMES)
    output:
        ibd='ibd/merged_ibd.tsv'
    conda: "../envs/evaluation.yaml"
    script:
        '../scripts/merge_rapid_ibd.py'

rule rapid_to_ersa:
    input:
        ibd=rules.merge_ibd_segments.output['ibd'],
        true_ibd=rules.simulate.output.seg
    params:
        cm_dir='cm',
        merge_gap='0.2',
        use_true_ibd=False
    output:
        ibd='ibd/ersa_ibd.tsv'
    conda: "../envs/evaluation.yaml"
    script:
        '../scripts/rapid_to_ersa.py'

rule ersa:
    input:
        ibd=rules.rapid_to_ersa.output['ibd']
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
        ERSA_L=5.0 # the average number of IBD segments in population
        ERSA_TH=1.5 # the average length of IBD segment
        ERSA_T=1.0 # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output} |& tee {log}
        """

rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        germline=rules.rapid_to_ersa.output['ibd'],
        ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/ersa_king.py"
