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
        bcftools filter {input} -r {wildcards.chrom} | bcftools norm --rm-dup none -O z -o vcf/imputed_chr{wildcards.chrom}.vcf.gz | tee {log}
        """


rule refined_ibd:
    input:
        vcf=rules.index_and_split.output,
        cmmap=beagle_map
    output:
        "ibd/chr{chrom}.ibd.gz"
    params:
        out="ibd/chr{chrom}"
    singularity:
        "docker://alexgenx/refined_ibd:latest"
    shell:
        """
            java -Xss5m -jar /app/refined-ibd.jar gt={input.vcf} map={input.cmmap} out={params.out}
        """

rule merge_matches:
    input:
         expand("ibd/chr{chrom}.ibd.gz", chrom=CHROMOSOMES)
    output:
         ibd="ibd/all.tsv"
    shell:
         "zcat ibd/*.ibd.gz > {output}"


rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        ibd=rules.merge_matches.output['ibd']
        #ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/ersa_king.py"
