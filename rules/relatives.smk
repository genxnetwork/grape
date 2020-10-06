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
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    log:
        "logs/vcf/index_and_split-{chrom}.log"
    benchmark:
        "benchmarks/vcf/index_and_split-{chrom}.txt"
    shell:
        """
        bcftools filter {input} -r {wildcards.chrom} -O z -o vcf/imputed_chr{wildcards.chrom}.vcf.gz | tee {log}
        """

rule convert_to_hap:
    input: rules.index_and_split.output
    output: "hap/imputed_chr{chrom}.hap"
    conda:
        "../envs/bcftools.yaml"
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    log:
        "logs/vcf/convert_to_hap-{chrom}.log"
    benchmark:
        "benchmarks/vcf/convert_to_hap-{chrom}.txt"
    shell:
        """
        bcftools convert vcf/imputed_chr{wildcards.chrom}.vcf.gz --hapsample hap/imputed_chr{wildcards.chrom} | tee {log}
        gunzip -f hap/imputed_chr{wildcards.chrom}.hap.gz | tee -a {log}
        """

rule convert_to_ped:
    input: rules.convert_to_hap.output
    output: "ped/imputed_chr{chrom}.ped"
    singularity:
        "docker://alexgenx/germline:stable"
    log:
        "logs/ped/convert_to_ped-{chrom}.log"
    benchmark:
        "benchmarks/ped/convert_to_ped-{chrom}.txt"
    shell:
        """
        impute_to_ped hap/imputed_chr{wildcards.chrom}.hap hap/imputed_chr{wildcards.chrom}.sample ped/imputed_chr{wildcards.chrom} | tee {log}
        """

# TODO: skip this step for now
rule split_map:
    input: rules.convert_to_ped.output
    output: expand("chr{i}.map", i=CHROMOSOMES)
    shell:
        """
        echo parallel -k -j 100% \
        grep ^{{}} {{}}.map '>' chr{{}}.map \
        ::: `seq 1 22`

        for i in {output}; do
            touch $i;
        done
        """

rule interpolate:
    input: rules.convert_to_ped.output
    output: "cm/chr{chrom}.cm.ped"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/cm/interpolate-{chrom}.log"
    benchmark:
        "benchmarks/cm/interpolate-{chrom}.txt"
    shell:
        """
        plink --file ped/imputed_chr{wildcards.chrom} --cm-map /media/ref/genetic_map_b37/genetic_map_chr{wildcards.chrom}_combined_b37.txt {wildcards.chrom} --recode --out cm/chr{wildcards.chrom}.cm | tee {log}
        """

rule germline:
    input: rules.interpolate.output
    output: "germline/chr{chrom}.germline.match"
    singularity:
        "docker://alexgenx/germline:stable"
    log:
        "logs/germline/germline-{chrom}.log"
    benchmark:
        "benchmarks/germline/germline-{chrom}.txt"
    shell:
        """
        germline -input ped/imputed_chr{wildcards.chrom}.ped cm/chr{wildcards.chrom}.cm.map -homoz -min_m 2.5 -err_hom 2 -err_het 1 -output germline/chr{wildcards.chrom}.germline | tee {log}
        # TODO: germline returns some length in BP instead of cM - clean up is needed
        grep -v MB germline/chr{wildcards.chrom}.germline.match > germline/chr{wildcards.chrom}.germline.match.clean && mv germline/chr{wildcards.chrom}.germline.match.clean germline/chr{wildcards.chrom}.germline.match
        """

rule ersa_params:
    input:
        king=rules.run_king.output,
        # TODO: wildcard violation
        # germline=rules.germline.output
        germline=expand("germline/chr{i}.germline.match", i=CHROMOSOMES)
    output: "ersa/params.txt"
    conda: "../envs/ersa_params.yaml"
    script: "../scripts/estimate_ersa_params.py"

rule merge_matches:
    input:
         expand("germline/chr{chrom}.germline.match", chrom=CHROMOSOMES)
    output:
          "germline/all.tsv"
    shell:
         "cat germline/*.match > {output}"

rule merge_ibd_segments:
    input:
        germline=rules.merge_matches.output[0]
    params:
        cm_dir='cm',
        merge_gap='0.6'
    output:
        ibd='germline/merged_ibd.tsv'
    conda: "../envs/evaluation.yaml"
    script:
        '../scripts/merge_ibd.py'

rule ersa:
    input:
        germline=rules.merge_ibd_segments.output['ibd'],
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
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.germline} -o {output} | tee {log}
        """


rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        germline=rules.merge_ibd_segments.output['ibd'],
        ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/ersa_king.py"
