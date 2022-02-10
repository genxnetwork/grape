CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', ]
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'nosex']


rule impute:
    input:
        rules.phase.output,
        refHaps=REF_HAPS
    output: temp("imputed/chr{chrom}.imputed.dose.vcf.gz")
    threads: 1
    log:
        "logs/impute/minimac4-{chrom}.log"
    benchmark:
        "benchmarks/impute/minimac4-{chrom}.txt"
    shell:
        """
            minimac4 \
            --refHaps {input.refHaps} \
            --haps phase/chr{wildcards.chrom}.phased.vcf.gz \
            --format GT,GP \
            --prefix imputed/chr{wildcards.chrom}.imputed \
            --minRatio 0.01 \
            --cpus {threads} |& tee {log}
        """

rule imputation_filter:
    input: rules.impute.output
    output: temp("imputed/chr{chrom}.imputed.dose.pass.vcf.gz")
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/impute/imputation_filter-{chrom}.log"
    benchmark:
        "benchmarks/impute/imputation_filter-{chrom}.txt"
    shell:
        """
        FILTER="'strlen(REF)=1 & strlen(ALT)=1'"

        bcftools view -i 'strlen(REF)=1 & strlen(ALT)=1' imputed/chr{wildcards.chrom}.imputed.dose.vcf.gz -v snps -m 2 -M 2 -O z -o imputed/chr{wildcards.chrom}.imputed.dose.pass.vcf.gz |& tee {log}
        """


rule merge_imputation_filter:
    input:
        expand("imputed/chr{i}.imputed.dose.pass.vcf.gz", i=CHROMOSOMES)
        #expand("phase/chr{chrom}.phased.vcf.gz", chrom=CHROMOSOMES)
        # TODO: wildcard violation
        # rules.imputation_filter.output
    output:
        "preprocessed/data.vcf.gz"
    params:
        list="vcf/imputed.merge.list",
        mode=config["mode"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf/merge_imputation_filter.log"
    benchmark:
        "benchmarks/vcf/merge_imputation_filter.txt"
    shell:
        """
        # for now just skip empty files
        true > {params.list} && \
        for i in {input}; do
            if [ -s $i ]
            then
                echo $i >> {params.list}
            else
                continue
            fi
        done

        bcftools concat -f {params.list} -O z -o {output} |& tee -a {log}
        bcftools index -f {output} |& tee -a {log}

        # check if there is a background data and merge it
        if [ -f "background/merged_imputed.vcf.gz" && {params.mode} = "client" ]; then
            mv {output} {output}.client
            bcftools merge --force-samples background/merged_imputed.vcf.gz {output}.client -O z -o {output} |& tee -a {log}
            bcftools index -f {output} |& tee -a {log}
        fi
        """

rule convert_imputed_to_plink:
    input: rules.merge_imputation_filter.output
    output: expand("plink/{i}.{ext}", i="merged_imputed", ext=PLINK_FORMATS)
    params:
        out = "plink/merged_imputed"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/convert_imputed_to_plink.log"
    benchmark:
        "benchmarks/plink/convert_imputed_to_plink.txt"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out} |& tee {log}
        """

# no need it bc it was done earlier in merge_imputation_filter
rule merge_convert_imputed_to_plink:
    input: rules.merge_imputation_filter.output
    output: expand("plink/{i}.{ext}", i="merged_imputed", ext=PLINK_FORMATS)
    params:
        background  = "background/merged_imputed",
        out         = "plink/client/merged_imputed"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/convert_imputed_to_plink.log"
    benchmark:
        "benchmarks/plink/convert_imputed_to_plink.txt"
    shell:
        """
        # please mind a merge step in merge_imputation_filter for germline
        plink --vcf {input} --make-bed --out {params.out} | tee {log}
        plink --bfile {params.background} --bmerge {params.out} --make-bed --out plink/merged_imputed |& tee {log}
        """
