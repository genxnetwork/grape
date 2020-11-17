CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']


rule phase:
    input:
        vcf="vcf/merged_mapped_sorted.vcf.gz",
        #idx="vcf/merged_mapped_sorted.vcf.gz.csi"
        vcfRef=vcfRef
    output: "phase/chr{chrom}.phased.vcf.gz"
    threads: 1
    singularity:
        "docker://biocontainers/bio-eagle:v2.4.1-1-deb_cv1"
    log:
        "logs/phase/eagle-{chrom}.log"
    benchmark:
        "benchmarks/phase/eagle-{chrom}.txt"
    shell:
        """
        /usr/bin/bio-eagle --vcfRef {input.vcfRef} \
        --numThreads {threads} \
        --vcfTarget {input.vcf}  \
        --geneticMapFile {GENETIC_MAP} \
        --chrom {wildcards.chrom} \
        --vcfOutFormat z \
        --outPrefix phase/chr{wildcards.chrom}.phased | tee {log}
        """

rule impute:
    input:
        rules.phase.output,
        refHaps=refHaps
    output: "imputed/chr{chrom}.imputed.dose.vcf.gz"
    threads: 1
    singularity:
        "docker://biocontainers/minimac4:v1.0.0-2-deb_cv1"
    log:
        "logs/impute/minimac4-{chrom}.log"
    benchmark:
        "benchmarks/impute/minimac4-{chrom}.txt"
    shell:
        """
        /usr/bin/minimac4 \
        --refHaps {input.refHaps} \
        --haps phase/chr{wildcards.chrom}.phased.vcf.gz \
        --format GT,GP \
        --prefix imputed/chr{wildcards.chrom}.imputed \
        --cpus {threads} | tee {log}
        """


rule imputation_filter:
    input: rules.impute.output
    output: "imputed/chr{chrom}.imputed.dose.pass.vcf.gz"
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
        FILTER="'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1'"

        bcftools view -i 'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1' imputed/chr{wildcards.chrom}.imputed.dose.vcf.gz -v snps -m 2 -M 2 -O z -o imputed/chr{wildcards.chrom}.imputed.dose.pass.vcf.gz | tee {log}
        """


rule merge_imputation_filter:
    input:
        expand("imputed/chr{i}.imputed.dose.pass.vcf.gz", i=CHROMOSOMES)
        #expand("phase/chr{chrom}.phased.vcf.gz", chrom=CHROMOSOMES)
        # TODO: wildcard violation
        # rules.imputation_filter.output
    output:
        "vcf/merged_imputed.vcf.gz"
    params:
        list="vcf/imputed.merge.list"
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

        # check if there is a background data
        # we need to merge separately from merge_convert_imputed_to_plink
        # because it used in index_and_split down to germline
        if [ -f "background/merged_imputed.vcf.gz" ]; then
            echo "background/merged_imputed.vcf.gz" >> {params.list}
        done

        bcftools concat -f {params.list} -O z -o {output} | tee -a {log}
        bcftools index -f {output} | tee -a {log}
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
        plink --vcf {input} --make-bed --out {params.out} | tee {log}
        """

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
        plink --bfile {params.background} --bmerge {params.out} --make-bed --out plink/merged_imputed | tee {log}
        """
