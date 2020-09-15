CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']


rule phase:
    input:
        vcf="vcf/merged_mapped_sorted.vcf.gz"
        #idx="vcf/merged_mapped_sorted.vcf.gz.csi"
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
        GENETIC_MAP=/media/ref/tables/genetic_map_hg19_withX.txt.gz

        /usr/bin/bio-eagle --vcfRef  /media/ref/1000genome/bcf/1000genome_chr{wildcards.chrom}.bcf \
        --numThreads {threads} \
        --vcfTarget {input.vcf}  \
        --geneticMapFile $GENETIC_MAP \
        --chrom {wildcards.chrom} \
        --vcfOutFormat z \
        --outPrefix phase/chr{wildcards.chrom}.phased | tee {log}
        """

rule impute:
    input: rules.phase.output
    output: "imputed/chr{chrom}.imputed.dose.vcf.gz"
    threads: workflow.cores
    singularity:
        "docker://biocontainers/minimac4:v1.0.0-2-deb_cv1"
    log:
        "logs/impute/minimac4-{chrom}.log"
    benchmark:
        "benchmarks/impute/minimac4-{chrom}.txt"
    shell:
        """
        /usr/bin/minimac4 \
        --refHaps /media/ref/Minimac/{wildcards.chrom}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
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
