
rule vcf_stats:
    input:
        vcf="vcf/merged_lifted.vcf.gz"
    output:
        stats="stats/lifted_vcf.txt",
        psc="stats/lifted_vcf.psc"
    params:
        samples="vcf/merged_lifted.vcf.samples"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
            bcftools query --list-samples {input.vcf} > {params.samples}
            bcftools stats -S {params.samples} {input.vcf} > {output.stats}
            # PSC means per-sample counts 
            cat {output.stats} | grep "^PSC" > {output.psc}  
        """

rule select_bad_samples:
    input:
        psc=rules.vcf_stats.output.psc
    output:
        bad_samples="vcf/lifted_vcf.badsamples",
        report="results/bad_samples_report.tsv"
    log: "logs/vcf/select_bad_samples.log"
    params:
        samples="vcf/merged_lifted.vcf.samples"
    conda:
        "../envs/evaluation.yaml"
    script:
        "../scripts/select_bad_samples.py"

rule plink_filter:
    input:
        vcf="vcf/merged_lifted.vcf.gz",
        bad_samples=rules.select_bad_samples.output.bad_samples
    output: temp(expand("plink/merged_filter.{ext}", ext=PLINK_FORMATS))
    conda:
        "../envs/plink.yaml"
    params:
        input   = "merged",
        out     = "merged_filter"
    log:
        "logs/plink/plink_filter.log"
    benchmark:
        "benchmarks/plink/plink_filter.txt"
    shell:
        """
        plink --vcf {input.vcf} --freqx --out plink/{params.out}
        plink --vcf {input.vcf} --remove {input.bad_samples} --geno 0.5 --maf 0.02 --hwe 0 --make-bed --keep-allele-order --out plink/{params.out} |& tee {log}
        """

rule pre_imputation_check:
    input:
        "plink/merged_filter.bim"
    params:
        SITE_1000GENOME
    output:
        "plink/merged_filter.bim.chr",
        "plink/merged_filter.bim.pos",
        "plink/merged_filter.bim.force_allele",
        "plink/merged_filter.bim.flip"
    log:
        "logs/plink/pre_imputation_check.log"
    benchmark:
        "benchmarks/plink/pre_imputation_check.txt"
    script:
        "../scripts/pre_imputation_check.py"

rule plink_clean_up:
    input:
        "plink/merged_filter.bim.chr",
        "plink/merged_filter.bim.pos",
        "plink/merged_filter.bim.force_allele",
        "plink/merged_filter.bim.flip",
        plink=expand("plink/merged_filter.{ext}", ext=PLINK_FORMATS)
    output:
        expand("{i}.{ext}", i="plink/merged_mapped", ext=PLINK_FORMATS)
    params:
        input = "plink/merged_filter",
        out = "plink/merged_mapped"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/plink_clean_up.log"
    benchmark:
        "benchmarks/plink/plink_clean_up.txt"
    shell:
        """
        # remove dublicates
        cut -f 2 {params.input}.bim | sort | uniq -d > plink/snp.dups
        plink --bfile {params.input}          --exclude       plink/snp.dups                  --make-bed --out plink/merged_filter_dub   |& tee -a {log}
        plink --bfile plink/merged_filter_dub --extract       plink/merged_filter.bim.chr     --make-bed --out plink/merged_extracted    |& tee -a {log}
        plink --bfile plink/merged_extracted  --flip          plink/merged_filter.bim.flip    --make-bed --out plink/merged_flipped      |& tee -a {log}
        plink --bfile plink/merged_flipped    --update-chr    plink/merged_filter.bim.chr     --make-bed --out plink/merged_chroped      |& tee -a {log}
        plink --bfile plink/merged_chroped    --update-map    plink/merged_filter.bim.pos     --make-bed --out {params.out}              |& tee -a {log}
        rm plink/merged_filter_dub.*
        rm plink/merged_extracted.*
        rm plink/merged_flipped.*
        rm plink/merged_chroped.*
        """

rule prepare_vcf:
    input: "plink/merged_mapped.bim"
    output:
        vcf="vcf/merged_mapped_sorted.vcf.gz",
        temp_vcf=temp("vcf/merged_mapped_regions.vcf.gz"),
        temp_vcf_csi=temp("vcf/merged_mapped_regions.vcf.gz.csi"),
        alleled=temp(expand('plink/merged_mapped_alleled.{ext}', ext=PLINK_FORMATS))
    params:
        input   = "plink/merged_mapped",
        vcf     = 'vcf/merged_mapped_sorted.vcf.gz'
    conda:
         "../envs/bcf_plink.yaml"
    log:
        plink="logs/plink/prepare_vcf.log",
        vcf="logs/vcf/prepare_vcf.log"
    benchmark:
        "benchmarks/plink/prepare_vcf.txt"
    shell:
        """
        plink --bfile {params.input} --a1-allele plink/merged_filter.bim.force_allele --make-bed --out plink/merged_mapped_alleled |& tee -a {log.plink}
        plink --bfile plink/merged_mapped_alleled --keep-allele-order --output-chr M --export vcf bgz --out vcf/merged_mapped_clean |& tee -a {log.vcf}
        bcftools sort vcf/merged_mapped_clean.vcf.gz -O z -o {output.temp_vcf} |& tee -a {log.vcf}
        bcftools index -f {output.temp_vcf} |& tee -a {log.vcf}
        # need to check output for the potential issues
        bcftools view {output.temp_vcf} --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -O z -o {output.vcf}
        bcftools norm --check-ref e -f {GRCH37_FASTA} vcf/merged_mapped_sorted.vcf.gz -O u -o /dev/null |& tee -a {log.vcf}
        bcftools index -f {output.vcf} | tee -a {log.vcf}
        """