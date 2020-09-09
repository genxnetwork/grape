
rule plink_filter:
    input: rules.liftover.output
    output: expand("plink/merged_filter.{ext}", ext=PLINK_FORMATS)
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
        # TODO: the old parameter --maf 1e-50 is too low, back to the "default" 0.05
        plink --vcf vcf/{input} --geno 0.5 --maf 0.05 --hwe 0 --make-bed --keep-allele-order --out plink/{params.out} | tee {log}
        """

rule pre_imputation_check:
    input: "plink/merged_filter.bim"
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
        "scripts/pre_imputation_check.py"

rule plink_clean_up:
    input:
        "plink/merged_filter.bim.chr",
        "plink/merged_filter.bim.pos",
        "plink/merged_filter.bim.force_allele",
        "plink/merged_filter.bim.flip"
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
        plink --bfile {params.input}         --extract       plink/merged_filter.bim.chr     --make-bed --out plink/merged_extracted    | tee -a {log}
        plink --bfile plink/merged_extracted --flip          plink/merged_filter.bim.flip    --make-bed --out plink/merged_flipped      | tee -a {log}
        plink --bfile plink/merged_flipped   --update-chr    plink/merged_filter.bim.chr     --make-bed --out plink/merged_chroped      | tee -a {log}
        plink --bfile plink/merged_chroped   --update-map    plink/merged_filter.bim.pos     --make-bed --out {params.out}              | tee -a {log}
        """

rule prepare_vcf:
    input: "plink/merged_mapped.bim"
    output:
        "vcf/merged_mapped_sorted.vcf.gz"
    params:
        input = "plink/merged_mapped"
    conda:
         "../envs/bcf_plink.yaml"
    log:
        plink="logs/plink/prepare_vcf.log",
        vcf="logs/vcf/prepare_vcf.log"
    benchmark:
        "benchmarks/plink/prepare_vcf.txt"
    shell:
        """
        GRCh37_fasta=/media/human_g1k_v37.fasta

        plink --bfile {params.input} --a1-allele plink/merged_filter.bim.force_allele --make-bed --out plink/merged_mapped_alleled | tee -a {log.plink}
        plink --bfile plink/merged_mapped_alleled --keep-allele-order --output-chr M --export vcf bgz --out vcf/merged_mapped_clean | tee -a {log.vcf}
        bcftools sort vcf/merged_mapped_clean.vcf.gz -O z -o vcf/merged_mapped_sorted.vcf.gz | tee -a {log.vcf}
        # need to check output for the potential issues
        # bcftools norm --check-ref e -f $GRCh37_fasta merged_mapped_sorted.vcf.gz -O u -o /dev/null
        bcftools index -f vcf/merged_mapped_sorted.vcf.gz | tee -a {log.vcf}
        """