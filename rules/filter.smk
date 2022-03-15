rule vcf_stats:
    input:
        vcf="vcf/{segment}.vcf.gz"
    output:
        stats="stats/{segment}_lifted_vcf.txt",
        psc="stats/{segment}_lifted_vcf.psc"
    params:
        samples="vcf/{segment}_merged_lifted.vcf.samples"
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
        bad_samples="vcf/{segment}_lifted_vcf.badsamples",
        report="results/{segment}_bad_samples_report.tsv"
    log: "logs/vcf/{segment}_select_bad_samples.log"
    params:
        samples="vcf/{segment}_merged_lifted.vcf.samples"
    conda:
        "../envs/evaluation.yaml"
    script:
        "../scripts/select_bad_samples.py"

rule plink_filter:
    input:
        vcf="vcf/{segment}_merged_lifted_id.vcf.gz",
        bad_samples=rules.select_bad_samples.output.bad_samples
    output:
        bed = temp("plink/{segment}_merged_filter.bed"),
        bim = temp("plink/{segment}_merged_filter.bim"),
        fam = temp("plink/{segment}_merged_filter.fam")
    conda:
        "../envs/plink.yaml"
    params:
        input   = "{segment}_merged",
        out     = "{segment}_merged_filter"
    log:
        "logs/plink/{segment}_plink_filter.log"
    benchmark:
        "benchmarks/plink/{segment}_plink_filter.txt"
    shell:
        """
        plink --vcf {input.vcf} --freqx --out plink/{params.out}
        plink --vcf {input.vcf} --remove {input.bad_samples} --make-bed --keep-allele-order --out plink/{params.out} |& tee {log}
        """

rule pre_imputation_check:
    input:
        "plink/{segment}_merged_filter.bim"
    params:
        SITE_1000GENOME
    output:
        temp("plink/{segment}_merged_filter.bim.chr"),
        temp("plink/{segment}_merged_filter.bim.pos"),
        "plink/{segment}_merged_filter.bim.force_allele",
        temp("plink/{segment}_merged_filter.bim.flip")
    log:
        "logs/plink/{segment}_pre_imputation_check.log"
    benchmark:
        "benchmarks/plink/{segment}_pre_imputation_check.txt"
    script:
        "../scripts/pre_imputation_check.py"

rule plink_clean_up:
    input:
        "plink/{segment}_merged_filter.bim.chr",
        "plink/{segment}_merged_filter.bim.pos",
        "plink/{segment}_merged_filter.bim.force_allele",
        "plink/{segment}_merged_filter.bim.flip",
        bed="plink/{segment}_merged_filter.bed",
        bim="plink/{segment}_merged_filter.bim",
        fam="plink/{segment}_merged_filter.fam"
    output:
        temp("plink/{segment}_merged_mapped.bim"),
        temp("plink/{segment}_merged_mapped.bed"),
        temp("plink/{segment}_merged_mapped.fam")
    params:
        input = "plink/{segment}_merged_filter",
        out = "plink/{segment}_merged_mapped"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/{segment}_plink_clean_up.log"
    benchmark:
        "benchmarks/plink/{segment}_plink_clean_up.txt"
    shell:
        """
        # remove dublicates
        cut -f 2 {params.input}.bim | sort | uniq -d > plink/{wildcards.segment}_snp.dups
        plink --bfile {params.input}          --exclude       plink/{wildcards.segment}_snp.dups                  --make-bed --out plink/{wildcards.segment}_merged_filter_dub   |& tee -a {log}
        plink --bfile {params.input}          --extract       plink/{wildcards.segment}_merged_filter.bim.freq    --make-bed --out plink/{wildcards.segment}_merged_freq         |& tee -a {log}
        plink --bfile plink/{wildcards.segment}_merged_freq       --extract       plink/{wildcards.segment}_merged_filter.bim.chr     --make-bed --out plink/{wildcards.segment}_merged_extracted    |& tee -a {log}
        plink --bfile plink/{wildcards.segment}_merged_extracted  --flip          plink/{wildcards.segment}_merged_filter.bim.flip    --make-bed --out plink/{wildcards.segment}_merged_flipped      |& tee -a {log}
        plink --bfile plink/{wildcards.segment}_merged_flipped    --update-chr    plink/{wildcards.segment}_merged_filter.bim.chr     --make-bed --out plink/{wildcards.segment}_merged_chroped      |& tee -a {log}
        plink --bfile plink/{wildcards.segment}_merged_chroped    --update-map    plink/{wildcards.segment}_merged_filter.bim.pos     --make-bed --out {params.out}              |& tee -a {log}
        rm plink/{wildcards.segment}_merged_filter_dub.*
        rm plink/{wildcards.segment}_merged_freq.*
        rm plink/{wildcards.segment}_merged_extracted.*
        rm plink/{wildcards.segment}_merged_flipped.*
        rm plink/{wildcards.segment}_merged_chroped.*
        """

rule prepare_vcf:
    input:
        bed="plink/{segment}_merged_mapped.bed",
        bim="plink/{segment}_merged_mapped.bim",
        fam="plink/{segment}_merged_mapped.fam"
    output:
        vcf=temp("vcf/{segment}_merged_mapped_sorted.vcf.gz"),
        temp_vcf=temp("vcf/{segment}_merged_mapped_regions.vcf.gz"),
        temp_vcf_csi=temp("vcf/{segment}_merged_mapped_regions.vcf.gz.csi"),
        bed=temp("plink/{segment}_merged_mapped_alleled.bed"),
        bim=temp("plink/{segment}_merged_mapped_alleled.bim"),
        fam=temp("plink/{segment}_merged_mapped_alleled.fam")
    params:
        input   = "plink/{segment}_merged_mapped",
        vcf     = "vcf/{segment}_merged_mapped_sorted.vcf.gz"
    conda:
         "../envs/bcf_plink.yaml"
    log:
        plink="logs/plink/{segment}_prepare_vcf.log",
        vcf="logs/vcf/{segment}_prepare_vcf.log"
    benchmark:
        "benchmarks/plink/{segment}_prepare_vcf.txt"
    shell:
        """
        plink --bfile {params.input} --a1-allele plink/{wildcards.segment}_merged_filter.bim.force_allele --make-bed --out plink/{wildcards.segment}_merged_mapped_alleled |& tee -a {log.plink}
        plink --bfile plink/{wildcards.segment}_merged_mapped_alleled --keep-allele-order --output-chr M --export vcf bgz --out vcf/{wildcards.segment}_merged_mapped_clean |& tee -a {log.vcf}
        mkdir vcf/temp_{wildcards.segment}
        bcftools sort -T vcf/temp_{wildcards.segment} vcf/{wildcards.segment}_merged_mapped_clean.vcf.gz -O z -o {output.temp_vcf} |& tee -a {log.vcf}
        bcftools index -f {output.temp_vcf} |& tee -a {log.vcf}
        # need to check output for the potential issues
        bcftools view {output.temp_vcf} --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -O z -o {output.vcf}
        bcftools norm --check-ref e -f {GRCH37_FASTA} vcf/{wildcards.segment}_merged_mapped_sorted.vcf.gz -O u -o /dev/null |& tee -a {log.vcf}
        bcftools index -f {output.vcf} | tee -a {log.vcf}
        """
