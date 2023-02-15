rule select_bad_samples:
    input:
        vcf='vcf/{batch}.vcf.gz'
    output:
        bad_samples='vcf/{batch}_lifted_vcf.badsamples',
        report='results/{batch}_bad_samples_report.tsv',
        outliers='results/{batch}_outliers.list',
        stats='stats/{batch}_lifted_vcf.txt',
        no_outliers_vcf=temp('vcf/{batch}_no_outliers.vcf.gz')
    log: 'logs/vcf/{batch}_select_bad_samples.log'
    params:
        samples='vcf/{batch}_merged_lifted.vcf.samples',
        missing_samples = config['missing_samples'],
        alt_hom_samples = config['alt_hom_samples'],
        het_samples = config['het_samples'],
        iqr_alpha = config['iqr_alpha'],
        psc='stats/{batch}_lifted_vcf.psc',
        keep_samples='stats/{batch}_keep_samples.list',

    conda:
        'evaluation'
    script:
        '../scripts/select_bad_samples.py'


rule plink_filter:
    input:
        vcf='vcf/{batch}_merged_lifted_id.vcf.gz',
        bad_samples=rules.select_bad_samples.output.bad_samples
    output:
        bed = temp('plink/{batch}_merged_filter.bed'),
        bim = temp('plink/{batch}_merged_filter.bim'),
        fam = temp('plink/{batch}_merged_filter.fam')
    conda:
        'plink'
    params:
        input   = '{batch}_merged',
        out     = '{batch}_merged_filter',
        batch   = '{batch}'
    log:
        'logs/plink/{batch}_plink_filter.log'
    benchmark:
        'benchmarks/plink/{batch}_plink_filter.txt'
    script:
        '../scripts/plink_filter.py'


rule pre_imputation_check:
    input:
        'plink/{batch}_merged_filter.bim'
    params:
        SITE_1000GENOME
    output:
        temp('plink/{batch}_merged_filter.bim.chr'),
        temp('plink/{batch}_merged_filter.bim.pos'),
        'plink/{batch}_merged_filter.bim.force_allele',
        temp('plink/{batch}_merged_filter.bim.flip')
    log:
        'logs/plink/{batch}_pre_imputation_check.log'
    benchmark:
        'benchmarks/plink/{batch}_pre_imputation_check.txt'
    script:
        '../scripts/pre_imputation_check.py'

rule plink_clean_up:
    input:
        'plink/{batch}_merged_filter.bim.chr',
        'plink/{batch}_merged_filter.bim.pos',
        'plink/{batch}_merged_filter.bim.force_allele',
        'plink/{batch}_merged_filter.bim.flip',
        bed='plink/{batch}_merged_filter.bed',
        bim='plink/{batch}_merged_filter.bim',
        fam='plink/{batch}_merged_filter.fam'
    output:
        temp('plink/{batch}_merged_mapped.bim'),
        temp('plink/{batch}_merged_mapped.bed'),
        temp('plink/{batch}_merged_mapped.fam')
    params:
        input = 'plink/{batch}_merged_filter',
        out = 'plink/{batch}_merged_mapped'
    conda:
        'plink'
    log:
        'logs/plink/{batch}_plink_clean_up.log'
    benchmark:
        'benchmarks/plink/{batch}_plink_clean_up.txt'
    shell:
        """
        # remove dublicates
        cut -f 2 {params.input}.bim | sort | uniq -d > plink/{wildcards.batch}_snp.dups
        plink --bfile {params.input}          --exclude       plink/{wildcards.batch}_snp.dups                  --make-bed --out plink/{wildcards.batch}_merged_filter_dub   |& tee -a {log}
        plink --bfile {params.input}          --extract       plink/{wildcards.batch}_merged_filter.bim.freq    --make-bed --out plink/{wildcards.batch}_merged_freq         |& tee -a {log}
        plink --bfile plink/{wildcards.batch}_merged_freq       --extract       plink/{wildcards.batch}_merged_filter.bim.chr     --make-bed --out plink/{wildcards.batch}_merged_extracted    |& tee -a {log}
        plink --bfile plink/{wildcards.batch}_merged_extracted  --flip          plink/{wildcards.batch}_merged_filter.bim.flip    --make-bed --out plink/{wildcards.batch}_merged_flipped      |& tee -a {log}
        plink --bfile plink/{wildcards.batch}_merged_flipped    --update-chr    plink/{wildcards.batch}_merged_filter.bim.chr     --make-bed --out plink/{wildcards.batch}_merged_chroped      |& tee -a {log}
        plink --bfile plink/{wildcards.batch}_merged_chroped    --update-map    plink/{wildcards.batch}_merged_filter.bim.pos     --make-bed --out {params.out}              |& tee -a {log}
        rm plink/{wildcards.batch}_merged_filter_dub.*
        rm plink/{wildcards.batch}_merged_freq.*
        rm plink/{wildcards.batch}_merged_extracted.*
        rm plink/{wildcards.batch}_merged_flipped.*
        rm plink/{wildcards.batch}_merged_chroped.*
        """

rule prepare_vcf:
    input:
        bed='plink/{batch}_merged_mapped.bed',
        bim='plink/{batch}_merged_mapped.bim',
        fam='plink/{batch}_merged_mapped.fam'
    output:
        vcf=temp('vcf/{batch}_merged_mapped_sorted.vcf.gz'),
        temp_vcf=temp('vcf/{batch}_merged_mapped_regions.vcf.gz'),
        temp_vcf_csi=temp('vcf/{batch}_merged_mapped_regions.vcf.gz.csi'),
        bed=temp('plink/{batch}_merged_mapped_alleled.bed'),
        bim=temp('plink/{batch}_merged_mapped_alleled.bim'),
        fam=temp('plink/{batch}_merged_mapped_alleled.fam')
    params:
        input   = 'plink/{batch}_merged_mapped',
        vcf     = 'vcf/{batch}_merged_mapped_sorted.vcf.gz'
    conda:
         'bcf_plink'
    log:
        plink='logs/plink/{batch}_prepare_vcf.log',
        vcf='logs/vcf/{batch}_prepare_vcf.log'
    benchmark:
        'benchmarks/plink/{batch}_prepare_vcf.txt'
    shell:
        """
        plink --bfile {params.input} --a1-allele plink/{wildcards.batch}_merged_filter.bim.force_allele --make-bed --out plink/{wildcards.batch}_merged_mapped_alleled |& tee -a {log.plink}
        plink --bfile plink/{wildcards.batch}_merged_mapped_alleled --keep-allele-order --output-chr M --export vcf bgz --out vcf/{wildcards.batch}_merged_mapped_clean |& tee -a {log.vcf}
        mkdir vcf/temp_{wildcards.batch}
        bcftools sort -T vcf/temp_{wildcards.batch} vcf/{wildcards.batch}_merged_mapped_clean.vcf.gz -O z -o {output.temp_vcf} |& tee -a {log.vcf}
        bcftools index -f {output.temp_vcf} |& tee -a {log.vcf}
        # need to check output for the potential issues
        bcftools view {output.temp_vcf} --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -O z -o {output.vcf}
        bcftools norm --check-ref e -f {GRCH37_FASTA} vcf/{wildcards.batch}_merged_mapped_sorted.vcf.gz -O u -o /dev/null |& tee -a {log.vcf}
        bcftools index -f {output.vcf} | tee -a {log.vcf}
        """
