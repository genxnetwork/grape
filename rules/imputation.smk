CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']


# TODO: extra step bc we already have plink in vcf/merged
rule convert_to_plink:
    input:
         vcf = config["vcf_file"]
    output: expand("plink/{i}.{ext}", i="merged", ext=PLINK_FORMATS_EXT)
    params:
        out = "plink/merged"
    conda:
        "../envs/plink.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out}
        """

rule plink_filter:
    input: rules.convert_to_plink.output
    output: expand("plink/{i}.{ext}", i="merged_filter", ext=PLINK_FORMATS_EXT)
    conda:
        "../envs/plink.yaml"
    params:
        input = "merged",
        out = "merged_filter"
    shell:
        """
        # TODO: the old parameter --maf 1e-50 is too low, back to the "default" 0.05
        plink --bfile plink/{params.input} --geno 0.1 --maf 0.05 --hwe 0 --make-bed --keep-allele-order --out plink/{params.out}
        """

rule pre_imputation_check:
    input: "plink/merged_filter.bim"
    output:
        "plink/merged_filter.bim.chr",
        "plink/merged_filter.bim.pos",
        "plink/merged_filter.bim.force_allele",
        "plink/merged_filter.bim.flip"
    script:
        "../scripts/pre_imputation_check.py"

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
    shell:
        """
        plink --bfile {params.input}         --extract       plink/merged_filter.bim.chr     --make-bed --out plink/merged_extracted
        plink --bfile plink/merged_extracted --flip          plink/merged_filter.bim.flip    --make-bed --out plink/merged_flipped
        plink --bfile plink/merged_flipped   --update-chr    plink/merged_filter.bim.chr     --make-bed --out plink/merged_chroped
        plink --bfile plink/merged_chroped   --update-map    plink/merged_filter.bim.pos     --make-bed --out {params.out}
        """

rule prepare_vcf:
    input: "plink/merged_mapped.bim"
    output:
        "vcf/merged_mapped_sorted.vcf.gz"
    params:
        input = "plink/merged_mapped"
    conda:
         "../envs/snakemake.yaml"
    shell:
        """
        GRCh37_fasta=/media/human_g1k_v37.fasta

        plink --bfile {params.input} --a1-allele plink/merged_filter.bim.force_allele --make-bed --out plink/merged_mapped_alleled
        plink --bfile plink/merged_mapped_alleled --keep-allele-order --output-chr M --export vcf bgz --out vcf/merged_mapped_clean
        bcftools sort vcf/merged_mapped_clean.vcf.gz -O z -o vcf/merged_mapped_sorted.vcf.gz
        # need to check output for the potential issues
        # bcftools norm --check-ref e -f $GRCh37_fasta merged_mapped_sorted.vcf.gz -O u -o /dev/null
        bcftools index -f vcf/merged_mapped_sorted.vcf.gz
        """

rule phase:
    input: "vcf/merged_mapped_sorted.vcf.gz"
    output: expand("phase/chr{i}.phased.vcf.gz", i=CHROMOSOMES)
    threads: workflow.cores
    singularity:
        "docker://biocontainers/bio-eagle:v2.4.1-1-deb_cv1"
    shell:
        """
        GENOTYPE_1000GENOME=/media/1000genome/bcf/1000genome_chr\$i.bcf
        # TODO: not sure what to do with that
        GENETIC_MAP=/media/tables/genetic_map_hg19_withX.txt.gz

        # TODO: do some check before because of ERROR: Target and ref have too few matching SNPs (M = 0)
        for i in `seq 1 22`; do
            if ! /usr/bin/bio-eagle --vcfRef  /media/1000genome/bcf/1000genome_chr$i.bcf \
            --numThreads {threads} \
            --vcfTarget {input}  \
            --geneticMapFile $GENETIC_MAP \
            --chrom $i \
            --vcfOutFormat z \
            --outPrefix phase/chr$i.phased
            then
                touch phase/chr$i.phased.vcf.gz
                continue
            fi
        done
        """

rule merge_phased:
    input:
        rules.phase.output
    output:
        "vcf/merged_phased.vcf.gz"
    params:
        list="vcf/phased.merge.list"
    conda:
        "../envs/bcftools.yaml"
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
        bcftools concat -f {params.list} -O z -o {output}
        """

rule convert_phased_to_plink:
    input: rules.merge_phased.output
    output: expand("plink/merged_phased.{ext}", ext=PLINK_FORMATS)
    params:
        out = "plink/merged_phased"
    conda:
        "../envs/plink.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out}
        """