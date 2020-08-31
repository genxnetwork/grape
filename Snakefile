# run as "snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all"
# TODO: add reporting https://github.com/tanaes/snarkmark/blob/master/rules/report.rule

import pandas as pd

configfile: "config.yaml"

samples_file    = config["samples_file"]
SAMPLES         = pd.read_table(samples_file).set_index("name", drop=False)
SAMPLES_PATH    = SAMPLES.path.values.tolist()
SAMPLES_NAME    = SAMPLES.name.values.tolist()

CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']

def get_samples_name(wildcards):
    return SAMPLES.loc[int(wildcards.sample), "name"] # mind the int index

def get_samples_path(wildcards):
    return SAMPLES.loc[int(wildcards.sample), "path"] # mind the int index

rule all:
    input:
        # "results/relatives.tsv"

rule convert_23andme_to_plink:
    input:
        get_samples_path
    output:
        expand("plink/{{sample}}.{ext}", ext=PLINK_FORMATS)
    params:
        get_samples_name
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --23file {input} {params} {params} i --output-chr M --make-bed --out plink/{params}
        """

rule merge_list:
    input:
        expand("plink/{sample}.{ext}", sample=SAMPLES_NAME, ext=PLINK_FORMATS)
    output:
        "plink/merge.list"
    shell:
        """
        true > {output} && for i in {SAMPLES_NAME}; do echo "plink/$i" >> {output}; done
        """

rule merge_to_vcf:
    input: rules.merge_list.output
    output:
        plink_clean     = expand("plink/{i}_clean.{ext}", i=SAMPLES_NAME, ext=PLINK_FORMATS),
        plink_merged    = expand("vcf/merged.{ext}", ext=PLINK_FORMATS),
        vcf             = "vcf/merged.vcf.gz",
        merge_list      = "plink/merge_clean.list"
    conda:
        "envs/plink.yaml"
    shell:
        """
        set +e
        plink --merge-list {input} --output-chr M --export vcf bgz --out vcf/merged
        exitcode=$?

        if [[ -f "vcf/merged.missnp" ]]; then
            for i in `cat {input}`;
                do plink --bfile $i --exclude vcf/merged.missnp --make-bed --out $i\_clean;
            done
            true > {output.merge_list} && for i in {SAMPLES_NAME}; do echo plink/$i\_clean >> {output.merge_list}; done
            plink --merge-list {output.merge_list} --output-chr M --export vcf bgz --out vcf/merged
            exit 0
        fi
        exit 1
        """

rule plink_filter:
    input: rules.merge_to_vcf.output["vcf"]
    output: expand("plink/merged_filter.{ext}", ext=PLINK_FORMATS)
    conda:
        "envs/plink.yaml"
    params:
        input   = "merged",
        out     = "merged_filter"
    shell:
        """
        # TODO: the old parameter --maf 1e-50 is too low, back to the "default" 0.05
        plink --bfile vcf/{params.input} --geno 0.1 --maf 0.05 --hwe 0 --make-bed --keep-allele-order --out plink/{params.out}
        """

rule pre_imputation_check:
    input: "plink/merged_filter.bim"
    output:
        "plink/merged_filter.bim.chr",
        "plink/merged_filter.bim.pos",
        "plink/merged_filter.bim.force_allele",
        "plink/merged_filter.bim.flip"
    script:
        "pre_imputation_check.py"

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
        "envs/plink.yaml"
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
        "envs/pipeline.yaml"
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
    output: "phase/chr{chrom}.phased.vcf.gz"
    threads: workflow.cores
    singularity:
        "docker://biocontainers/bio-eagle:v2.4.1-1-deb_cv1"
    shell:
        """
        GENETIC_MAP=/media/tables/genetic_map_hg19_withX.txt.gz

        /usr/bin/bio-eagle --vcfRef  /media/1000genome/bcf/1000genome_chr{wildcards.chrom}.bcf \
        --numThreads {threads} \
        --vcfTarget {input}  \
        --geneticMapFile $GENETIC_MAP \
        --chrom {wildcards.chrom} \
        --vcfOutFormat z \
        --outPrefix phase/chr{wildcards.chrom}.phased
        """
    
rule impute:
    input: rules.phase.output
    output: "imputed/chr{chrom}.imputed.dose.vcf.gz"
    threads: workflow.cores
    singularity:
        "docker://biocontainers/minimac4:v1.0.0-2-deb_cv1"
    shell:
        """
        /usr/bin/minimac4 \
        --refHaps /media/Minimac/{wildcards.chrom}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
        --haps phase/chr{wildcards.chrom}.phased.vcf.gz \
        --format GT,GP \
        --prefix imputed/chr{wildcards.chrom}.imputed \
        --cpus {threads}
        """

rule imputation_filter:
    input: rules.impute.output
    output: "imputed/chr{chrom}.imputed.dose.pass.vcf.gz"
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        FILTER="'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1'"

        bcftools view -i 'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1' imputed/chr{wildcards.chrom}.imputed.dose.vcf.gz -v snps -m 2 -M 2 -O z -o imputed/chr{wildcards.chrom}.imputed.dose.pass.vcf.gz
        """

rule merge_imputation_filter:
    input:
        expand("imputed/chr{i}.imputed.dose.pass.vcf.gz", i=CHROMOSOMES)
        # TODO: wildcard violation
        # rules.imputation_filter.output
    output:
        "vcf/merged_imputed.vcf.gz"
    params:
        list="vcf/imputed.merge.list"
    conda:
        "envs/bcftools.yaml"
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
        bcftools index -f {output}
        """

rule convert_imputed_to_plink:
    input: rules.merge_imputation_filter.output
    output: expand("plink/{i}.{ext}", i="merged_imputed", ext=PLINK_FORMATS)
    params:
        out = "plink/merged_imputed"
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out}
        """

rule run_king:
    input: rules.convert_imputed_to_plink.output
    output: "king/merged_imputed_king.seg"
    params:
        input = "plink/merged_imputed",
        out = "king/merged_imputed_king"
    threads: workflow.cores
    singularity:
        "docker://lifebitai/king:latest"
    shell:
        """
        # TODO: add cores
        KING_DEGREE=4

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out}
        """

rule index_and_split:
    input: rules.merge_imputation_filter.output
    output: "vcf/imputed_chr{chrom}.vcf.gz"
    conda:
        "envs/bcftools.yaml"
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    shell:
        """
        bcftools filter {input} -r {wildcards.chrom} -O z -o vcf/imputed_chr{wildcards.chrom}.vcf.gz;
        """

rule convert_to_hap:
    input: rules.index_and_split.output
    output: "hap/imputed_chr{chrom}.hap"
    conda:
        "envs/bcftools.yaml"
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    shell:
        """
        bcftools convert vcf/imputed_chr{wildcards.chrom}.vcf.gz --hapsample hap/imputed_chr{wildcards.chrom} && gunzip -f hap/imputed_chr{wildcards.chrom}.hap.gz;
        """

rule convert_to_ped:
    input: rules.convert_to_hap.output
    output: "ped/imputed_chr{chrom}.ped"
    singularity:
        "docker://alexgenx/germline:stable"
    shell:
        """
        impute_to_ped hap/imputed_chr{wildcards.chrom}.hap hap/imputed_chr{wildcards.chrom}.sample ped/imputed_chr{wildcards.chrom}
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
    input: rules.convert_to_hap.output
    output: "cm/chr{chrom}.cm.ped"
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --file ped/imputed_chr{wildcards.chrom} --cm-map /media/genetic_map_b37/genetic_map_chr{wildcards.chrom}\_combined_b37.txt {wildcards.chrom} --recode --out cm/chr{wildcards.chrom}.cm
        """

rule germline:
    input: rules.interpolate.output
    output: "germline/chr{chrom}.germline.match"
    singularity:
        "docker://alexgenx/germline:stable"
    shell:
        """
        germline -input ped/imputed_chr{wildcards.chrom}.ped cm/chr{wildcards.chrom}.cm.map -homoz -min_m 2.5 -err_hom 2 -err_het 1 -output germline/chr{wildcards.chrom}.germline
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
    conda: "envs/ersa_params.yaml"
    script: "scripts/estimate_ersa_params.py"

rule merge_matches:
    input:
         expand("germline/chr{chrom}.germline.match", chrom=CHROMOSOMES)
    output:
          "germline/all.tsv"
    shell:
         "cat germline/*.match > {output}"


rule ersa:
    # TODO: look for conda/container/pip
    # TODO: failed for chr8 bc of MB instead of Cm but why??? need to double check in interpolate
    # TODO: grep -v MB germline/chr8.germline.match > germline/chr8.germline.match.clean
    input:
        germline=rules.merge_matches.output,
        estimated=rules.ersa_params.output
    output: "ersa/relatives.tsv"
    singularity:
        "docker://alexgenx/ersa:stable"
    shell:
        """
        #ERSA_L=13.7
        #ERSA_TH=3.197
        #ERSA_T=2.5
        source {input.estimated}
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.germline} -o {output}
        """

rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        germline=rules.merge_matches.output,
        ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "envs/evaluation.yaml"
    script: "scripts/evaluate.py"
