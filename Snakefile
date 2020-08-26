#configfile: "config.yaml"
# run as snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all

# TODO: add config
# TODO: add reporting https://github.com/tanaes/snarkmark/blob/master/rules/report.rule

CHROMOSOMES     = [str(i) for i in range(1, 23)]
DATASETS        = [i for i in range(1, 41)] # include #40
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']
PLINK           = "plink/*.bed"

rule all:
    input:
        "ersa/ersa.chr1"

rule convert_23andme_to_plink:
    input:
        expand("input/{i}.txt", i=DATASETS)
    output:
        sorted(expand("plink/{i}.{ext}", i=DATASETS, ext=PLINK_FORMATS))
    conda:
        # TODO: separate for plink and bcftools
        "envs/plink.yaml"
    shell:
        """
        for i in {DATASETS}; do
            plink --23file input/$i.txt $i $i i --output-chr M --make-bed --out plink/$i
        done
        """

rule merge_list:
    input:
        rules.convert_23andme_to_plink.output
    output:
        "plink/merge.list"
    shell:
        """
        true > {output} && for i in {DATASETS}; do echo "plink/$i" >> {output}; done
        """

rule merge_to_vcf:
    input: rules.merge_list.output # sorted(expand("plink/{i}.{ext}", i=DATASETS, ext=PLINK_FORMATS))
    output:
        plink_clean = sorted(expand("plink/{i}_clean.{ext}", i=DATASETS, ext=PLINK_FORMATS)),
        vcf = "vcf/merged.vcf.gz"
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
            true > {input} && for i in {DATASETS}; do echo plink/$i\_clean >> {input}; done
            plink --merge-list {input} --output-chr M --export vcf bgz --out vcf/merged
            exit 0
        fi
        exit 1
        """

# TODO: extra step bc we already have plink in vcf/merged
rule convert_to_plink:
    input: rules.merge_to_vcf.output["vcf"]
    output: expand("plink/{i}.{ext}", i="merged", ext=PLINK_FORMATS_EXT)
    params:
        out = "plink/merged"
    conda:
        "envs/plink.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out}
        """

rule plink_filter:
    input: rules.convert_to_plink.output
    output: expand("plink/{i}.{ext}", i="merged_filter", ext=PLINK_FORMATS_EXT)
    conda:
        "envs/plink.yaml"
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
    
rule impute:
    input: rules.phase.output
    output: expand("imputed/chr{i}.imputed.dose.vcf.gz", i=CHROMOSOMES)
    threads: workflow.cores
    singularity:
        "docker://biocontainers/minimac4:v1.0.0-2-deb_cv1"
    shell:
        """
        # TODO: hardcoded for now
        MINIMAC3_M3VCF=/media/Minimac/\$i.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz

        for i in `seq 1 22`; do
            if ! /usr/bin/minimac4 \
            --refHaps /media/Minimac/$i.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
            --haps phase/chr$i.phased.vcf.gz \
            --format GT,GP \
            --prefix imputed/chr$i.imputed \
            --cpus {threads}
            then
                touch imputed/chr$i.imputed.dose.vcf.gz
                continue
            fi
        done
        """

rule imputation_filter:
    input: rules.impute.output
    output: expand("imputed/chr{i}.imputed.dose.pass.vcf.gz", i=CHROMOSOMES)
    threads: workflow.cores
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        FILTER="'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1'"

        for i in `seq 1 22`; do
            if ! bcftools view --threads {threads} -i 'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1' imputed/chr$i.imputed.dose.vcf.gz -v snps -m 2 -M 2 -O z -o imputed/chr$i.imputed.dose.pass.vcf.gz
            then
                touch imputed/chr$i.imputed.dose.pass.vcf.gz
                continue
            fi
        done
        """

rule merge_imputation_filter:
    input:
        rules.imputation_filter.output
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
    output: expand("vcf/imputed_chr{i}.vcf.gz", i=CHROMOSOMES)
    conda:
        "envs/bcftools.yaml"
    threads: workflow.cores
    shell:
        """
        bcftools index -f {input}

        for i in `seq 1 22`; do
            bcftools filter {input} -r $i --threads {threads} -O z -o vcf/imputed_chr$i.vcf.gz;
        done
        """

rule convert_to_hap:
    input: rules.index_and_split.output
    output: expand("hap/imputed_chr{i}.hap", i=CHROMOSOMES)
    conda:
        "envs/bcftools.yaml"
    threads: workflow.cores
    shell:
        """
        for i in `seq 1 22`; do
            bcftools convert vcf/imputed_chr$i.vcf.gz --threads {threads} --hapsample hap/imputed_chr$i && gunzip -f hap/imputed_chr$i.hap.gz;
        done
        """

rule convert_to_ped:
    input: rules.convert_to_hap.output
    output: expand("ped/imputed_chr{i}.ped", i=CHROMOSOMES)
    singularity:
        "docker://alexgenx/germline:stable"
    shell:
        """
        for i in `seq 1 22`; do
            impute_to_ped hap/imputed_chr$i.hap hap/imputed_chr$i.sample ped/imputed_chr$i
        done
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
    output: expand("cm/chr{i}.cm.ped", i=CHROMOSOMES)
    conda:
        "envs/plink.yaml"
    shell:
        """
        for i in `seq 1 22`; do
            plink --file ped/imputed_chr$i --cm-map /media/genetic_map_b37/genetic_map_chr$i\_combined_b37.txt $i --recode --out cm/chr$i.cm
        done
        """

# TODO: does not work as conda do not know why
rule germline:
    input: rules.interpolate.output
    output: expand("germline/chr{i}.germline.match", i=CHROMOSOMES)
    singularity:
        "docker://alexgenx/germline:stable"
    shell:
        """
        for i in `seq 1 22`; do
            germline -input ped/imputed_chr$i.ped cm/chr$i.cm.map -homoz -min_m 2.5 -err_hom 2 -err_het 1 -output germline/chr$i.germline;
            # TODO: germline returns some length in BP instead of cM - clean up is needed
            grep -v MB germline/chr$i.germline.match > germline/chr$i.germline.match.clean && mv germline/chr$i.germline.match.clean germline/chr$i.germline.match
        done
        """

rule ersa:
    # TODO: look for conda/container/pip
    # TODO: failed for chr8 bc of MB instead of Cm but why??? need to double check in interpolate
    # TODO: grep -v MB germline/chr8.germline.match > germline/chr8.germline.match.clean
    input: rules.germline.output
    output: expand("ersa/ersa.chr{i}", i=CHROMOSOMES)
    singularity:
        "docker://alexgenx/ersa:stable"
    shell:
        """
        ERSA_L=13.7
        ERSA_TH=3.197
        ERSA_T=2.5

        for i in `seq 1 22`; do
            ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH germline/chr$i.germline.match -o ersa/ersa.chr$i
        done
        """
