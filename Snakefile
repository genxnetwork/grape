#configfile: "config.yaml"

CHROMOSOMES     = [str(i) for i in range(1, 23)]
DATASETS        = [i for i in range(1,40)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']
PLINK           = "plink/*.bed"

rule all:
    input:
        "vcf/merged.vcf.gz"

rule convert_23andme_to_plink:
    input:
        "input/{sample}.txt"
    output:
        sorted(expand("plink/{{sample}}.{ext}", ext=PLINK_FORMATS))
    conda:
        "envs/pipeline.yaml"
    benchmark:
        "benchmarks/convert_23andme_to_plink_{sample}.txt"
    shell:
        """
        plink --23file {input} {wildcards.sample} {wildcards.sample} i --output-chr M --make-bed --out plink/{wildcards.sample}
        """

rule merge_list:
    input:
        sorted(expand("plink/{i}.{ext}", i=DATASETS, ext=PLINK_FORMATS))
    output:
        "merge.list"
    params:
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
        "envs/pipeline.yaml"
    shell:
        """
        set +e
        plink --merge-list {input} --output-chr M --export vcf bgz --out merged
        exitcode=$?
        echo $exitcode

        if [[ -f "merged.missnp" ]]; then
            for i in `cat {input}`;
                do plink --bfile $i --exclude merged.missnp --make-bed --out $i\_clean;
            done
            true > {input} && for i in {DATASETS}; do echo plink/$i\_clean >> {input}; done
            plink --merge-list {input} --output-chr M --export vcf bgz --out vcf/merged
            exit 0
        fi
        exit 1
        """

rule convert_to_plink:
    input: rules.merge_to_vcf.output["vcf"]
    output: expand("plink/{i}.{ext}", i="merged", ext=PLINK_FORMATS_EXT)
    params:
        out = "plink/merged"
    conda:
        "envs/pipeline.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out}
        """

rule plink_filter:
    input: rules.convert_to_plink.output
    output: expand("plink/{i}.{ext}", i="merged_filter", ext=PLINK_FORMATS_EXT)
    conda:
        "envs/pipeline.yaml"
    params:
        input = "merged",
        out = "merged_filter"
    shell:
        """
        plink --bfile plink/{params.input} --geno 0.1 --maf 1e-50 --hwe 0 --make-bed --keep-allele-order --out plink/{params.out}
        """

rule pre_imputation_check:
    input: "plink/merged_filter.bim"
    output:
        "plink/merged_filter.chr",
        "plink/merged_filter.pos",
        "plink/merged_filter.force_allele",
        "plink/merged_filter.flip"
    script:
        "pre_imputation_check.py"


rule plink_clean_up:
    input:
        "plink/merged_filter.bim"
    output:
        expand("{i}.{ext}", i="plink/merged_mapped", ext=PLINK_FORMATS)
    params:
        input = "plink/merged_filter",
        out = "plink/merged_mapped"
    conda:
        "envs/pipeline.yaml"
    shell:
        """
        plink --bfile {params.input}                --extract       plink/merged_filter.chr     --make-bed --out plink/merged_mapped_extracted
        plink --bfile plink/merged_mapped_extracted --flip          plink/merged_filter.flip    --make-bed --out plink/merged_mapped_flipped
        plink --bfile plink/merged_mapped_flipped   --update-chr    plink/merged_filter.chr     --make-bed --out plink/merged_mapped_chroped
        plink --bfile plink/merged_mapped_chroped   --update-map    plink/merged_filter.pos     --make-bed --out {params.out}
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

        plink --bfile {params.input} --a1-allele plink/merged_filter.force_allele --make-bed --out plink/merged_mapped_alleled
        plink --bfile plink/merged_mapped_alleled --keep-allele-order --output-chr M --export vcf bgz --out vcf/merged_mapped_clean
        bcftools sort vcf/merged_mapped_clean.vcf.gz -O z -o vcf/merged_mapped_sorted.vcf.gz
        # need to check output for the potential issues
        # bcftools norm --check-ref e -f $GRCh37_fasta merged_mapped_sorted.vcf.gz -O u -o /dev/null
        bcftools index -f vcf/merged_mapped_sorted.vcf.gz
        """

rule phase:
    input: "vcf/merged_mapped_sorted.vcf.gz"
    output: expand("phase/chr{i}.phased.vcf.gz", i=CHROMOSOMES)
    # conda:
    #     "envs/pipeline.yaml"
    container:
        "docker://indraniel/eagle-241"
    shell:
        """
        EAGLE=/bin/eagle
        GENOTYPE_1000GENOME=/etc/genome
        GENETIC_MAP=/etc/map
        
        echo parallel -k -j 100% \
        eagle --vcfRef  $GENOTYPE_1000GENOME \
        --vcfTarget {input}  \
        --geneticMapFile $GENETIC_MAP \
        --chrom {{1}} \
        --vcfOutFormat z \
        --outPrefix phase/chr{{1}}.phased \
        ::: `seq 1 22`
        """
    
rule impute:
    input: rules.phase.output
    output: expand("chr{i}.imputed.dose.vcf.gz", i=CHROMOSOMES)
    shell:
        """
        MINIMAC3=/bin/minimac
        MINIMAC3_M3VCF=/etc/minimac

        echo parallel -k -j 100% \
        $MINIMAC3 --refHaps $MINIMAC3_M3VCF \
        --haps {{1}}.chr{{1}}.phased.vcf.gz \
        --format GT,GP \
        --prefix chr{{1}}.imputed \
        ::: `ls *.phased.vcf.gz* | awk -F'[_.]' '{{print $1}}'`

        for i in {output}; do
            touch $i;
        done
        """

rule imputation_filter:
    input: rules.impute.output
    output: expand("chr{i}.dose.pass.vcf.gz", i=CHROMOSOMES)
    shell:
        """
        BCFTOOLS=/bin/bcftools
        
        for i in `ls *imputed.dose.vcf.gz`;
            do echo $BCFTOOLS view -i 'R2>0.3 & strlen(REF)=1 & strlen(ALT)=1' $i -v snps -m 2 -M 2 -O z -o ${{i%.*.*}}.pass.vcf.gz;
        done

        for i in {output}; do
            touch $i;
        done        
        """

rule merge_imputation_filter:
    input: rules.imputation_filter.output
    output: "merged_imputed.vcf.gz"
    shell:
        """
        IMPUTELIST=impute.list
        BCFTOOLS=/bin/bcftools

        ls *dose.pass.vcf.gz > $IMPUTELIST
        echo $BCFTOOLS concat -f $IMPUTELIST -O z -o merged_imputed.vcf.gz

        touch merged_imputed.vcf.gz
        """

rule convert_imputed_to_plink:
    input: rules.merge_imputation_filter.output
    output: expand("{i}.{ext}", i="merged_imputed", ext=PLINK_FORMATS)
    params:
        out = "merged_imputed"
    shell:
        """
        PLINK=/bin/plink

        echo $PLINK --vcf {input} --make-bed --out {params.out}
        touch {params.out}.bed && touch {params.out}.bim && touch {params.out}.fam
        """

rule run_king:
    input: rules.convert_imputed_to_plink.output
    output: "merged_imputed_king"
    params:
        input = "merged_imputed"
    shell:
        """
        KING=/bin/king
        KING_DEGREE=4

        echo $KING -b {params.input}.bed --ibdseg --degree $KING_DEGREE --prefix {params.input}\_king
        touch {params.input}\_king
        """

rule index_and_split:
    input: rules.merge_imputation_filter.output
    output: expand("imputed_chr{i}.vcf.gz", i=CHROMOSOMES)
    shell:
        """
        BCFTOOLS=/bin/bcftools
        
        echo $BCFTOOLS index -f {input}

        echo parallel -k -j 100% \
        $BCFTOOLS filter {input} -r {{1}} -O z -o imputed_chr{{1}}.vcf.gz \
        ::: `seq 1 22`

        for i in {output}; do
            touch $i;
        done
        """

rule convert_to_hap:
    input: rules.merge_imputation_filter.output
    output: expand("chr{i}.hap", i=CHROMOSOMES)
    shell:
        """
        BCFTOOLS=/bin/bcftools

        echo parallel -k -j 100% \
        $BCFTOOLS convert imputed_chr{{1}}.vcf.gz --hapsample imputed_chr{{1}} '&&' gunzip -f chr{{1}}.hap.gz \
        ::: `seq 1 22`

        for i in {output}; do
            touch $i;
        done
        """

rule convert_to_ped:
    input: rules.convert_to_hap.output
    output: expand("chr{i}.ped", i=CHROMOSOMES)
    shell:
        """
        IMPUTE_2_PED=/bin/impute_to_ped

        echo parallel -k -j 100% \
        $IMPUTE_2_PED chr{{1}}.hap chr{{1}}.sample chr{{1}} \
        ::: `seq 1 22`

        for i in {output}; do
            touch $i;
        done
        """

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
    input: rules.split_map.output
    output: expand("chr{i}.cm", i=CHROMOSOMES)
    shell:
        """
        PLINK=/bin/plink

        echo parallel -k -j 100% \
        $PLINK --file chr{{1}} --cm-map genetic_map_b37/genetic_map_chr{{1}}_combined_b37.txt {{1}} --recode --out chr{{1}}.cm \
        ::: `seq 1 22`

        for i in {output}; do
            touch $i;
        done
        """

rule germline:
    input: rules.interpolate.output
    output: expand("chr{i}.germline", i=CHROMOSOMES)
    shell:
        """
        GERMLINE=/bin/germline

        echo parallel -k -j 100% \
        $GERMLINE -input chr{{}}.ped chr{{}}.cm.map -homoz -min_m 2.5 -err_hom 2 -err_het 1 -output chr{{}}.germline \
        ::: `seq 1 22`

        for i in {output}; do
            touch $i;
        done
        """

rule ersa:
    input: rules.germline.output
    output: "ersa"
    shell:
        """
        ERSA_L=13.7
        ERSA_TH=3.197
        ERSA_T=2.5

        ERSA=/bin/ersa

        echo parallel -k -j 100% \
        $ERSA --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH chr{{}}.germline -o ersa{{}} \
        ::: `seq 1 22`

        touch {output}
        """
