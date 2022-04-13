CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', ]
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'nosex']


rule impute:
    input:
        rules.phase.output,
        refHaps=REF_HAPS
    output: temp("imputed/{batch}_chr{chrom}.imputed.dose.vcf.gz")
    log:
        "logs/impute/{batch}_minimac4-{chrom}.log"
    benchmark:
        "benchmarks/impute/{batch}_minimac4-{chrom}.txt"
    shell:
        """
            minimac4 \
            --refHaps {input.refHaps} \
            --haps phase/{wildcards.batch}_chr{wildcards.chrom}.phased.vcf.gz \
            --format GT,GP \
            --prefix imputed/{wildcards.batch}_chr{wildcards.chrom}.imputed \
            --minRatio 0.01 |& tee {log}
        """

rule imputation_filter:
    input: rules.impute.output
    output: temp("imputed/{batch}_chr{chrom}.imputed.dose.pass.vcf.gz")
    # TODO: because "The option is currently used only for the compression of the output stream"
    # threads: workflow.cores
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/impute/{batch}_imputation_filter-{chrom}.log"
    benchmark:
        "benchmarks/impute/{batch}_imputation_filter-{chrom}.txt"
    shell:
        """
        FILTER="'strlen(REF)=1 & strlen(ALT)=1'"

        bcftools view -i 'strlen(REF)=1 & strlen(ALT)=1' imputed/{wildcards.batch}_chr{wildcards.chrom}.imputed.dose.vcf.gz -v snps -m 2 -M 2 -O z -o imputed/{wildcards.batch}_chr{wildcards.chrom}.imputed.dose.pass.vcf.gz |& tee {log}
        """


imputed_dose_pass = ["chr{i}.imputed.dose.pass.vcf.gz".format(i=chr) for chr in CHROMOSOMES]
if NUM_BATCHES > 1:
    imputed_dose_pass_batch = ["imputed/{batch}_" + line for line in imputed_dose_pass]
    data_vcf_gz = ["preprocessed/{batch}_data.vcf.gz"]
    card = "{batch}"
else:
    imputed_dose_pass_batch = ["imputed/batch1_" + line for line in imputed_dose_pass]
    data_vcf_gz = ["preprocessed/data.vcf.gz"]
    card = "batch1"
rule merge_imputation_filter:
    input:
        imputed_dose_pass_batch
        # TODO: wildcard violation
        # rules.imputation_filter.output
    output:
        data_vcf_gz
    params:
        list="vcf/"+ card +"_imputed.merge.list",
        mode=config["mode"],
        merged_imputed = "background/"+ card +"_merged_imputed.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf/"+ card +"_merge_imputation_filter.log"
    benchmark:
        "benchmarks/vcf/"+ card +"_merge_imputation_filter.txt"
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

        bcftools concat -f {params.list} -O z -o {output} |& tee -a {log}
        bcftools index -f {output} |& tee -a {log}

        # check if there is a background data and merge it
        if [ -f "{params.merged_imputed}" && {params.mode} = "client" ]; then
            mv {output} {output}.client
            bcftools merge --force-samples {params.merged_imputed} {output}.client -O z -o {output} |& tee -a {log}
            bcftools index -f {output} |& tee -a {log}
        fi
        """


rule convert_imputed_to_plink:
    input:
        rules.merge_imputation_filter.output
    output:
        "plink/{batch}_merged_imputed.bed",
        "plink/{batch}_merged_imputed.bim",
        "plink/{batch}_merged_imputed.fam"
    params:
        out = "plink/{batch}_merged_imputed"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/{batch}_convert_imputed_to_plink.log"
    benchmark:
        "benchmarks/plink/{batch}_convert_imputed_to_plink.txt"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out} |& tee {log}
        """

# no need it bc it was done earlier in merge_imputation_filter
rule merge_convert_imputed_to_plink:
    input:
        rules.merge_imputation_filter.output
    output:
        "plink/{batch}_merged_imputed.bed",
        "plink/{batch}_merged_imputed.bim",
        "plink/{batch}_merged_imputed.fam"
    params:
        background  = "background/{batch}_merged_imputed",
        out         = "plink/client/{batch}_merged_imputed"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/{batch}_convert_imputed_to_plink.log"
    benchmark:
        "benchmarks/plink/{batch}_convert_imputed_to_plink.txt"
    shell:
        """
        # please mind a merge step in merge_imputation_filter for germline
        plink --vcf {input} --make-bed --out {params.out} | tee {log}
        plink --bfile {params.background} --bmerge {params.out} --make-bed --out plink/{batch}_merged_imputed |& tee {log}
        """
