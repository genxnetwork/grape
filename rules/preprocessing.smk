rule recode_vcf:
    input: vcf = 'input.vcf.gz'
    output: vcf = 'vcf/merged_recoded.vcf.gz'
    conda: "../envs/bcftools.yaml"
    shell:
        """
        rm -f chr_name_conv.txt
        for i in {{1..22}} X Y XY MT
        do
            echo "chr$i $i" >> chr_name_conv.txt
        done
        
        bcftools annotate --rename-chrs chr_name_conv.txt {input.vcf} | bcftools view -m2 -M2 -v snps -t "^X,Y,XY,MT" -O z -o {output.vcf} 
        """

if need_remove_imputation:
    rule remove_imputation:
        input:
            vcf=rules.recode_vcf.output['vcf']
        output:
            vcf='vcf/imputation_removed.vcf.gz'
        log: "logs/vcf/remove_imputation.log"
        script: '../scripts/remove_imputation.py'
else:
    rule copy_vcf:
        input:
            vcf=rules.recode_vcf.output['vcf']
        output:
            vcf='vcf/imputation_removed.vcf.gz'
        shell:
            """
                cp {input.vcf} {output.vcf}
            """

if assembly == "hg38":
    rule liftover:
        input:
            vcf='vcf/imputation_removed.vcf.gz'
        output:
            vcf="vcf/merged_lifted.vcf.gz"
        singularity:
            "docker://genxnetwork/picard:stable"
        log:
            "logs/liftover/liftover.log"
        params:
            mem_gb = _mem_gb_for_ram_hungry_jobs()
        resources:
            mem_mb = _mem_gb_for_ram_hungry_jobs() * 1024
        shell:
             """
                java -Xmx{params.mem_gb}g -jar /picard/picard.jar LiftoverVcf WARN_ON_MISSING_CONTIG=true MAX_RECORDS_IN_RAM=25000 I={input.vcf} O={output.vcf} CHAIN={LIFT_CHAIN} REJECT=vcf/rejected.vcf.gz R={GRCH37_FASTA} |& tee -a {log}
             """
else:
    rule copy_liftover:
        input:
            vcf='vcf/imputation_removed.vcf.gz'
        output:
            vcf="vcf/merged_lifted.vcf.gz"
        shell:
            """
                cp {input.vcf} {output.vcf}
            """

rule recode_snp_ids:
    input:
        vcf="vcf/merged_lifted.vcf.gz"
    output:
        vcf="vcf/merged_lifted_id.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
            bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' {input.vcf} -O z -o {output.vcf}
        """

include: "../rules/filter.smk"

if need_phase:
    include: "../rules/phasing.smk"
else:
    rule copy_phase:
        input:
            vcf="vcf/merged_mapped_sorted.vcf.gz"
        output:
            vcf="phase/merged_phased.vcf.gz"
        shell:
            """
                cp {input.vcf} {output.vcf}
            """

if need_imputation:
    include: "../rules/imputation.smk"
else:
    rule copy_imputation:
        input:
            vcf="phase/merged_phased.vcf.gz"
        output:
            vcf="preprocessed/data.vcf.gz"
        shell:
             """
                cp {input.vcf} {output.vcf}
             """

rule convert_mapped_to_plink:
    input:
        vcf="preprocessed/data.vcf.gz"
    output:
        bed="preprocessed/data.bed",
        fam="preprocessed/data.fam",
        bim="preprocessed/data.bim"
    params:
        out = "preprocessed/data"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/convert_mapped_to_plink.log"
    benchmark:
        "benchmarks/plink/convert_mapped_to_plink.txt"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out} |& tee {log}
        """

rule ibis_mapping:
    input:
        bim=rules.convert_mapped_to_plink.output['bim']
    params:
        input = "preprocessed/data",
        genetic_map_GRCh37 = expand(GENETIC_MAP_GRCH37, chrom=CHROMOSOMES)
    conda:
        "../envs/ibis.yaml"
    output:
        "preprocessed/data_mapped.bim"
    log:
        "logs/ibis/run_ibis_mapping.log"
    benchmark:
        "benchmarks/ibis/run_ibis_mapping.txt"
    shell:
        """
        (add-map-plink.pl -cm {input.bim} {params.genetic_map_GRCh37}> {output}) |& tee -a {log}
        """
