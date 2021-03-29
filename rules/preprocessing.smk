
rule recode_vcf:
    input: vcf='input.vcf'
    output: vcf = "vcf/merged_recoded.vcf.gz"
    log: "logs/plink/recode_vcf.log"
    conda: "../envs/plink.yaml"
    shell: "plink --vcf {input.vcf} --chr 1-22 --snps-only just-acgt --output-chr M --not-chr XY,MT --export vcf bgz --out vcf/merged_recoded |& tee {log}"


if need_remove_imputation:
    rule remove_imputation:
        input:
            vcf=rules.recode_vcf.input['vcf']
        output:
            vcf='vcf/imputation_removed.vcf.gz'
        log: "logs/vcf/remove_imputation.log"
        script: '../../scripts/remove_imputation.py'
else:
    rule copy_vcf:
        input:
            vcf=rules.recode_vcf.input['vcf']
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