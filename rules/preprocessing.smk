NUM_BATCHES        = config["num_batches"]

rule get_lists:
    input:
        vcf="input.vcf.gz"
    output:
        temp(expand("vcf/segment{i}.txt",i=list(range(1,int(NUM_BATCHES) + 1))))
    params:
        num_batches=NUM_BATCHES
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools query --list-samples input.vcf.gz >> vcf/samples.txt
        total_lines=$(wc -l < vcf/samples.txt)
        num_files={params.num_batches}
        ((lines_per_file = (total_lines + num_files - 1) / num_files))
        split -l $lines_per_file vcf/samples.txt vcf/segment --additional-suffix=.txt --numeric-suffixes=1
        for file in $(find vcf -name "segment0[1-9].txt") 
        do 
            new=$(echo "$file" | sed 's/0//g')
            mv "$file" "$new"
        done
        """


rule split_into_segments:
    input:
        vcf="input.vcf.gz",
        samples="vcf/{segment}.txt"
    output:
        vcf=temp("vcf/{segment}.vcf.gz")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -S {input.samples} {input.vcf} -O z -o {output.vcf} --force-samples
        """

rule recode_vcf:
    input: vcf='vcf/{segment}.vcf.gz'
    output: vcf=temp('vcf/{segment}_merged_recoded.vcf.gz')
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
            vcf=temp('vcf/{segment}_imputation_removed.vcf.gz')
        log: "logs/vcf/remove_imputation{segment}.log"
        script: '../scripts/remove_imputation.py'
else:
    rule copy_vcf:
        input:
            vcf=rules.recode_vcf.output['vcf']
        output:
            vcf=temp('vcf/{segment}_imputation_removed.vcf.gz')
        shell:
            """
                cp {input.vcf} {output.vcf}
            """

if assembly == "hg38":
    rule liftover:
        input:
            vcf='vcf/{segment}_imputation_removed.vcf.gz'
        output:
            vcf=temp("vcf/{segment}_merged_lifted.vcf.gz")
        singularity:
            "docker://genxnetwork/picard:stable"
        log:
            "logs/liftover/liftover{segment}.log"
        params:
            mem_gb=_mem_gb_for_ram_hungry_jobs()
        resources:
            mem_mb=_mem_gb_for_ram_hungry_jobs() * 1024
        shell:
            """
               java -Xmx{params.mem_gb}g -jar /picard/picard.jar LiftoverVcf WARN_ON_MISSING_CONTIG=true MAX_RECORDS_IN_RAM=5000 I={input.vcf} O={output.vcf} CHAIN={LIFT_CHAIN} REJECT=vcf/chr{segment}_rejected.vcf.gz R={GRCH37_FASTA} |& tee -a {log}
            """
else:
    rule copy_liftover:
        input:
            vcf='vcf/{segment}_imputation_removed.vcf.gz'
        output:
            vcf=temp("vcf/{segment}_merged_lifted.vcf.gz")
        shell:
            """
                cp {input.vcf} {output.vcf}
            """

rule recode_snp_ids:
    input:
        vcf="vcf/{segment}_merged_lifted.vcf.gz"
    output:
        vcf=temp("vcf/{segment}_merged_lifted_id.vcf.gz")
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
            vcf="vcf/{segment}_merged_mapped_sorted.vcf.gz"
        output:
            vcf=temp("phase/{segment}_merged_phased.vcf.gz")
        shell:
            """
                cp {input.vcf} {output.vcf}
            """

if need_imputation:
    include: "../rules/imputation.smk"
else:
    rule copy_imputation:
        input:
            vcf="phase/{segment}_merged_phased.vcf.gz"
        output:
            vcf=temp("preprocessed/{segment}_data.vcf.gz")
        shell:
            """
               cp {input.vcf} {output.vcf}
            """

rule convert_mapped_to_plink:
    input:
        vcf="preprocessed/{segment}_data.vcf.gz"
    output:
        bed=temp("preprocessed/{segment}_data.bed"),
        fam=temp("preprocessed/{segment}_data.fam"),
        bim=temp("preprocessed/{segment}_data.bim")
    params:
        out="preprocessed/{segment}_data"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/plink/convert_mapped_to_plink{segment}.log"
    benchmark:
        "benchmarks/plink/convert_mapped_to_plink{segment}.txt"
    shell:
        """
        plink --vcf {input} --make-bed --out {params.out} |& tee {log}
        """

rule ibis_mapping:
    input:
        bim=rules.convert_mapped_to_plink.output['bim']
    params:
        input="preprocessed/segement{segment}_data",
        genetic_map_GRCh37=expand(GENETIC_MAP_GRCH37,chrom=CHROMOSOMES)
    conda:
        "../envs/ibis.yaml"
    output:
        "preprocessed/{segment}_data_mapped.bim"
    log:
        "logs/ibis/run_ibis_mapping{segment}.log"
    benchmark:
        "benchmarks/ibis/run_ibis_mapping{segment}.txt"
    shell:
        """
        (add-map-plink.pl -cm {input.bim} {params.genetic_map_GRCh37}> {output}) |& tee -a {log}
        """

#seg = glob_wildcards("preprocessed/{segment}_data.bed")
rule index_vcf:
    input:
        segments_vcf="preprocessed/{segment}_data.vcf.gz"
    output:
        segments_vcf_index = temp("preprocessed/{segment}_data.vcf.gz.csi")
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools index -f {input.segments_vcf}
        """
rule merge_vcf:
    input:
        segments_vcf_index = expand("preprocessed/segment{s}_data.vcf.gz.csi",s=list(range(1,int(NUM_BATCHES) + 1))),
        segments_vcf = expand("preprocessed/segment{s}_data.vcf.gz",s=list(range(1,int(NUM_BATCHES) + 1)))
    output:
        vcf = "preprocessed/data.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
            bcftools merge --merge id {input.segments_vcf} -O z -o {output.vcf}
        """
rule merge_bed:
    input:
        segments_bim_mapped=expand("preprocessed/segment{s}_data_mapped.bim",s=list(range(1,int(NUM_BATCHES) + 1))),
        segments_bed=expand("preprocessed/segment{s}_data.bed",s=list(range(1,int(NUM_BATCHES) + 1))),
        segments_fam=expand("preprocessed/segment{s}_data.fam",s=list(range(1,int(NUM_BATCHES) + 1)))
    output:
        bed="preprocessed/data.bed",
        fam="preprocessed/data.fam",
        bim_mapped="preprocessed/data_mapped.bim"
    params:
        seg=expand("preprocessed/segment{s}_data",s=list(range(1,int(NUM_BATCHES) + 1)))
    conda:
        "../envs/plink.yaml"
    shell:
        """
        for file in $(find preprocessed -name "*_mapped.bim") 
        do 
            new=$(echo "$file" | sed 's/_mapped//g')
            mv "$file" "$new"
        done
        
        for file in {params.seg}
        do
            echo "$file.bed $file.bim $file.fam" >> files_list.txt
        done
        
        plink --merge-list files_list.txt --make-bed --out preprocessed/data 
        mv preprocessed/data.bim preprocessed/data_mapped.bim
        """

rule remove_mapping:
    input:
        bim_mapped = "preprocessed/data_mapped.bim"
    output:
        bim = "preprocessed/data.bim"
    conda:
        "../envs/remove_map.yaml"
    script:
        "../scripts/remove_mapping.py"
