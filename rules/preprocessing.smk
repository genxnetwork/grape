NUM_BATCHES = config['num_batches']
BATCHES = list(range(1,int(NUM_BATCHES) + 1))
with open('pass_batches.list', 'w') as list:
    for b in BATCHES:
        list.write(f'batch{b}\n')



if NUM_BATCHES > 1:
    rule get_lists:
        input:
            vcf='input.vcf.gz'
        output:
            temp(expand('vcf/batch{i}.txt',i=BATCHES))
        params:
            num_batches=NUM_BATCHES
        conda:
            '../envs/bcftools.yaml'
        shell:
            '''
            bcftools query --list-samples input.vcf.gz >> vcf/samples.txt
            total_lines=$(wc -l < vcf/samples.txt)
            num_files={params.num_batches}
            ((lines_per_file = (total_lines + num_files - 1) / num_files))
            split -l $lines_per_file vcf/samples.txt vcf/batch --additional-suffix=.txt --numeric-suffixes=1
            for file in $(find vcf -name 'batch0[1-9].txt')
            do
                new=$(echo "$file" | sed "s/0//g")
                mv "$file" "$new"
            done
            '''


    rule split_into_batches:
        input:
            vcf='input.vcf.gz',
            samples='vcf/{batch}.txt'
        output:
            vcf='vcf/{batch}.vcf.gz'
        conda:
            '../envs/bcftools.yaml'
        shell:
            '''
            bcftools view -S {input.samples} {input.vcf} -O z -o {output.vcf} --force-samples
            '''
else:
    rule copy_batch:
        input:
            vcf = 'input.vcf.gz'
        output:
            vcf=temp('vcf/batch1.vcf.gz')
        shell:
            '''
            cp {input.vcf} vcf/batch1.vcf.gz
            '''


if need_remove_imputation:
    rule remove_imputation:
        input:
            vcf='vcf/{batch}.vcf.gz'
        output:
            vcf=temp('vcf/{batch}_imputation_removed.vcf.gz')
        log: 'logs/vcf/remove_imputation{batch}.log'
        script: '../scripts/remove_imputation.py'
else:
    rule copy_vcf:
        input:
            vcf='vcf/{batch}.vcf.gz'
        output:
            vcf=temp('vcf/{batch}_imputation_removed.vcf.gz')
        shell:
            '''
                cp {input.vcf} {output.vcf}
            '''


rule recode_vcf:
    input:
        vcf='vcf/{batch}_imputation_removed.vcf.gz'
    output:
        vcf=temp('vcf/{batch}_merged_recoded.vcf.gz')
    conda:
        '../envs/bcftools.yaml'
    shell:
        '''
        rm -f chr_name_conv.txt
        for i in {{1..22}} X Y XY MT
        do
            echo "chr$i $i" >> chr_name_conv.txt
        done
        bcftools annotate --rename-chrs chr_name_conv.txt {input.vcf} | bcftools view -m2 -M2 -v snps -t "^X,Y,XY,MT" -O z -o {output.vcf}
        '''


if assembly == 'hg38':
    rule liftover:
        input:
            vcf='vcf/{batch}_imputation_removed.vcf.gz'
        output:
            vcf=temp('vcf/{batch}_merged_lifted.vcf.gz')
        conda:
            '../envs/liftover.yaml'
        log:
            'logs/liftover/liftover{batch}.log'
        params:
            mem_gb=_mem_gb_for_ram_hungry_jobs()
        resources:
            mem_mb=_mem_gb_for_ram_hungry_jobs() * 1024
        shell:
            '''
               java -Xmx{params.mem_gb}g -jar /picard/picard.jar LiftoverVcf WARN_ON_MISSING_CONTIG=true MAX_RECORDS_IN_RAM=5000 I={input.vcf} O={output.vcf} CHAIN={LIFT_CHAIN} REJECT=vcf/chr{batch}_rejected.vcf.gz R={GRCH37_FASTA} |& tee -a {log}
            '''
else:
    rule copy_liftover:
        input:
            vcf='vcf/{batch}_imputation_removed.vcf.gz'
        output:
            vcf=temp('vcf/{batch}_merged_lifted.vcf.gz')
        shell:
            '''
                cp {input.vcf} {output.vcf}
            '''


rule recode_snp_ids:
    input:
        vcf='vcf/{batch}_merged_lifted.vcf.gz'
    output:
        vcf=temp('vcf/{batch}_merged_lifted_id.vcf.gz')
    conda:
        '../envs/bcftools.yaml'
    shell:
        '''
            bcftools annotate --set-id "%CHROM:%POS:%REF:%FIRST_ALT" {input.vcf} -O z -o {output.vcf}
        '''


include: '../rules/filter.smk'

if need_phase:
    include: '../rules/phasing.smk'
else:
    rule copy_phase:
        input:
            vcf='vcf/{batch}_merged_mapped_sorted.vcf.gz'
        output:
            vcf=temp('phase/{batch}_merged_phased.vcf.gz')
        shell:
            '''
                cp {input.vcf} {output.vcf}
            '''


if need_imputation:
    include: '../rules/imputation.smk'
else:
    if NUM_BATCHES > 1:
        vcf_output = 'preprocessed/{batch}_data.vcf.gz'
        vcf_input = 'phase/{batch}_merged_phased.vcf.gz'
    else:
        vcf_output = 'preprocessed/data.vcf.gz'
        vcf_input = 'phase/batch1_merged_phased.vcf.gz'

    checkpoint copy_imputation:
        input:
            vcf=vcf_input
        output:
            vcf=vcf_output
        shell:
            '''
               cp {input.vcf} {output.vcf}
            '''


if NUM_BATCHES > 1:
    checkpoint convert_mapped_to_plink:
        input:
            vcf='preprocessed/{batch}_data.vcf.gz'
        output:
            bed='preprocessed/{batch}_data.bed',
            fam='preprocessed/{batch}_data.fam',
            bim='preprocessed/{batch}_data.bim'
        params:
            out='preprocessed/{batch}_data'
        conda:
            '../envs/plink.yaml'
        log:
            'logs/plink/convert_mapped_to_plink{batch}.log'
        benchmark:
            'benchmarks/plink/convert_mapped_to_plink{batch}.txt'
        shell:
            '''
            plink --vcf {input} --make-bed --out {params.out} |& tee {log}
            '''


    rule ibis_mapping:
        input:
            bim=rules.convert_mapped_to_plink.output['bim']
        params:
            genetic_map_GRCh37=expand(GENETIC_MAP_GRCH37,chrom=CHROMOSOMES)
        conda:
            '../envs/ibis.yaml'
        output:
            'preprocessed/{batch}_data_mapped.bim'
        log:
            'logs/ibis/run_ibis_mapping{batch}.log'
        benchmark:
            'benchmarks/ibis/run_ibis_mapping{batch}.txt'
        shell:
            '''
            (add-map-plink.pl -cm {input.bim} {params.genetic_map_GRCh37}> {output}) |& tee -a {log}
            '''


    rule index_vcf:
        input:
            batches_vcf='preprocessed/{batch}_data.vcf.gz'
        output:
            batches_vcf_index=temp('preprocessed/{batch}_data.vcf.gz.csi')
        conda:
            '../envs/bcftools.yaml'
        shell:
            '''
            bcftools index -f {input.batches_vcf}
            '''

    def get_merge_vcf_input(wildcards):
        with open('pass_batches.list', 'r') as list:
            batches_left = []
            for line in list:
                batches_left.append(line.strip('\n'))
        index = ['preprocessed/{s}_data.vcf.gz.csi'.format(s=batch) for batch in batches_left]
        vcf = ['preprocessed/{s}_data.vcf.gz'.format(s=batch) for batch in batches_left]
        return index + vcf

    rule merge_vcf:
        input:
            get_merge_vcf_input
        output:
            vcf='preprocessed/data.vcf.gz'
        threads:
            workflow.cores
        params:
            batches_vcf=expand('batch{s}_data.vcf.gz',s=BATCHES),
            vcf='data.vcf.gz'
        conda:
            '../envs/bcftools.yaml'
        shell:
            '''
            rm complete_vcf_list.txt || true
            for FILE in {input}
            do
                if [[ $FILE == *.gz ]]
                then
                    echo $FILE >> complete_vcf_list.txt
                fi
            done
            bcftools merge --threads {threads} --file-list complete_vcf_list.txt --force-samples -O z -o {output.vcf}
            rm complete_vcf_list.txt
            '''


    def get_merge_bed_input(wildcards):
        with open('pass_batches.list', 'r') as list:
            batches_left = []
            for line in list:
                batches_left.append(line.strip('\n'))
        bim = ['preprocessed/{s}_data_mapped.bim'.format(s=batch) for batch in batches_left]
        bed = ['preprocessed/{s}_data.bed'.format(s=batch) for batch in batches_left]
        fam = ['preprocessed/{s}_data.fam'.format(s=batch) for batch in batches_left]
        return bed + bim + fam

    rule merge_bed:
        input:
            get_merge_bed_input
        output:
            bed='preprocessed/data.bed',
            fam='preprocessed/data.fam',
            bim_mapped='preprocessed/data_mapped.bim'
        threads:
            workflow.cores
        conda:
            '../envs/plink.yaml'
        shell:
            '''
            for file in $(find preprocessed -name '*_mapped.bim')
            do
                new=$(echo "$file" | sed "s/_mapped//g")
                mv "$file" "$new"
            done

            rm files_list.txt || true
            for file in {input}
            do
                if [[ $file == *.fam ]]
                then
                    echo ${{file%.*}} >> files_list.txt
                fi
            done

            plink --merge-list files_list.txt --make-bed --out preprocessed/data
            mv preprocessed/data.bim preprocessed/data_mapped.bim
            rm files_list.txt
            '''


    rule remove_mapping:
        input:
            bim_mapped = 'preprocessed/data_mapped.bim'
        output:
            bim='preprocessed/data.bim'
        conda:
            '../envs/remove_map.yaml'
        script:
            '../scripts/remove_mapping.py'
else:
    rule single_batch_convert_mapped_to_plink:
        input:
            vcf='preprocessed/data.vcf.gz'
        output:
            bed='preprocessed/data.bed',
            fam='preprocessed/data.fam',
            bim='preprocessed/data.bim'
        params:
            out='preprocessed/data'
        conda:
            '../envs/plink.yaml'
        log:
            'logs/plink/convert_mapped_to_plink_batch1.log'
        benchmark:
            'benchmarks/plink/convert_mapped_to_plink_batch1.txt'
        shell:
            '''
            plink --vcf {input} --make-bed --out {params.out} |& tee {log}
            '''


    rule single_batch_ibis_mapping:
        input:
            bim=rules.single_batch_convert_mapped_to_plink.output['bim']
        params:
            genetic_map_GRCh37=expand(GENETIC_MAP_GRCH37,chrom=CHROMOSOMES)
        conda:
            '../envs/ibis.yaml'
        output:
            'preprocessed/data_mapped.bim'
        log:
            'logs/ibis/run_ibis_mapping_batch1.log'
        benchmark:
            'benchmarks/ibis/run_ibis_mapping_batch1.txt'
        shell:
            '''
            (add-map-plink.pl -cm {input.bim} {params.genetic_map_GRCh37}> {output}) |& tee -a {log}
            '''
