rule phase:
        input:
            vcf='vcf/{batch}_merged_mapped_sorted.vcf.gz',
            #idx='vcf/merged_mapped_sorted.vcf.gz.csi'
            vcfRef=REF_VCF
        output: temp('phase/{batch}_chr{chrom}.phased.vcf.gz')
        log:
            'logs/phase/{batch}_eagle-{chrom}.log'
        benchmark:
            'benchmarks/phase/{batch}_eagle-{chrom}.txt'
        shell:
            '''
            eagle --vcfRef {input.vcfRef} \
            --vcfTarget {input.vcf}  \
            --geneticMapFile {GENETIC_MAP} \
            --chrom {wildcards.chrom} \
            --vcfOutFormat z \
            --pbwtIters 2 \
            --Kpbwt 20000 \
            --outPrefix phase/{wildcards.batch}_chr{wildcards.chrom}.phased |& tee {log}
            '''


phase = ['chr{i}.phased.vcf.gz'.format(i=chr) for chr in CHROMOSOMES]
phase_batch = [ 'phase/{batch}_' + line for line in phase]
rule merge_phased:
    input:
        phase_batch
    output:
        'phase/{batch}_merged_phased.vcf.gz'
    params:
        list='vcf/{batch}_phased.merge.list',
        mode=config['mode']
    conda:
        '../envs/bcftools.yaml'
    log:
        'logs/vcf/{batch}_merged_phased.log'
    benchmark:
        'benchmarks/vcf/{batch}_merged_phased.txt'
    shell:
        '''
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
        if [ -f "background/{wildcards.batch}_merged_imputed.vcf.gz" && {params.mode} = "client" ]; then
            mv {output} {output}.client
            bcftools merge --force-samples background/{wildcards.batch}_merged_imputed.vcf.gz {output}.client -O z -o {output} |& tee -a {log}
            bcftools index -f {output} |& tee -a {log}
        fi
        '''
