rule phase:
        input:
            vcf="vcf/merged_mapped_sorted.vcf.gz",
            #idx="vcf/merged_mapped_sorted.vcf.gz.csi"
            vcfRef=REF_VCF
        output: temp("phase/chr{chrom}.phased.vcf.gz")
        threads: 1
        singularity:
            "docker://biocontainers/bio-eagle:v2.4.1-1-deb_cv1"
        log:
            "logs/phase/eagle-{chrom}.log"
        benchmark:
            "benchmarks/phase/eagle-{chrom}.txt"
        shell:
            """
            /usr/bin/bio-eagle --vcfRef {input.vcfRef} \
            --numThreads {threads} \
            --vcfTarget {input.vcf}  \
            --geneticMapFile {GENETIC_MAP} \
            --chrom {wildcards.chrom} \
            --vcfOutFormat z \
            --pbwtIters 2 \
            --Kpbwt 20000 \
            --outPrefix phase/chr{wildcards.chrom}.phased |& tee {log}
            """

rule merge_phased:
    input:
        expand("phase/chr{i}.phased.vcf.gz", i=CHROMOSOMES)
    output:
        "phase/merged_phased.vcf.gz"
    params:
        list="vcf/phased.merge.list",
        mode=config["mode"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf/merged_phased.log"
    benchmark:
        "benchmarks/vcf/merged_phased.txt"
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
        if [ -f "background/merged_imputed.vcf.gz" && {params.mode} = "client" ]; then
            mv {output} {output}.client
            bcftools merge --force-samples background/merged_imputed.vcf.gz {output}.client -O z -o {output} |& tee -a {log}
            bcftools index -f {output} |& tee -a {log}
        fi
        """