rule index_input:
    input: "preprocessed/data.vcf.gz"
    output: "preprocessed/data.vcf.gz.csi"
    conda: "../envs/bcftools.yaml"
    shell: "bcftools index -f {input}"

rule convert_mapped_to_plink:
    input:
        vcf="preprocessed/data.vcf.gz",
        index="preprocessed/data.vcf.gz.csi"
    output:
        plink=expand("plink/{i}.{ext}", i="merged_ibis", ext=PLINK_FORMATS),
        vcf122=temp("vcf/merged_mapped_sorted_22.vcf.gz")
    params:
        out = "plink/merged_ibis"
    conda:
        "../envs/bcf_plink.yaml"
    log:
        "logs/plink/convert_mapped_to_plink.log"
    benchmark:
        "benchmarks/plink/convert_mapped_to_plink.txt"
    shell:
        """
        # leave only chr1..22 because we need to map it later
        bcftools view {input} --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -O z -o {output.vcf122}
        plink --vcf {output.vcf122} --make-bed --out {params.out} |& tee {log}
        """

rule run_king:
    input: rules.convert_mapped_to_plink.output
    output:
        king="king/data.seg",
        kinship="king/data.kin",
        kinship0="king/data.kin0",
        segments="king/data.segments.gz"
    params:
        input="plink/merged_ibis",
        out="king/data",
        kin="king/data"
    threads: workflow.cores
    singularity:
        "docker://lifebitai/king:latest"
    log:
        "logs/king/run_king.log"
    benchmark:
        "benchmarks/king/run_king.txt"
    shell:
        """
        KING_DEGREE=3

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out} |& tee {log}
        king -b {params.input}.bed --cpus {threads} --kinship --degree 4 --prefix {params.kin} |& tee -a {log}
        
        # we need at least an empty file for the downstream analysis
        if [ ! -f "{output.king}" ]; then
            touch {output.king}
        fi
        if [ ! -f "{output.kinship}" ]; then
            touch {output.kinship}
        fi
        if [ ! -f "{output.kinship0}" ]; then
            touch {output.kinship0}
        fi
        if [ ! -f "{output.segments}" ]; then
            touch {output.segments}
        fi
        """

rule ibis_mapping:
    input:
        rules.convert_mapped_to_plink.output
    params:
        input = "plink/merged_ibis"
    singularity:
        "docker://genxnetwork/ibis:stable"
    output:
        "plink/merged_ibis_mapped.bim"
    log:
        "logs/ibis/run_ibis_mapping.log"
    benchmark:
        "benchmarks/ibis/run_ibis_mapping.txt"
    shell:
        """
        (add-map-plink.pl -cm {params.input}.bim {GENETIC_MAP_GRCH37}> plink/merged_ibis_mapped.bim) |& tee -a {log}
        """

rule ibis:
    input:
        rules.ibis_mapping.output
    params:
        input = "plink/merged_ibis"
    singularity:
        "docker://genxnetwork/ibis:stable"
    output:
        ibd     = "ibis/merged_ibis.seg"
    log:
        "logs/ibis/run_ibis.log"
    benchmark:
        "benchmarks/ibis/run_ibis.txt"
    threads: workflow.cores
    shell:
        """
        ibis {params.input}.bed {input} {params.input}.fam -t {threads} -mt 560 -mL 7 -ibd2 -mL2 3 -hbd -f ibis/merged_ibis |& tee -a {log}

        # cat {output.ibd} | awk '{{sub(":", "_", $1); sub(":", "_", $2); print $1, $1 "\t" $2, $2 "\t" $3 "\t" $4, $5 "\t" 0, 0 "\t" $10 "\t" $9 "\t" "cM" "\t" 0 "\t" 0 "\t" 0}};' > {output}
        """

rule transform_ibis_segments:
    input:
        ibd=rules.ibis.output.ibd
    output:
        germline = "ibd/merged_ibd.tsv"
    log:
        "logs/ibis/transform_ibis_segments.log"
    conda:
        "../envs/evaluation.yaml"
    script:
        "../scripts/transform_ibis_segments.py"

rule split_map:
    input:
        bim = "plink/merged_ibis_mapped.bim"
    output: expand("cm/chr{chrom}.cm.map", chrom=CHROMOSOMES)
    params:
        cm_dir='cm'
    conda:
        "../envs/evaluation.yaml"
    script:
        "../scripts/split_map.py"


rule ersa:
    input:
        ibd=rules.transform_ibis_segments.output['germline'],
        cm=expand("cm/chr{chrom}.cm.map", chrom=CHROMOSOMES) # it does not really need it, just to invoke split_map
    output:
        "ersa/relatives.tsv"
    singularity:
        "docker://genxnetwork/ersa:stable"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    shell:
        """
        ERSA_L=0.5 # the average number of IBD segments in population
        ERSA_TH=6.00 # the average length of IBD segment in population
        ERSA_T=7.0 # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -ci --alpha 0.01 --dmax 14 -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output}  |& tee {log}
        """


rule merge_king_ersa:
    input:
        king=rules.run_king.output['king'],
        king_segments=rules.run_king.output['segments'],
        ibd=rules.transform_ibis_segments.output['germline'],
        ersa=rules.ersa.output[0],
        kinship=rules.run_king.output['kinship'],
        kinship0=rules.run_king.output['kinship0']
    params:
        cm_dir='cm'
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    log: "logs/merge/merge-king-ersa.log"
    script: "../scripts/merge_king_ersa.py"
