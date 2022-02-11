rule run_king:
    input: rules.convert_mapped_to_plink.output
    output:
        king="king/data.seg",
        kinship="king/data.kin",
        kinship0="king/data.kin0",
        segments="king/data.segments.gz"
    params:
        input = "preprocessed/data",
        out = "king/data",
        kin = "king/data"
    threads: workflow.cores
    conda:
        "../envs/king.yaml"
    log:
        "logs/king/run_king.log"
    benchmark:
        "benchmarks/king/run_king.txt"
    shell:
        """
        KING_DEGREE=3

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree $KING_DEGREE --prefix {params.out} |& tee {log}
        king -b {params.input}.bed --cpus {threads} --kinship --degree $KING_DEGREE --prefix {params.kin} |& tee -a {log}

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

rule index_input:
    input: "preprocessed/data.vcf.gz"
    output: "preprocessed/data.vcf.gz.csi"
    conda: "../envs/bcftools.yaml"
    shell: "bcftools index -f {input}"

rule index_and_split:
    input:
        vcf="preprocessed/data.vcf.gz",
        index="preprocessed/data.vcf.gz.csi"
    output: "vcf/imputed_chr{chrom}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    # TODO: because "The option is currently used only for the compression of the output stream"
    log:
        "logs/vcf/index_and_split-{chrom}.log"
    benchmark:
        "benchmarks/vcf/index_and_split-{chrom}.txt"
    shell:
        """
        bcftools filter -r {wildcards.chrom} {input.vcf} | bcftools norm --rm-dup none -O z -o vcf/imputed_chr{wildcards.chrom}.vcf.gz |& tee {log}
        """


rule vcf_to_ped:
    # --a2-allele <uncompressed VCF filename> 4 3 '#'
    input:
        vcf=rules.index_and_split.output
    output:
        ped="ped/imputed_chr{chrom}.ped",
        map_file="ped/imputed_chr{chrom}.map"
    params:
        zarr="zarr/imputed_chr{chrom}.zarr"
    conda:
        "../envs/vcf_to_ped.yaml"
    log:
        "logs/ped/vcf_to_ped-{chrom}.log"
    benchmark:
        "benchmarks/vcf/vcf_to_ped-{chrom}.txt"
    script:
        "../scripts/vcf_to_ped.py"


rule interpolate:
    input:
        rules.vcf_to_ped.output,
        cmmap=CMMAP
    output: "cm/chr{chrom}.cm.map"
    conda:
        "../envs/plink.yaml"
    log:
        "logs/cm/interpolate-{chrom}.log"
    benchmark:
        "benchmarks/cm/interpolate-{chrom}.txt"
    shell:
        """
        plink --file ped/imputed_chr{wildcards.chrom} --keep-allele-order --cm-map {input.cmmap} {wildcards.chrom} --recode --out cm/chr{wildcards.chrom}.cm |& tee {log}
        """

rule germline:
    input: rules.interpolate.output
    output: "germline/chr{chrom}.germline.match"
    log:
        "logs/germline/germline-{chrom}.log"
    benchmark:
        "benchmarks/germline/germline-{chrom}.txt"
    shell:
        """
        germline -input ped/imputed_chr{wildcards.chrom}.ped cm/chr{wildcards.chrom}.cm.map -min_m 2.5 -err_hom 2 -err_het 1 -output germline/chr{wildcards.chrom}.germline |& tee {log}
        # TODO: germline returns some length in BP instead of cM - clean up is needed
        set +e
        grep -v MB germline/chr{wildcards.chrom}.germline.match > germline/chr{wildcards.chrom}.germline.match.clean && mv germline/chr{wildcards.chrom}.germline.match.clean germline/chr{wildcards.chrom}.germline.match
        set -e
        """

rule merge_matches:
    input:
         expand("germline/chr{chrom}.germline.match", chrom=CHROMOSOMES)
    output:
          "germline/all.tsv"
    shell:
         "cat germline/*.match > {output}"

rule merge_ibd_segments:
    input:
        germline=rules.merge_matches.output[0]
    params:
        cm_dir='cm',
        merge_gap='0.6',
        use_true_ibd=use_simulated_ibd,
        true_ibd='pedsim/simulated/data.seg' # it is in the params because in the case of true data we do not have this information
    output:
        ibd='ibd/merged_ibd.tsv'
    conda: "../envs/evaluation.yaml"
    script:
        '../scripts/merge_ibd.py'

rule ersa:
    input:
        ibd=rules.merge_ibd_segments.output['ibd']
    output:
        "ersa/relatives.tsv"
    conda:
        "../envs/ersa.yaml"
    log:
        "logs/ersa/ersa.log"
    benchmark:
        "benchmarks/ersa/ersa.txt"
    shell:
        """
        ERSA_L=2.0 # the average number of IBD segments in population
        ERSA_TH=1.5 # the average length of IBD segment
        ERSA_T=1.0 # min length of segment to be considered in segment aggregation
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.ibd} -o {output} |& tee {log}
        """


rule merge_king_ersa:
    input:
        king=rules.run_king.output['king'],
        king_segments=rules.run_king.output['segments'],
        ibd=rules.merge_ibd_segments.output['ibd'],
        ersa=rules.ersa.output[0],
        kinship=rules.run_king.output['kinship'],
        kinship0=rules.run_king.output['kinship0']
    params:
        cm_dir='cm'
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    log: "logs/merge/merge-king-ersa.log"
    script: "../scripts/merge_king_ersa.py"
