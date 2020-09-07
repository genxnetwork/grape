CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']

rule run_king:
    input: rules.convert_phased_to_plink.output
    output: "king/merged_imputed_king.seg"
    params:
        input = "plink/merged_phased",
        out = "king/merged_imputed_king"
    threads: workflow.cores
    singularity:
        "docker://lifebitai/king:latest"
    shell:
        """
        # TODO: add cores
        KING_DEGREE=4

        king -b {params.input}.bed --cpus {threads} --ibdseg --degree KING_DEGREE --prefix {params.out}
        """

rule index_and_split:
    input: rules.merge_phased.output
    output: expand("vcf/imputed_chr{i}.vcf.gz", i=CHROMOSOMES)
    conda:
        "../envs/bcftools.yaml"
    threads: workflow.cores
    shell:
        """
        bcftools index -f {input}

        for i in `seq 1 22`; do
            bcftools filter {input} -r $i --threads {threads} -O z -o vcf/imputed_chr$i.vcf.gz;
        done
        """

rule convert_to_hap:
    input: rules.index_and_split.output
    output: expand("hap/imputed_chr{i}.hap", i=CHROMOSOMES)
    conda:
        "../envs/bcftools.yaml"
    threads: workflow.cores
    shell:
        """
        for i in `seq 1 22`; do
            bcftools convert vcf/imputed_chr$i.vcf.gz --threads {threads} --hapsample hap/imputed_chr$i && gunzip -f hap/imputed_chr$i.hap.gz;
        done
        """

rule convert_to_ped:
    input: rules.convert_to_hap.output
    output: expand("ped/imputed_chr{i}.ped", i=CHROMOSOMES)
    singularity:
        "docker://alexgenx/germline:stable"
    shell:
        """
        for i in `seq 1 22`; do
            impute_to_ped hap/imputed_chr$i.hap hap/imputed_chr$i.sample ped/imputed_chr$i
        done
        """

# TODO: skip this step for now
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
    input: rules.convert_to_ped.output
    output: expand("cm/chr{i}.cm.ped", i=CHROMOSOMES)
    conda:
        "../envs/plink.yaml"
    shell:
        """
        for i in `seq 1 22`; do
            plink --file ped/imputed_chr$i --cm-map /media/genetic_map_b37/genetic_map_chr$i\_combined_b37.txt $i --recode --out cm/chr$i.cm
        done
        """

rule germline:
    input: rules.interpolate.output
    output: expand("germline/chr{i}.germline.match", i=CHROMOSOMES)
    singularity:
        "docker://alexgenx/germline:stable"
    shell:
        """
        for i in `seq 1 22`; do
            germline -input ped/imputed_chr$i.ped cm/chr$i.cm.map -homoz -min_m 2.5 -err_hom 4 -err_het 2 -output germline/chr$i.germline;
            # TODO: germline returns some length in BP instead of cM - clean up is needed
            grep -v MB germline/chr$i.germline.match > germline/chr$i.germline.match.clean && mv germline/chr$i.germline.match.clean germline/chr$i.germline.match
        done
        """

rule ersa_params:
    input:
        king=rules.run_king.output,
        germline=rules.germline.output
    output: "ersa/params.txt"
    conda: "../envs/ersa_params.yaml"
    script: "../scripts/estimate_ersa_params.py"

rule merge_matches:
    input:
         expand("germline/chr{chrom}.germline.match", chrom=CHROMOSOMES)
    output:
          "germline/all.tsv"
    shell:
         "cat germline/*.match > {output}"


rule ersa:
    # TODO: look for conda/container/pip
    # TODO: failed for chr8 bc of MB instead of Cm but why??? need to double check in interpolate
    # TODO: grep -v MB germline/chr8.germline.match > germline/chr8.germline.match.clean
    input:
        germline=rules.merge_matches.output,
        estimated=rules.ersa_params.output
    output: "ersa/relatives.tsv"
    singularity:
        "docker://alexgenx/ersa:stable"
    shell:
        """
        #ERSA_L=13.7
        #ERSA_TH=3.197
        #ERSA_T=2.5
        source {input.estimated}
        ersa --avuncular-adj -t $ERSA_T -l $ERSA_L -th $ERSA_TH {input.germline} -o {output}
        """

rule merge_king_ersa:
    input:
        king=rules.run_king.output,
        germline=rules.merge_matches.output,
        ersa=rules.ersa.output
    output: "results/relatives.tsv"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/ersa_king.py"

rule evaluate_accuracy:
    input:
        rel=rules.merge_king_ersa.output[0],
        fam=rules.simulate.output['fam']
    output:
        "results/accuracy.png"
    conda: "../envs/evaluation.yaml"
    script: "../scripts/evaluate.py"