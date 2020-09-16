# run as "snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all"
# TODO: add reporting https://github.com/tanaes/snarkmark/blob/master/rules/report.rule

import pandas as pd

#configfile: "config.yaml"

samples_file    = config["samples_file"]
SAMPLES         = pd.read_table(samples_file).set_index("name", drop=False)
SAMPLES_PATH    = SAMPLES.path.values.tolist()
SAMPLES_NAME    = SAMPLES.name.values.tolist()

CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']

def get_samples_name(wildcards):
    return SAMPLES.loc[int(wildcards.sample), "name"] # mind the int index

def get_samples_path(wildcards):
    return SAMPLES.loc[int(wildcards.sample), "path"] # mind the int index

# include: "rules/report.rule"

rule all:
    input:
        "results/relatives.tsv",
        # rules.report_benchmark_summary.output

include: "rules/preprocessing.smk"
include: "rules/filter.smk"
include: "rules/imputation.smk"
include: "rules/relatives.smk"