# run as "snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all"
# TODO: add reporting https://github.com/tanaes/snarkmark/blob/master/rules/report.rule

import pandas as pd
from os.path import join


configfile: "config.yaml"

REF_DIR         = config["reference"]["ref_dir"]
GRCh37_fasta    = join(REF_DIR, config["reference"]["GRCh37_fasta"])
GENETIC_MAP     = join(REF_DIR, config["reference"]["GENETIC_MAP"])
vcfRef          = join(REF_DIR, config["reference"]["vcfRef"])
refHaps         = join(REF_DIR, config["reference"]["refHaps"])
lift_chain      = join(REF_DIR, config["reference"]["lift_chain"])
cmmap           = join(REF_DIR, config["reference"]["cmmap"])
SITE_1000GENOME = join(REF_DIR, config["reference"]["SITE_1000GENOME"])
hapmap_ped      = join(REF_DIR, config["reference"]["hapmap_ped"])
hapmap_mp       = join(REF_DIR, config["reference"]["hapmap_mp"])
hapmap_fam      = join(REF_DIR, config["reference"]["hapmap_fam"])
hd_genotype_chip= join(REF_DIR, config["reference"]["hd_genotype_chip"])
pedsim_map      = join(REF_DIR, config["reference"]["pedsim_map"])

samples_file    = config["samples_file"]
SAMPLES         = pd.read_table(samples_file).set_index("name", drop=False)
SAMPLES_PATH    = SAMPLES.path.values.tolist()
SAMPLES_NAME    = SAMPLES.name.values.tolist()

CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']

print('VARABLES:')
print(REF_DIR, cmmap, samples_file)

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