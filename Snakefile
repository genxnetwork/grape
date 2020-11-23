# run as "snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all"
# TODO: add reporting https://github.com/tanaes/snarkmark/blob/master/rules/report.rule

import pandas as pd
import os
from os.path import join


configfile: "config.yaml"

use_rapid       = config["use_rapid"]
is_client       = config["mode"] == "client"
use_simulated_ibd = config["use_simulated_ibd"] if "use_simulated_ibd" in config else False
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

print('execution mode is: ' + config["mode"])
print('dirs:')
print(os.getcwd())
print(os.listdir('.'))

CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam', 'log']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'log', 'nosex']


# include: "rules/report.rule"

# if config["mode"] == "all":
#     ruleorder: convert_imputed_to_plink > merge_convert_imputed_to_plink
# else:
#     ruleorder: merge_convert_imputed_to_plink > convert_imputed_to_plink

rule all:
    input:
        "results/relatives.tsv",
        # rules.report_benchmark_summary.output

#include: "rules/preprocessing.smk"
include: "rules/filter.smk"
include: "rules/imputation.smk"
if use_rapid:
    include: "rules/relatives_rapid.smk"
else:
    include: "rules/relatives.smk"