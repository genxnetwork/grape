# run as "snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all"
# TODO: add reporting https://github.com/tanaes/snarkmark/blob/master/rules/report.rule

from os.path import join


configfile: "config.yaml"

flow               = config["flow"]
is_client          = config["mode"] == "client"
use_simulated_ibd  = config["use_simulated_ibd"] if "use_simulated_ibd" in config else False

REF_DIR            = config["ref_dir"]
GRCH37_FASTA       = join(REF_DIR, config["reference"]["GRCh37_fasta"])
GENETIC_MAP        = join(REF_DIR, config["reference"]["GENETIC_MAP"])
GENETIC_MAP_GRCH37 = join(REF_DIR, config["reference"]["genetic_map_GRCh37"])
REF_VCF            = join(REF_DIR, config["reference"]["vcfRef"])
REF_HAPS           = join(REF_DIR, config["reference"]["refHaps"])
LIFT_CHAIN         = join(REF_DIR, config["reference"]["lift_chain"])
CMMAP              = join(REF_DIR, config["reference"]["cmmap"])
SITE_1000GENOME    = join(REF_DIR, config["reference"]["SITE_1000GENOME"])
HAPMAP_PED         = join(REF_DIR, config["reference"]["hapmap_ped"])
HAPMAP_MP          = join(REF_DIR, config["reference"]["hapmap_mp"])
HAPMAP_FAM         = join(REF_DIR, config["reference"]["hapmap_fam"])
HD_GENOTYPE_CHIP   = join(REF_DIR, config["reference"]["hd_genotype_chip"])
PEDSIM_MAP         = join(REF_DIR, config["reference"]["pedsim_map"])

print('execution mode is: ' + config["mode"])

CHROMOSOMES     = [str(i) for i in range(1, 23)]
PLINK_FORMATS   = ['bed', 'bim', 'fam']
PLINK_FORMATS_EXT   = ['bed', 'bim', 'fam', 'nosex']

# include: "rules/report.rule"

# if config["mode"] == "all":
#     ruleorder: convert_imputed_to_plink > merge_convert_imputed_to_plink
# else:
#     ruleorder: merge_convert_imputed_to_plink > convert_imputed_to_plink

rule all:
    input:
        "results/relatives.tsv",
        # rules.report_benchmark_summary.output


if flow == 'ibis':
    include: "../../rules/relatives_ibis.smk"
elif flow == 'germline':
    include: "../../rules/relatives.smk"
elif flow == 'ibis_king':
    include: "../../rules/relatives_ibis_king.smk"
