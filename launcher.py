import snakemake
import datetime
import psutil
import argparse
import shutil
import os
from random import randint


from inspect import getsourcefile
from weight.ibd_segments_weigher import IBDSegmentsWeigher



# Returns an integer value for total available memory, in GB.
def total_memory_gb():
    n_bytes = psutil.virtual_memory().total
    return int(n_bytes / (1024 * 1024 * 1024))


def get_parser_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'command',
        default='find',
        help="""What pipeline should do, possible values are find, preprocess, simulate.
                preprocess converts hg38 per-sample 23andme input files to the one single vcf in hg37;
                find detects relatives in vcf file;
                simulate generates pedigree and vcf file with distant relatives from 1000 genomes CEU(CEPH) population using pedsim;
                reference downloads and preprocess all the references to the --ref-directory;
                For running the main pipeline you can provide .vcf file and use find or use preprocess with 23andme inputs"""
    )

    parser.add_argument(
        "--configfile",
        default="config.yaml",
        help="Snakemake YAML config file path")

    parser.add_argument(
        "--snakefile",
        default="",
        help="Snakemake file path")

    parser.add_argument(
        "--cores",
        default=max(1, psutil.cpu_count() - 1),
        type=int,
        help="Number of CPU cores to use in this pipeline run (default %(default)s)")

    parser.add_argument(
        "--num-batches",
        default=1,
        type=int,
        help="Number of batches to split input vcf.gz file into (default is 1). Splitting can speed up the preprocessing stage for the large datasets.")

    parser.add_argument(
        "--memory",
        default=max(1, total_memory_gb() - 1),
        type=int,
        help="Total memory (in GB) allowed for use by the Snakemake scheduler (default %(default)s)")

    parser.add_argument(
        "--client",
        help="The client mode assumes that background data is already been processed and placed in the background directory",
        action="store_true")

    parser.add_argument(
        "--real-run",
        help="If this argument is present, Snakemake will run the pipeline instead of dry-run, it is False by default",
        action="store_true")

    parser.add_argument(
        "--unlock",
        help="If this argument is present, Snakemake will simply unlock working directory without launch anything, it is False by default",
        action="store_true")

    parser.add_argument(
        "--directory",
        default=".",
        help="Snakemake working directory")

    parser.add_argument(
        "--ref-directory",
        default="",
        help="Snakemake directory with references. If empty, than it is equal to the /media/ref from config.yaml")

    parser.add_argument(
        '--input',
        default='input',
        help='Directory where input files are located'
    )

    parser.add_argument(
        '--samples',
        default='samples.tsv',
        help='File with the list of all samples'
    )

    parser.add_argument(
        '--assembly',
        default='hg38',
        help='Name of genome assembly. Default is hg38, the only other possible value is hg37. Hg38 data will be lifted to hg37'
    )

    parser.add_argument(
        '--remove-imputation',
        action='store_true',
        help='If present, preprocess workflow will remove all lines containing "IMPUTED" from input vcf file'
    )

    parser.add_argument(
        '--impute',
        action='store_true',
        help='If present, preprocess workflow will impute input vcf file. You MUST also add --phase to cmd in this case'
    )

    parser.add_argument(
        '--phase',
        action='store_true',
        help='If present, preprocess workflow will phase input vcf file'
    )

    parser.add_argument(
        '--vcf-file',
        default='input.vcf',
        help='Path to the input vcf file for "vcf" command only'
    )

    parser.add_argument(
        '--rule',
        default=None,
        help='Rule which will be rerun forcefully'
    )

    parser.add_argument(
        '--until',
        default=None,
        help='Rule which will be the last to run'
    )

    parser.add_argument(
        '--target',
        default=['all'],
        nargs='*',
        help='Target rules, snakemake will run only rules that lead to the input of this rules'
    )

    parser.add_argument(
        '--flow',
        default='ibis',
        help='How to find ibd segments: values are ibis, ibis-king, germline-king, default is ibis'
    )

    parser.add_argument(
        '--zero-seg-count',
        default=0.5,
        type=float,
        help="""
            Average count of IBD segments in two unrelated individuals in population.
            Smaller values of 0.1, 0.2 tend to give more distant matches than default 0.5.
            """
    )

    parser.add_argument(
        '--zero-seg-len',
        default=5.0,
        type=float,
        help="""
            Average length of IBD segment in two unrelated individuals in population.
            Smaller values of tend to give more distant matches than default 5.0
            """
    )

    parser.add_argument(
        '--alpha',
        default=0.01,
        type=float,
        help="""
            ERSA P-value limit for testing for an existence of an relationship.
            Values of 0.02-0.05 tend to give more distant matches that default 0.01.
            """
    )

    parser.add_argument(
        '--ibis-seg-len',
        default=7.0,
        type=float,
        help="""
                Minimum length of IBD segment for ibis.
                Smaller values of it tend to give more distant matches than default 7.0 and more false-positives.
            """
    )

    parser.add_argument(
        '--ibis-min-snp',
        default=500,
        type=int,
        help="""
                Minimum number of SNPs in IBD segment.
                Smaller values of it tend to give more distant matches than default 500 and more false-positives.
            """
    )
    
    parser.add_argument(
        '--rapid-error-rate',
        default=0.0025,
        type=float,
        help="""
                Estimation of per-SNP genotyping error rate.
                Smaller values lead to the smaller genotype window size for the RaPID. 
            """
    )
        
    parser.add_argument(
        '--rapid-min-snp',
        default=500,
        type=int,
        help="""
                Minimum number of SNPs in IBD segment.
                Smaller values of it tend to give more distant matches than default 500 and more false-positives.
            """
    )
            
    parser.add_argument(
        '--rapid-num-runs',
        default=10,
        type=int,
        help="""
                Number of RaPID random projections.
                More projections should give better accuracy. 
            """
    )
                
    parser.add_argument(
        '--rapid-num-success',
        default=2,
        type=int,
        help="""
                Number of segment matches in the --rapid-num-runs random projections.
                Larger values should give less false-positives.
            """
    )
    
    parser.add_argument(
        '--rapid-seg-len',
        default=5.0,
        type=float,
        help="""
                Minimum length of found IBD segments.
                Smaller values of it tend to give more distant matches than default 5.0 and more false-positives.
            """
    )
    
    parser.add_argument(
        '--stat-file',
        default='stat_file.txt',
        help='File for writing statistics'
    )

    parser.add_argument(
        "--sim-params-file",
        default="params/Relatives.def",
        help="Snakemake YAML config file path")

    parser.add_argument(
        "--sim-samples-file",
        default="ceph_unrelated_all.tsv",
        help="List of samples from 1000genomes for pedsim to use as founders. You can choose only from 'ceph_unrelated_all.tsv', 'all.tsv'")

    parser.add_argument(
        '--conda-prefix',
        default='/tmp',
        help='Conda prefix for environments'
    )

    parser.add_argument(
        '--use-bundle',
        action='store_true',
        default=False,
        help='Download all references as single file'
    )

    parser.add_argument(
        '--chip',
        default='background.vcf.gz',
        help='Path to chip file'
    )

    parser.add_argument(
        '--weight-mask',
        help='Mask of weights used to re-weight IBD segments length while using ERSA algorithm'
            'for `ibis` and `ibis-king` flows'
    )

    parser.add_argument(
        '--missing-samples',
        default=15.0,
        type=float,
        help='Upper bound of missing SNPs (%). Samples with higher values are removed from the relatedness detection analysis.')

    parser.add_argument(
        '--alt-hom-samples',
        default=1.0,
        type=float,
        help='Lower bound of homozygous alternative SNPs (%). Samples with lower values are removed from the relatedness detection analysis.')

    parser.add_argument(
        '--het-samples',
        default=5.0,
        type=float,
        help='Lower bound of heterozygous SNPs (%). Samples with lower values are removed from the relatedness detection analysis.')

    parser.add_argument(
        '--iqr-alpha',
        default=float("inf"),
        type=float,
        help='IQR multiplier for outliers detection by amount of missing samples. Samples below 1st quantile - IQR*alpha and above 3rd quantile + IQR*alpha are removed from the relatedness detection analysis.')

    parser.add_argument(
        '--seed',
        default=randint(0, 10**7),
        type=int,
        help='Random seed for Ped-sim pedigree simulation. The default value is randomly generated.')

    args = parser.parse_args()

    valid_commands = [
        'preprocess',
        'find',
        'simulate',
        'reference',
        'bundle',
        'compute-weight-mask',
        'simbig',
        'remove_relatives'
    ]

    if args.command not in valid_commands:
        raise RuntimeError(f'command {args.command} not in list of valid commands: {valid_commands}')

    if args.impute and not args.phase:
        raise ValueError('If --impute is present, then --phase must also be present')

    if args.command != 'reference' and args.use_bundle:
        raise ValueError('--use-bundle option only available for reference downloading')

    if args.num_batches > args.cores:
        raise ValueError('Number of batches is bigger than number cores, please change --num-batches value to be lower or equal --cores')

    if any((i < 0 or i > 100) for i in (args.het_samples, args.missing_samples, args.alt_hom_samples)):
        raise ValueError('Percentage cannot be higher than 100 or lower than 0')

    return args


def copy_file(working_dir, file_path):
    samples_name = os.path.split(file_path)[-1]
    samples_path = os.path.join(working_dir, samples_name)
    if not os.path.exists(samples_path):
        shutil.copy(file_path, os.path.join(working_dir, samples_name))


def copy_input(input_dir, working_dir, samples_file):
    input_name = os.path.split(input_dir)[-1]
    dest_path = os.path.join(working_dir, input_name)
    if not os.path.exists(dest_path):
        shutil.copytree(input_dir, dest_path)

    samples_name = os.path.split(samples_file)[-1]
    samples_path = os.path.join(working_dir, samples_name)
    if not os.path.exists(samples_path):
        shutil.copy(samples_file, os.path.join(working_dir, samples_name))


if __name__ == '__main__':

    args = get_parser_args()

    print(args)

    print()

    start_time = datetime.datetime.now()

    # in case when launcher.py is executed outside the Snakemake dir
    current_path = os.path.dirname(getsourcefile(lambda: 0))

    if not os.path.exists(args.directory):
        os.makedirs(args.directory)

    if args.command == 'simulate':
        copy_input(
            os.path.join(current_path, 'workflows/pedsim/params'),
            args.directory, os.path.join(current_path, 'workflows/pedsim/', args.sim_samples_file)
        )
        # for some reason launching with docker from command line
        # sets root directory for 'configfile' directive in bundle.Snakefile as snakemake.workdir
        # therefore config.yaml must be in snakemake.workdir
        shutil.copy(
            os.path.join(current_path, 'workflows/pedsim/config.yaml'),
            os.path.join(args.directory, 'config.yaml')
        )

    if args.command == 'simbig':
        copy_input(
            os.path.join(current_path, 'workflows/simbig/params'),
            args.directory, os.path.join(current_path, 'workflows/simbig/', args.sim_samples_file)
        )
        # for some reason launching with docker from command line
        # sets root directory for 'configfile' directive in bundle.Snakefile as snakemake.workdir
        # therefore config.yaml must be in snakemake.workdir
        shutil.copy(
            os.path.join(current_path, 'workflows/simbig/config.yaml'),
            os.path.join(args.directory, 'config.yaml')
        )

    if args.command == 'preprocess':
        shutil.copy(args.vcf_file, os.path.join(args.directory, 'input.vcf.gz'))

    if args.command in ['preprocess', 'find', 'reference', 'bundle', 'remove_relatives', 'compute-weight-mask']:
        if args.directory != '.':
            shutil.copy(os.path.join(current_path, 'config.yaml'), os.path.join(args.directory, 'config.yaml'))

    snakefiles = {
        'preprocess': 'workflows/preprocess2/Snakefile',
        'find': 'workflows/find/Snakefile',
        'simulate': 'workflows/pedsim/Snakefile',
        'reference': 'workflows/reference/Snakefile',
        'bundle': 'workflows/bundle/Snakefile',
        'simbig': 'workflows/simbig/Snakefile',
        'remove_relatives': 'workflows/remove_relatives/Snakefile',
        'compute-weight-mask': 'workflows/weight/Snakefile'
    }

    if args.client:
        background_path = os.path.join(args.directory, 'background/merged_imputed.vcf.gz')
        if not os.path.exists(background_path):
            raise RuntimeError(f'Background data is missing for the client mode')

    if not args.snakefile:
        if args.command == 'reference' and args.use_bundle:
            snakefile = os.path.join(current_path, snakefiles['bundle'])
        else:
            snakefile = os.path.join(current_path, snakefiles[args.command])
    else:
        snakefile = args.snakefile

    if 'CONDA_ENVS_PATH' not in os.environ:
        os.environ['CONDA_ENVS_PATH'] = '/tmp/envs'
    if 'CONDA_PKGS_DIRS' not in os.environ:
        os.environ['CONDA_PKGS_DIRS'] = '/tmp/conda/pkgs'

    print(os.environ)

    config_dict = {'mode': 'client'} if args.client is not None else {}
    config_dict['sim_params_file'] = args.sim_params_file
    config_dict['sim_samples_file'] = args.sim_samples_file
    config_dict['assembly'] = args.assembly
    config_dict['mem_gb'] = args.memory
    config_dict['num_batches'] = args.num_batches
    if args.ref_directory != '':
        config_dict['ref_dir'] = args.ref_directory
    if args.chip:
        config_dict['chip'] = args.chip
    if args.flow not in ['ibis', 'ibis-king', 'germline-king', 'rapid']:
        raise ValueError(f'--flow can be one of the ["ibis", "ibis-king", "germline-king", "rapid"] and not {args.flow}')
    config_dict['flow'] = args.flow
    if args.command in ['preprocess', 'simulate', 'reference', 'bundle', 'simbig', 'remove_relatives']:
        config_dict['remove_imputation'] = args.remove_imputation
        config_dict['impute'] = args.impute
        config_dict['phase'] = args.phase

    config_dict['zero_seg_len'] = args.zero_seg_len
    config_dict['zero_seg_count'] = args.zero_seg_count
    config_dict['alpha'] = args.alpha
    config_dict['ibis_seg_len'] = args.ibis_seg_len
    config_dict['ibis_min_snp'] = args.ibis_min_snp

    config_dict['rapid_error_rate'] = args.rapid_error_rate
    config_dict['rapid_min_snp'] = args.rapid_min_snp
    config_dict['rapid_num_runs'] = args.rapid_num_runs
    config_dict['rapid_num_success'] = args.rapid_num_success
    config_dict['rapid_seg_len'] = args.rapid_seg_len
    config_dict['missing_samples'] = args.missing_samples
    config_dict['alt_hom_samples'] = args.alt_hom_samples
    config_dict['het_samples'] = args.het_samples
    config_dict['iqr_alpha'] = args.iqr_alpha

    config_dict['seed'] = args.seed

    if args.weight_mask:
        config_dict['weight_mask'] = os.path.join(args.directory, args.weight_mask)
        config_dict['ersa_r'] = IBDSegmentsWeigher.from_json_mask_file(config_dict['weight_mask']) \
            .adjusted_expected_recombination_number

    if not snakemake.snakemake(
            snakefile=snakefile,
            configfiles=[args.configfile or 'config.yaml'],
            config=config_dict,
            workdir=args.directory if args.command != 'reference' else args.ref_directory,
            cores=args.cores,
            unlock=args.unlock,
            printshellcmds=True,
            dryrun=(not args.real_run),
            targets=args.target,
            stats=args.stat_file,
            forcerun=[args.rule] if args.rule is not None else [],
            until=[args.until] if args.until is not None else [],
            use_conda=True,
            conda_prefix=args.conda_prefix,
            envvars=['CONDA_ENVS_PATH', 'CONDA_PKGS_DIRS'],
            keepgoing=True
    ):
        raise ValueError("Pipeline failed see Snakemake error message for details")

    print(args)

    end_time = datetime.datetime.now()
    print("--- Pipeline running time: %s ---" % (str(end_time - start_time)))
