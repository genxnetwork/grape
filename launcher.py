import snakemake
import datetime
import psutil
import argparse
import shutil
import os


def get_parser_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('command',
                        default='find',
                        help="""What pipeline should do, possible values are find, preprocess, simulate, hapmap.
                        preprocess converts hg38 per-sample 23andme input files to the one single vcf in hg37;
                        find detects relatives in vcf file;
                        simulate generates pedigree and vcf file with distant relatives from 1000 genomes CEU(CEPH) population using pedsim;
                        hapmap extracts CEU(CEPH) data for running find;
                        For running the main pipeline you can provide .vcf file and use find or use preprocess with 23andme inputs""")

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
        '--stat-file',
        default='stat_file.txt',
        help='File for writing statistics'
    )

    parser.add_argument(
        "--sim-params-file",
        default="params/Relatives.def",
        help="Snakemake YAML config file path")

    # --singularity-prefix /tmp --singularity-args='-B /media:/media -B /tmp:/tmp -W /tmp' --conda-prefix /tmp
    parser.add_argument(
        '--singularity-prefix',
        default='/tmp',
        help='Directory where snakemake will put singularity images'
    )
    parser.add_argument(
        '--singularity-args',
        default='-B /media:/media -B /tmp:/tmp -W /tmp',
        help='Additional singularity arguments'
    )
    parser.add_argument(
        '--conda-prefix',
        default='/tmp',
        help='Conda prefix for environments'
    )

    return parser.parse_args()


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

    if not os.path.exists(args.directory):
        os.makedirs(args.directory)

    valid_commands = ['preprocess', 'find', 'simulate', 'hapmap']
    if args.command not in valid_commands:
        raise RuntimeError(f'command {args.command} not in list of valid commands: {valid_commands}')

    if args.command == 'preprocess':
        copy_input(args.input, args.directory, args.samples)

    if args.command == 'simulate':
        copy_input('workflows/pedsim/params', args.directory, 'workflows/pedsim/ceph_unrelated_all.tsv')
        # for some reason launching with docker from command line
        # sets root directory for 'configfile' directive in Snakefile as snakemake.workdir
        # therefore config.yaml must be in snakemake.workdir
        shutil.copy('workflows/pedsim/config.yaml', os.path.join(args.directory, 'config.yaml'))

    if args.command in ['preprocess', 'find']:
        shutil.copy('config.yaml', os.path.join(args.directory, 'config.yaml'))

    snakefiles = {
        'preprocess': 'workflows/preprocess/Snakefile',
        'find': 'Snakefile',
        'simulate': 'workflows/pedsim/Snakefile',
        'hapmap': 'workflows/hapmap/Snakefile'
    }

    if not args.snakefile:
        snakefile = snakefiles[args.command]
    else:
        snakefile = args.snakefile

    if 'CONDA_ENVS_PATH' not in os.environ:
        os.environ['CONDA_ENVS_PATH'] = '/tmp/envs'
    if 'CONDA_PKGS_DIRS' not in os.environ:
        os.environ['CONDA_PKGS_DIRS'] = '/tmp/conda/pkgs'

    print(os.environ)
    print()
    print(os.getcwd())
    print()
    print(os.listdir('.'))
    if not snakemake.snakemake(
            snakefile=snakefile,
            #configfiles=[args.configfile],
            config={'sim_params_file': args.sim_params_file},
            workdir=args.directory,
            cores=args.cores,
            unlock=args.unlock,
            printshellcmds=True,
            dryrun=(not args.real_run),
            targets=['all'],
            stats=args.stat_file,
            forcerun=[args.rule] if args.rule is not None else [],
            until=[args.until] if args.until is not None else [],
            use_conda=True,
            conda_prefix=args.conda_prefix,
            use_singularity=True,
            singularity_prefix=args.singularity_prefix,
            singularity_args=args.singularity_args,
            envvars=['CONDA_ENVS_PATH', 'CONDA_PKGS_DIRS']
    ):
        raise ValueError("Pipeline failed see Snakemake error message for details")

    print(args)

    end_time = datetime.datetime.now()
    print("--- Pipeline running time: %s ---" % (str(end_time - start_time)))