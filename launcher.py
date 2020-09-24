import snakemake
import datetime
import psutil
import argparse
import shutil
import os


def get_parser_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--configfile",
        default="config.yaml",
        help="Snakemake YAML config file path")

    parser.add_argument(
        "--snakefile",
        default="Snakefile",
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

    parser.add_argument(
        '--simulate',
        action='store_true',
        help='If this argument is present, simulate data and run pipeline on it. You need to provide correct Snakefile from workflows/pedsim/Snakefile'
    )

    parser.add_argument(
        '--hapmap',
        action='store_true',
        help='If this argument is present, run pipeline on HapMap CEU data. You need to provide correct Snakefile from workflows/hapmap/Snakefile'
    )

    return parser.parse_args()




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
        os.mkdir(args.directory)

    if not (args.simulate or args.hapmap):
        copy_input(args.input, args.directory, args.samples)

    if 'CONDA_ENVS_PATH' not in os.environ:
        os.environ['CONDA_ENVS_PATH'] = '/tmp/envs'
    if 'CONDA_PKGS_DIRS' not in os.environ:
        os.environ['CONDA_PKGS_DIRS'] = '/tmp/conda/pkgs'

    print(os.environ)

    if not snakemake.snakemake(
            snakefile=args.snakefile,
            configfiles=[args.configfile],
            workdir=args.directory,
            cores=args.cores,
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