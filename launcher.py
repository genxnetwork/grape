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
        "--dry-run",
        help="If this argument is present, Snakemake will do a dry run of the pipeline, it is True by default",
        action="store_false")

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
        '--stat_file',
        default='stat_file.txt',
        help='File for writing statistics'
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

    start_time = datetime.datetime.now()
    copy_input(args.input, args.directory, args.samples)
    # Please mind dryrun = True!
    if not snakemake.snakemake(
            snakefile=args.snakefile,
            configfiles=[args.configfile],
            workdir=args.directory,
            cores=args.cores,
            printshellcmds=True,
            dryrun=args.dry_run,
            targets=['all'],
            stats=args.stat_file,
            forcerun=[args.rule] if args.rule is not None else [],
            until=[args.until] if args.until is not None else [],
            use_conda=True,
            use_singularity=True,
            singularity_prefix='/media/singulariry_cache',
            singularity_args='-B /media:/media'
    ):
        raise ValueError("Pipeline failed see Snakemake error message for details")

    print(args)

    end_time = datetime.datetime.now()
    print("--- Pipeline running time: %s ---" % (str(end_time - start_time)))