# GenX Relatives Detection Pipeline


## Real data in 23andme format and Grch38 assembly

### Requirements

1. Docker
2. samples.tsv tab-separated file with format and one line for each sample:

```text
name	path
1	input/1.txt
2	input/2.txt
```
3. Folder with inputs in 23andme format. One sample per file.
4. Folder with all references with path /media/ref

### How to run full pipeline

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \ 
-e CONDA_ENVS_PATH=/tmp/envs -e CONDA_PKGS_DIRS=/tmp/conda/pkgs genx_relatives:latest \
launcher.py --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--singularity-prefix /tmp --singularity-args='-B /media:/media -B /tmp:/tmp -W /tmp' --conda-prefix /tmp --real-run
```

## Evaluation on Simulated Data

### Requirements

1. Docker
2. Folder with all references with path /media/ref

### How to run simulation

Use option --simulate and pass different workflow description in workflows/pedsim/Snakefile. 
Options --input and --samples are not needed in this case.

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \ 
-e CONDA_ENVS_PATH=/tmp/envs -e CONDA_PKGS_DIRS=/tmp/conda/pkgs genx_relatives:latest \
launcher.py --directory /media/pipeline_data/simulation \
--singularity-prefix /tmp --singularity-args='-B /media:/media -B /tmp:/tmp -W /tmp' --conda-prefix /tmp \
--real-run --simulate -s workflows/pedsim/Snakefile
```

## Evaluation on Hapmap Data

### Requirements

1. Docker
2. Folder with all references with path /media/ref

### How to run hapmap

Use option --hapmap and pass different workflow description in workflows/hapmap/Snakefile. 
Options --input and --samples are not needed in this case.

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \ 
-e CONDA_ENVS_PATH=/tmp/envs -e CONDA_PKGS_DIRS=/tmp/conda/pkgs genx_relatives:latest \
launcher.py --directory /media/pipeline_data/hapmap \
--singularity-prefix /tmp --singularity-args='-B /media:/media -B /tmp:/tmp -W /tmp' --conda-prefix /tmp \
--real-run --hapmap -s workflows/hapmap/Snakefile

```