# GenX Relatives Detection Pipeline


## Atlas

### Requirements

1. Docker
2. samples.tsv tab-separated file with format and one line for each sample:
```text
name	path
1	input/1.txt
2	input/2.txt
```
3. Folder with inputs in 23andme format. One sample per file.
4. Folder with all references

### How to run full pipeline

```shell script

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \ 
launcher.py --samples /media/ref/samples.tsv --input /media/ref/input --directory /tmp/pipeline-real-run-1 \
--singularity-prefix /tmp --singularity-args -B /tmp:/tmp -W /tmp --conda-prefix /tmp --real-run
```

## Evaluation on Simulated Data

### Requirements

1. Docker
2. Folder with all references

### How to run simulation

Use option --simulate and pass different workflow description in workflows/pedsim/Snakefile. 
Options --input and --samples are not needed in this case.

```shell script
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \ 
launcher.py --directory /tmp/pipeline-dry-run-1 --singularity-prefix /tmp --singularity-args -B /tmp:/tmp -W /tmp --conda-prefix /tmp \
--real-run --simulate -s workflows/pedsim/Snakefile
```

## Evaluation on Hapmap Data

### Requirements

1. Docker
2. Folder with all references

### How to run hapmap

TODO: not implemented yet
