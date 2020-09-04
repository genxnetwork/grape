# GenX relatives detection pipeline

### Snakemake launch

    snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all

## Visualization of the DAG

    snakemake --dag all | dot -Tsvg > dag.svg

```shell script
snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all
```

### Force-launch single rule

```shell script
snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -R somerule --until somerule
```

### Simulate data

```shell script
snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p -s workflows/pedsim/Snakefile
```

### Launch using docker container

```shell script

cp -r input /media/pipeline_data/atlas40/
cp samples.tsv /media/pipeline_data/atlas40/

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro alexgenx/snakemake:latest /bin/bash
# this is inside docker container

snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media/singulariry_cache --singularity-args="-B /media:/media" -p --configfile config.yaml --directory /media/pipeline_data/atlas40 -n

```

### References

1. 