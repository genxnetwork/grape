# GenX relatives detection pipeline

### Snakemake launch

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

