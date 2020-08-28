# GenX relatives detection pipeline

### Snakemake launch

```shell script
snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all
```