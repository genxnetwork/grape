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

### Build snakemake docker container

```shell script
docker build -t snakemake_t -f containers/snakemake/Dockerfile -m 4GB .
docker tag snakemake_t:latest alexgenx/snakemake:latest
docker push alexgenx/snakemake:latest
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

1. For lifting:
    chain = '/media/hg38ToHg19.over.chain.gz'
    ref = '/media/human_g1k_v37.fasta', size = 3G
2. For phasing:
    map = '/media/tables/genetic_map_hg19_withX.txt.gz', size = 51M
    bcf = '/media/1000genome/bcf/1000genome_chr{1..22}.bcf', size = 14G 
3. For imputation:
    m3vcf = '/media/Minimac/{1..22}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz', size = 3.2G 
4. For interpolation:
    map = '/media/genetic_map_b37/genetic_map_chr{1..22}_combined_b37.txt', size = 120M 
5. For imputation check:
    tab = '/media/1000genome/allele_info/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.only_rs.biallelic.tab', size = 12G
    