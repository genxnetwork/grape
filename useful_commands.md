# GenX relatives detection pipeline

### Launcher

```text

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro grape:latest \
launcher.py find --samples /media/samples.tsv --input /media/input --directory /media/data

```

### Snakemake launch

```text

snakemake --cores all --use-conda -p all

```

### Visualization of the DAG

```text


    snakemake --dag all | dot -Tsvg > dag.svg

```

```text

snakemake --cores all --use-conda -p all

```

### Force-launch single rule

```text

snakemake --cores all --use-conda -R somerule --until somerule

```

### Simulate data

```text

snakemake --cores all --use-conda -p -s workflows/pedsim/Snakefile

```

### Build snakemake docker container

```text

docker build -t alexgenx/snakemake:latest -f containers/snakemake/Dockerfile -m 4GB .
docker push alexgenx/snakemake:latest

```

### Launch using docker container

```text

snakemake --cores all --use-conda -p --configfile config.yaml --directory /media/data -n

```

### How do I trigger re-runs for rules with updated input files

```text

snakemake -n -R `snakemake --list-input-changes`

```

### Clean-up Snakemake working dir

If you want to use the same working directory for different input files it is better to clean-up it first

```text

snakemake --delete-all-outputs --cores 1

```

### References

1. For lifting:
    chain = '/media/ref/hg38ToHg19.over.chain.gz'
    ref = '/media/ref/human_g1k_v37.fasta', size = 3G
2. For phasing:
    map = '/media/ref/tables/genetic_map_hg19_withX.txt.gz', size = 51M
    bcf = '/media/ref/1000genome/bcf/1000genome_chr{1..22}.bcf', size = 14G 
3. For imputation:
    m3vcf = '/media/ref/Minimac/{1..22}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz', size = 3.2G 
4. For interpolation:
    map = '/media/ref/genetic_map_b37/genetic_map_chr{1..22}_combined_b37.txt', size = 120M 
5. For imputation check:
    tab = '/media/ref/1000genome/allele_info/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.only_rs.biallelic.tab', size = 12G
    