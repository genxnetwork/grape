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


### Pipeline checking 

To check if snakemake correctly sees input files do not pass --real-run to the launcher.py:

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data 
```

### How to run full pipeline

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
 --real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
 --real-run
```

## Evaluation on Simulated Data

### Requirements

1. Docker
2. Folder with all references with path /media/ref

### How to run simulation

Use command simulate. Options --input and --samples are not needed in this case.

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py simulate --directory /media/pipeline_data/simulation --real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --directory /media/pipeline_data/simulation --real-run 
```

## Evaluation on Hapmap Data

### Requirements

1. Docker
2. Folder with all references with path /media/ref

### How to run hapmap

Use command hapmap for preparing hapmap ceu data. Options --input and --samples are not needed in this case.
Then, use find

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py hapmap --directory /media/pipeline_data/hapmap --real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --directory /media/pipeline_data/hapmap --real-run
```