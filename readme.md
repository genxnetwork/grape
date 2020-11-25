[TOC]

## Degree of kinship estimation pipeline

### Description

This project is intended to implement best-practices of estimation of recent shared ancestry in a production-ready way.

#### Overview

The pipeline is implemented with the Snakemake workflow management system. All internal tools needed for execution are wrapped in Singularity containers or isolated in a Conda environment. We recommend using it with the provided Docker image, which already possesses all the needed dependencies.

You should have the following datasets in place in order to run the pipeline:

* input samples in one of the formats:
  * 23andme (multiple) along with the samples desctiption (see bellow for details)
  * VCF (single)
* reference genome and associated files (genetic map, lift chain, sites, etc)

#### Information about stages

The main worklfow steps

1. Preprocessing: we remove all multiallelic variants and indels.
2. Liftover: we use picard tools and lift data from hg38 to hg37.
3. Phasing: Eagle 2.4.1 and 1000 Genomes reference panel.
4. Imputation: Minimac4 and 1000 Genomes reference panel.
5. Close Relatives: KING IBD search.
6. IBD Search: Germline with merging closely located IBD segments together.
7. Distant Relatives: ERSA with default params estimated on CEU founders.
8. Merge: KING degree has priority over ersa degree for close relatives (degrees 1-3), otherwise, we take ERSA output.

The visualisation of execution graph: [svg](https://bitbucket.org/genxglobal/genx-relatives-snakemake/raw/077f33cfdd421ae17b5c02a3a5f8eb34bd20e1fd/dag.svg).

Multi-core parallelization is highly utilized due to the ability to split input data by each sample/chromosome.

### Installation

1. Clone the repository:

```text

git clone git@bitbucket.org:genxglobal/genx-relatives-snakemake.git

```

2. Download reference datasets (in our examples is `/media/ref`) folder using the following credentials:

```text

sftp genx-reference-sftp@20.54.91.13
Password: b2mR4wQpJJdeKdsW

```

3. (Optional) Compile Funnel from https://github.com/messwith/funnel with go 1.12+ and make. Then one can just use bin/funnel binary.  
In this Funnel fork we simply added the ‘--privileged’ flag to all task docker commands because w/o this flag  singularity containers do not work inside docker.

### Usage

At the moment, pipeline supports 23andme (separate files) and VCF (single) as an input.

#### Input data format: 23andme

The information about samples for analysis should be provided as a path to a tab-separated text file (samples.tsv).

Input data is expected in 23andMe format, one file for each sample:

| name | path |
| --- | --- |
| 1 | input/1.txt |
| 2 | input/2.txt | 

#### Input data format: vcf

Another option is using gzipped vcf file format. If vcf file is in hg38 assembly, then you can just use `vcf` command 
with the `--vcf-file <path>` option. If vcf file is in hg37, you should pass `--assembly hg37` to the `vcf` command.   

#### Console execution

##### Pipeline checking

First, it is suggested to run pipeline in dry-run mode (w/o --real-run flag) to check if Snakemake correctly sees input files:

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data

```

##### How to run full pipeline

Add --real-run flag after succesfull dry-run for the production run

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

# if input data is in 23andme format
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--real-run

# if input data is in vcf format the execution is in two steps
# first, prepare the input vcf file
# use --assembly hg37 if vcf file is in hg37 and not in hg38 
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py vcf --vcf-file /media/ref/input.vcf.gz --directory /media/pipeline_data/real-data \
--real-run

# second, now we can find relatives

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--real-run

```

#### Execution by scheduler

The pipeline can be executed using lightweight scheduler [Funnel](https://ohsu-comp-bio.github.io/funnel/), which implements [Task Execution Schema](https://github.com/ga4gh/task-execution-schemas) and developed by [GA4GH](https://github.com/ga4gh/wiki/wiki).  
  
During execution, incoming data for analysis can be obtained in several ways: locally, FTP, HTTPS, S3, Google, etc.  
The resulting files can be uploaded in the same ways. It is possible to add another feature such as writing to the database, sending to the REST service.  
The scheduler itself can work in various environments from a regular VM to a Kubernetes cluster with resource quotas support.  

More information: https://ohsu-comp-bio.github.io/funnel/docs/  

How to execute dry-run (sample output):

```text 

# Firstly, if server is not running, start it
/path/to/funnel server run

# Then, use funnel as client
/path/to/funnel task create examples/snakemake-dry-23andme.json

```

How to execute operational run (sample output):

```text 

/path/to/funnel task create examples/snakemake-real-23andme.json

```

#### Standalone version (not recommended)

It is possible to run the pipeline using standalone version. First, you need to clone the repository and setup the references as described above.

The main idea of using Docker containers that you no need configure your execution environment manually. You need to do it before you can run the pipeliune in the standalone mode.

Assuming that you use Ubuntu 18 the following steps are needed to run the pipeline:

1. Docker

https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04

It is recomended to move docker storage from /var to some other partition with enought free space:
https://www.guguweb.com/2019/02/07/how-to-move-docker-data-directory-to-another-location-on-ubuntu/

2. Singularity

You need to compile at least version 3.x from https://github.com/hpcng/singularity/releases/
Please note that you need Go compiler in order to do so https://golang.org/dl/

3. Conda

Snakemake pipeline can use Singularity containers (same as Docker but working from user space) as well as Conda wrapped tools for the virtualization of the execution steps.

https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-18-04

4. Snakemake

https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

5. Setup env vars for temp / cache directories for Singularity, ex:

export SINGULARITY_TMPDIR=/media/tmp
export SINGULARITY_CACHEDIR=/media/tmp

You can also pass this values using Snakemake:
--singularity-prefix DIR
--singularity-args ARGS

6. Launch

```text
snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all -n
```

Please mind '-n' flag for dry-run

### Usefull commands

Please see [useful_commands.md](useful_commands.md).

### Evaluation on Simulated Data

Pedigree simulation is performed on European populations from 1KG using the `pedsim` package.  
Pedsim can use sex-specific genetic maps and randomly assigns the sex of each parent (or uses user-specified sexes) when using such maps.  
Sex-specific map is preferrable because men and women have different recombination rates.  
Founders for the pedigree simulation are selected from 1000genomes HD genotype chip data, CEU population.  
CEU data consists of trios and we select no more than one member of each trio as founder.  

Visualization of structure of simulated pedigree is given below:

![pedigree](https://bitbucket.org/genxglobal/genx-relatives-snakemake/downloads/pedsim.png)

#### How to run simulation

Use command simulate. Options --input and --samples are not needed in this case.

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py simulate --directory /media/pipeline_data/simulation --real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --directory /media/pipeline_data/simulation --real-run
```

#### Evaluation results


![results](https://bitbucket.org/genxglobal/genx-relatives-snakemake/downloads/accuracy_merged2.png)

The pipeline shows good accuracy for degrees from 1 to 6.  
The results for degrees from 7 to 8 been improved by merging small IBD segments together if they are located near each other, thus,  
it detects more than a half of pairs with those degrees of kinship.

### Evaluation on Hapmap Data

Use command HapMap for preparing Hapmap CEU data. Options --input and --samples are not needed in this case.
Then, use find

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py hapmap --directory /media/pipeline_data/hapmap --real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --directory /media/pipeline_data/hapmap --real-run
```

#### Evaluation results

HapMap has information only about close relatives presented. The pipeline determines them with 100% accuracy.


### Credits

License: GNU GPL v3
