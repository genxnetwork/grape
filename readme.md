## Degree of kinship estimation pipeline
### Description

The project intends to implement best-practices of estimation of recent shared ancestry in a production-ready way.

### Stack

The pipeline is implemented with the Snakemake framework. All used components are wrapped in Singularity containers or isolated in a Conda environment.

The visualisation of execution graph: [svg](https://bitbucket.org/genxglobal/genx-relatives-snakemake/raw/077f33cfdd421ae17b5c02a3a5f8eb34bd20e1fd/dag.svg).


Multi-core parallelization is highly utilized due to the ability to split input data by each sample/chromosome.

Information about stages:

1. Preprocessing: we remove all multiallelic variants and indels.

2. Liftover: we use picard tools and lift data from hg38 to hg37.
3. Phasing: Eagle 2.4.1 and 1000 Genomes reference panel.
4. Imputation: Minimac4 and 1000 Genomes reference panel.
5. Close Relatives: KING IBD search.
6. IBD Search: Germline with merging closely located IBD segments together.
7. Distant Relatives: ERSA with default params estimated on CEU founders.
8. Merge: KING degree has priority over ersa degree for close relatives (degrees 1-3), otherwise, we take ERSA output.


### Installation

Clone the repository.

Download reference datasets to `/media/ref` folder using the following script:

```text

sftp genx-reference-sftp@20.54.91.13
Password: b2mR4wQpJJdeKdsW

```

Compile Funnel from https://github.com/messwith/funnel with go 1.12+ and make. Then one can just use bin/funnel binary.  
This Funnel fork simply adds the ‘--privileged’ flag to all task docker commands.  
Without ‘--privileged’ singularity containers do not work inside docker.

### Usage
#### Input data format

The information about samples for analysis should be provided as a path to a tab-separated text file (samples.tsv).

Input data is expected in 23andMe format, one file for each sample:

| name | path |
| --- | --- |
| 1 | input/1.txt |
| 2 | input/2.txt | 

#### Console execution

##### Pipeline checking

To check if snakemake correctly sees input files do not pass --real-run to the launcher.py:

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data
```

##### How to run full pipeline

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--real-run
```

#### Execution by scheduler
The pipeline can be executed using lightweight scheduler [Funnel](https://ohsu-comp-bio.github.io/funnel/), which implements [Task Execution Schema](https://github.com/ga4gh/task-execution-schemas) developed by [GA4GH](https://github.com/ga4gh/wiki/wiki).  
  
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
