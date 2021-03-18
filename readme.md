## Degree of kinship estimation pipeline
### Description

The project intends to implement best-practices of estimation of recent shared ancestry in a production-ready way.

Main features:

1. It can handle input in hg37 and hg38.
2. Implements phasing and imputation pipeline with GERMLINE + ERSA recent shared ancestry estimation
3. It has a very fast alternative pipeline without phasing and imputation. It uses IBIS + ERSA.
4. It has a special simulation workflow for the accuracy analysis.
5. It is fully containerized in Docker.
6. Fast pipeline workflow with `--flow ibis` option can process 2000 samples in a few minutes.

### Stack

The pipeline is implemented with the Snakemake framework. All used components are wrapped in Singularity containers or isolated in a Conda environment.

The visualisation of execution graph: [svg](https://bitbucket.org/genxglobal/genx-relatives-snakemake/raw/077f33cfdd421ae17b5c02a3a5f8eb34bd20e1fd/dag.svg).

Multi-core parallelization is highly utilized due to the ability to split input data by each sample/chromosome.

Information about stages:

**Germline workflow**

1. Preprocessing: we remove all multiallelic variants and indels.
2. Liftover: we use picard tools and lift data from hg38 to hg37.
3. Phasing: Eagle 2.4.1 and 1000 Genomes reference panel.
4. Imputation: Minimac4 and 1000 Genomes reference panel.
5. Close Relatives: KING IBD search.
6. IBD Search: Germline with merging closely located IBD segments together.
7. Distant Relatives: ERSA with default params estimated on CEU founders.
8. Merge: KING degree has priority over ersa degree for close relatives (degrees 1-3), otherwise, we take ERSA output.

**Fast IBIS workflow**

1. Preprocessing: we remove all multiallelic variants and indels.
2. Liftover: we use picard tools and lift data from hg38 to hg37.
3. Close Relatives: KING IBD search.
4. IBD Search: IBIS.
5. Distant Relatives: ERSA with default params estimated on CEU founders.
6. Merge: KING degree has priority over ersa degree for close relatives (degrees 1-3), otherwise, we take ERSA output.


### Installation

1. Clone the repository.

2. Build docker container 
```   
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .
```

3. Download all needed references to the `--ref-directory` of your choice. 
```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py reference  --ref-directory /media/ref

```

4. (Optional) Compile Funnel from https://github.com/messwith/funnel with go 1.12+ and make. Then one can just use bin/funnel binary.  
This Funnel fork simply adds the ‘--privileged’ flag to all task docker commands.  
Without ‘--privileged’ singularity containers do not work inside docker.

### Usage

#### Reference downloading

Firstly, one needs to download all needed references to the `--ref-directory` of your choice.
These references will take up to 40GB of disk space. `--ref-directory` argument can now be used in all subsequent commands
```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py reference  --ref-directory /media/ref

```
#### Input data format: vcf

One option is using gzipped vcf file format. If vcf file is in hg38 assembly, then you can just use `vcf` command
with the `--vcf-file <path>` option. In this case, input file will be lifted. 
If vcf file is in hg37, you should pass `--assembly hg37` to the `vcf` command.

#### Input data format: 23andme

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
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--ref-directory /media/ref
```

##### How to run the full pipeline with phasing and imputation

With input in .vcf.gz format:

```text

docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

# use --assembly hg37 if vcf file is in hg37 and not in hg38 
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py vcf --vcf-file /media/ref/input.vcf.gz --directory /media/pipeline_data/real-data --ref-directory /media/ref \
--real-run

# now we can find relatives
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --directory /media/pipeline_data/real-data --ref-directory /media/ref \
--real-run
```

With input in 23andme format:

```text
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--ref-directory /media/ref --real-run

# now we can find relatives
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--ref-directory /media/ref --real-run
```

##### Very fast relatives detection using king and ibis

You should just add ```--flow ibis``` to the ```find``` command of ```launcher.py```.

```text
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --samples /media/ref/samples.tsv --input /media/ref/input --directory /media/pipeline_data/real-data \
--flow ibis --ref-directory /media/ref --real-run 
```

In this case, nothing will be phased or imputed. Slight loss of accuracy is possible for degrees 8-10, 
especially if you use different chips in the same batch


#### Description of output file

Output file is in .tsv file format, and contains one line for each detected pair of relatives.

```text
id1	id2	king_degree	king_relation	shared_genome_proportion	kinship	ersa_degree	final_degree	total_seg_len	seg_count
HGDP00274_HGDP00274	HGDP00315_HGDP00315	3	3	0.1216		4	3	764.7229547063728	77
HGDP00274_HGDP00274	HGDP00319_HGDP00319			0.031042068715083804		5	5	222.26121200000003	23.0
```

 * `id1` - ID of first sample in a pair of relatives.
 * `id2` - ID of second sample in a pair of relatives, `id1` is always less than `id2` by the rules of string comparison in python.
 * `king_degree` - Numeric degree of relationship estimated by KING. 0 means duplicates or MZ twins, 
   1 means parent-offspring (PO), 2 can be either full siblings (FS), half siblings and grandmother/grandfather with a granddaughter/grandson.
   3 is aunt/uncle with a niece/nephew, as described in table in https://en.wikipedia.org/wiki/Coefficient_of_relationship.
   If `king_degree` exists, then `final_degree` will be equal to `king_degree`. 
 * `king_relation` - further differentiation for first 2 degrees of KING. `king_degree` 1 means PO - parent-offspring, 
   also KING detects FS in some of second degrees.
 * `shared_genome_proportion` is the approximate fraction of genome shared by two individuals. 
   It should be approximately 0.5 for PO and FS, 0.25 for grandmother-granddaugher and half-siblings.
 * `kinship` is the KING kinship coefficient.
 * `ersa_degree` is the degree estimated from IBD segments by ERSA, it is used for the `final_degree` in the cases where `king_degree` does not exist.
 * `final_degree` is simply `king_degree` for close relatives up to 3rd degree and `ersa_degree` for distant relatives.
 * `total_seg_len` is the total length of all IBD segments, for the first 3 degrees it is calculated using KING IBD data, 
   for the 4th+ degrees it is calculated using IBID or Germline IBD data.
 * `seg_count` is the total number of all IBD segments found by KING for the first 3 degrees and found by IBIS\Germline for the 4th+ degrees.



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

7. Useful commands

Please see useful_commands.md.

### Known limitations

1. It is known that in some small, isolated populations IBD sharing is very high. 
   Therefore, our pipeline will overestimate the relationship degree for them. 
   It is not recommended to mix standard populations like CEU and small populations as isolated native ones. 
   This problem is discussed in https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034267 .
2. If one intends to analyze diverse datasets, it is recommended to impute them with the same pipeline. 
   It can be done with our pipeline using `--until` option:
    ```text
    
    docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .
    
    # use --assembly hg37 if vcf file is in hg37 and not in hg38 
    docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
    launcher.py vcf --vcf-file /media/ref/input.vcf.gz --directory /media/pipeline_data/imputed_data \
    --real-run
    
    # now we can find relatives
    docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
    launcher.py find --directory /media/pipeline_data/imputed_data \
    --until merge_imputation_filter --real-run
    ```
    Then, one can grab file `/media/pipeline_data/imputed_data/vcf/merged_imputed.vcf.gz`


### Evaluation on Simulated Data

Pedigree simulation is performed on European populations from 1KG using the `pedsim` package.  
Pedsim can use sex-specific genetic maps and randomly assigns the sex of each parent (or uses user-specified sexes) when using such maps.  
Sex-specific map is preferable because men and women have different recombination rates.  
Founders for the pedigree simulation are selected from 1000genomes HD genotype chip data, CEU population.  
CEU data consists of trios, and we select no more than one member of each trio as founder.  

Visualization of structure of simulated pedigree is given below:

![pedigree](https://bitbucket.org/genxglobal/genx-relatives-snakemake/downloads/pedsim.png)

#### How to run simulation

Use command simulate. Options --input and --samples are not needed in this case.

```text
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB .

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py simulate --directory /media/pipeline_data/simulation --ref-directory /media/ref --real-run
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
launcher.py hapmap --directory /media/pipeline_data/hapmap --ref-directory /media/ref --real-run

docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --directory /media/pipeline_data/hapmap --ref-directory /media/ref --real-run
```

#### Evaluation results

HapMap has information only about close relatives presented. The pipeline determines them with 100% accuracy.


### Credits

License: GNU GPL v3
