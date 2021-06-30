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

![image](https://user-images.githubusercontent.com/19895289/123947586-35c90900-d9a9-11eb-976a-b98b1a25ba5d.png)

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

Pipeline has three steps: reference downloading, preprocessing data and finding relatives. 
Reference downloading should only be done once. Before using it, one needs to build a docker container:

```
docker build -t genx_relatives:latest -f containers/snakemake/Dockerfile -m 8GB . 
```

#### Reference downloading

Firstly, one needs to download all needed references to the `--ref-directory` of your choice.
These references will take up to **40GB** of disk space. `--ref-directory` argument can now be used in all subsequent commands
```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py reference --real-run --ref-directory /media/ref --directory /media/ref
```

Please note that reference post-processing is a computationally intensive task and can take a lot of time. Consider using `--cores all` (or number) flag because the default behavior uses only 1 core.

#### Preprocessing

Preprocessing features:
    - Lifting to hg37 if input file is in hg38, invoked by option `--assembly hg38`. 
    - Optional removing of imputed SNPs. In our tests 600K SNPs is enough and further imputation does not improve quality.
     If your dataset has less then 300K SNPs, imputation can improve results.
     Search of relatives requires substantially more time with increased number of SNPs.  
     Invoked by option `--remove-imputation`. Currently it removes all SNPs with `IMPUTED` in it.
    - Optional phasing, invoked by `--phase`. It is required for searching for relatives with Germline `--flow germline`, 
    however we do not recommend using Germline for now.
    - Optional imputation, invoked by `--impute`. One must also use `--phase` with this.
    
Option `--ref-directory` required to point to the downloaded references, option `--directory` points to the working directory
where results and some files from intermediate steps will be saved.
   
**Important: add `--real-run` to all of these commands if you want a real launch and not just building of computational graph**
 
Simple command when input file is in hg37 already and removing imputation, phasing and imputation is not required:
 
```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --ref-directory /media/ref --vcf-file /media/input.vcf.gz --directory /media/pipeline_data/real-data 
```

Command with lifting from hg38 to hg37:

```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --ref-directory /media/ref --vcf-file /media/input.vcf.gz --directory /media/pipeline_data/real-data \
--assembly hg38  
```

Command with lifting from hg38 to hg37, phasing and imputation:

```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py preprocess --ref-directory /media/ref --vcf-file /media/input.vcf.gz --directory /media/pipeline_data/real-data \
--assembly hg38 --phase --impute 
```

#### Finding relatives using IBIS

After preprocessing, there will be a file `data.vcf.gz` in directory /media/pipeline_data/real-data/preprocessed/ .
One have to use the same `--directory` value for both preprocessing and invoking search of relatives using `find` command.
 `--vcf-file` is ignored. Option `--flow ibis` invokes fast IBD estimation using IBIS software.  

```
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest \
launcher.py find --ref-directory /media/ref --directory /media/pipeline_data/real-data --flow ibis 
```


#### Description of output file

Output file is in .tsv file format, and contains one line for each detected pair of relatives.

```text
id1      id2     king_degree king_relation shared_genome_proportion kinship kinship_degree ersa_degree ersa_lower_bound ersa_upper_bound shared_ancestors final_degree  total_seg_len        seg_count
g1-b1-i1 g2-b1-i1     1           PO                0.4996          0.2493      1.0             1               1               1               0.0             1       3359.9362945470302      40
g1-b1-i1 g2-b2-i1     1           PO                0.4997          0.2459      1.0             1               1               1               0.0             1       3362.012253715002       40
g1-b1-i1 g2-b3-i1     1           PO                0.4999          0.2467      1.0             1               1               1               0.0             1       3363.150814131464       40
g1-b1-i1 g3-b1-i1     2           2                 0.2369          0.1163      2.0             2               2               2               1.0             2       1644.634182188072       60
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
   It should be approximately 0.5 for PO and FS, 0.25 for grandmother-granddaughter and half-siblings.
   For the first 3 degrees it is calculated as total len of IBD2 segments + half of total length of IBD1 segments. 
   For 4th+ degrees it is simply half of total length of IBD1 segments.  
 * `kinship` is the KING kinship coefficient.
 * `ersa_degree` is the degree estimated from IBD segments by ERSA, it is used for the `final_degree` in the cases where `king_degree` does not exist.
 * `ersa_lower_bound` is the lower bound degree estimation of ERSA using confidence interval 0.99,
  i.e. with probability (1-0.99)/2=0.005 degree will be lower than `ersa_lower_bound`.
 * `ersa_upper_bound` is the upper bound degree estimation of ERSA using confidence interval 0.99,
  i.e. with probability (1-0.99)/2=0.005 degree will be higher than `ersa_upper_bound`.
 * `shared_ancestors` is the most likeliest number of shared ancestors, if it is 0, then one relative is a direct descendant of the other, 
 if 1 then they probably have one common ancestor, i.e. half siblings, if 2 then they have common mother and father, for example.  
 * `final_degree` is simply `king_degree` for close relatives up to 3rd degree and `ersa_degree` for distant relatives.
 * `total_seg_len` is the total length of all IBD segments, for the first 3 degrees it is calculated using KING IBD data, 
   for the 4th+ degrees it is calculated using IBID or Germline IBD data.
 * `seg_count` is the total number of all IBD segments found by KING for the first 3 degrees and found by IBIS\Germline for the 4th+ degrees.

#### Description of some launcher parameters

`--zero-seg-count`, default=0.5
Average count of IBD segments in two unrelated individuals in population. 
Smaller values of 0.1, 0.2 tend to give more distant matches than default 0.5.

`--zero-seg-len`, default=5.0
Average length of IBD segment in two unrelated individuals in population. 
Smaller values of tend to give more distant matches than default 5.0

`--alpha`, default=0.01
ERSA P-value limit for testing for an existence of an relationship.
Values of 0.02-0.05 tend to give more distant matches that default 0.01. 

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

```bash
export SINGULARITY_TMPDIR=/media/tmp
export SINGULARITY_CACHEDIR=/media/tmp
```

You can also pass this values using Snakemake:
`--singularity-prefix DIR`
`--singularity-args ARGS`

6. Launch

```text
snakemake --cores all --use-conda --use-singularity --singularity-prefix=/media --singularity-args="-B /media:/media" -p all -n
```

Please mind `-n` flag for dry-run

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

```bash
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

```bash
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
