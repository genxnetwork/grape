# List of test cases

### Accuracy check via simulation: IBIS workflow

**Description:**

We are checking accuracy of IBIS workflow with `simulate` command. 
It takes roughly an hour. 

**Command:**

```console
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py simulate --ref-directory /media/ref --cores 8 --directory /media/data --flow ibis --assembly hg37 --real-run
```

**Desired result:**

Recall of almost 100\% for the first 3 degrees, around 95\% for 4-6 degrees, and better than 0\% for 7-9 degrees.

### Accuracy check via simulation: KING workflow

**Description:**

We are checking accuracy of IBIS+KING workflow with `simulate` command and option `--flow ibis-king`. 
It takes roughly an hour. 

**Command:**

```console
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py simulate --ref-directory /media/ref --cores 8 --directory /media/data --flow ibis-king --assembly hg37 --real-run
```

**Desired result:**

Recall of almost 100\% for the first 3 degrees, around 95\% for 4-6 degrees, and better than 0\% for 7-9 degrees.

### Khazar Dataset

**Description:**

We are checking how GRAPE performs on the dataset with 314 samples with 350K SNPs in hg37 assembly. 
We perform preprocessing and search for SNPs. 

**Command:**

```console
# download data
wget https://evolbio.ut.ee/khazar/new_data_in_paper.{bed,bim,fam}

# replace _ in ids with -, because otherwise plink fails in our preprocessing 
sed -i 's/_/-/g' new_data_in_paper.fam

# convert data to vcf.gz
plink --bfile new_data_in_paper --recode vcf-iid bgz --out /media/data/khazar/khazar314 

# run preprocessing
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py preprocess --ref-directory /media/ref --cores 8 --directory /media/runs/khazar --vcf-file /media/data/khazar/khazar314.vcf.gz --assembly hg37 --real-run

# run relationship inference
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py find --ref-directory /media/ref --cores 8 --directory /media/runs/khazar --flow ibis --real-run
```

**Desired result:**

`result.csv` file with approximately 56 relatives.

### Performance test on 10K AADR dataset of IBIS workflow

**Description:**

We need to be sure that GRAPE is working well with big datasets. We chose 10K dataset, because testing on 100K takes too long.
Dataset is taken from Allen Ancient DNA Resource: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data.
It has 10379 unique individuals (6442 ancient, 3937 present-day) with 1,233,013 SNPs. 

**Command:**

```console
cd /media/data/aadr

# download dataset
wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V50/V50.0/SHARE/public.dir/v50.0_1240K_public.tar

# unpack dataset
tar -xvf v50.0_1240K_public.tar

# install converter from packed ancestry map format to plink ped
conda install -c bioconda eigensoft

# convert to plink ped format
convertf -p par.PACKEDANCESTRYMAP.PED

# convert to vcf format
plink --ped v50.0_1240k_public.ped --map v50.0_1240k_public.pedsnp --alleleACGT --recode vcf-iid bgz --out aadr

# remove underscores from sample ids
bcftools query --list-samples aadr.vcf.gz > aadr.samples
bcftools reheader aadr.vcf.gz -s aadr.samples | bcftools view -O z -o aadr.reheaded.vcf.gz

# run preprocessing
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py preprocess --ref-directory /media/ref --cores 8 --directory /media/data/aadr --vcf-file /media/data/aadr/aadr.reheaded.vcf.gz --assembly hg37 --real-run

# run relationship inference
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py find --ref-directory /media/ref --cores 8 --directory /media/runs/aadr --flow ibis --real-run
```

**Desired result:**

`result.csv` file with no relatives between ancient samples and modern samples. Also, the running time of `preprocess` and `find` command. 
