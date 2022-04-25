# Performance test on 10K AADR dataset of IBIS workflow

| Info | Description |
|:--|:--|
| Test case ID  |   |
| Test priority  |   |
| Module name  |   |
| Test executed by  |   |
| Test execution date  |   |
| Description  | We need to be sure that GRAPE is working well with big datasets. We chose 10K dataset, because testing on 100K takes too long. Dataset is taken from Allen Ancient DNA Resource: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data. It has 10379 unique individuals (6442 ancient, 3937 present-day) with 1,233,013 SNPs.  |

## Pre-conditions

- /media is a parent directory to all data we use in the pipeline;
- /media/ref is a reference directory;
- /media/data is a working pipeline directory.

## Dependencies

## Steps


### Test Step № 1

go to aadr

```bash
c0d /media/data/aadr
```

#### Test result


#### Status

Success

#### Notes

### Test Step № 2

download dataset with the command

```bash
wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V50/V50.0/SHARE/public.dir/v50.0_1240K_public.tar
```

#### Test result

#### Status

Success

#### Notes

### Test Step № 3

unpack dataset with the command

```bash
tar -xvf v50.0_1240K_public.tar
```

#### Test result



#### Status

Success

#### Notes

### Test Step № 4

install converter from packed ancestry map format to plink ped with the command

```bash
conda install -c bioconda eigensoft
```

#### Test result


#### Status

Success

#### Notes

### Test Step № 5

convert to plink ped format with the command

```bash
convertf -p par.PACKEDANCESTRYMAP.PED
```

#### Test result



#### Status

Success

#### Notes


### Test Step № 6

convert to vcf format with the command

```bash
plink --ped v50.0_1240k_public.ped --map v50.0_1240k_public.pedsnp --alleleACGT --recode vcf-iid bgz --out aadr
```

#### Test result


#### Status

Success

#### Notes

### Test Step № 7

remove underscores from sample ids

```bash
bcftools query --list-samples aadr.vcf.gz > aadr.samples
bcftools reheader aadr.vcf.gz -s aadr.samples | bcftools view -O z -o aadr.reheaded.vcf.gz
```


#### Test result


#### Status

Success

#### Notes

### Test Step № 8

run preprocessing

```bash
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py preprocess --ref-directory /media/ref --cores 8 --directory /media/data/aadr --vcf-file /media/data/aadr/aadr.reheaded.vcf.gz --assembly hg37 --real-run
```

#### Test result


#### Status

Success

#### Notes

### Test Step № 9

run relationship inference

```bash
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py find --ref-directory /media/ref --cores 8 --directory /media/runs/aadr --flow ibis --real-run
```

#### Test result


#### Status

Success

#### Notes


## Post-conditions
