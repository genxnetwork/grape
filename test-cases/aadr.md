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

Go to aadr

```bash
cd /media/data/aadr
```

### Test Step № 2

Download dataset with the command

```bash
wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V50/V50.0/SHARE/public.dir/v50.0_1240K_public.tar
```

### Test Step № 3

Unpack dataset with the command

```bash
tar -xvf v50.0_1240K_public.tar
```


### Test Step № 4

Install converter from packed ancestry map format to plink ped with the command

```bash
conda install -c bioconda eigensoft
```
or
```bash
sudo apt install eigensoft
```

### Test Step № 5

Convert to plink ped format with the command

```bash
convertf -p par.PACKEDANCESTRYMAP.PED
```

par.PACKEDANCESTRYMAP.PED:
```bash
genotypename:    v50.0_1240k_public.geno
snpname:         v50.0_1240k_public.snp
indivname:       v50.0_1240k_public.ind
outputformat:    PED
genotypeoutname: v50.0_1240k_public.ped
snpoutname:      v50.0_1240k_public.pedsnp
indivoutname:    v50.0_1240k_public.pedind
```



### Test Step № 6

Convert to vcf format with the command

```bash
plink --ped v50.0_1240k_public.ped --map v50.0_1240k_public.pedsnp \
--alleleACGT --recode vcf-iid bgz --out aadr
```

### Test Step № 7

Remove underscores from sample `ids`

```bash
bcftools query --list-samples aadr.vcf.gz > aadr.samples
cat aadr.samples | while read line; do echo "${line//_}" >> aadr_clean.samples; done
bcftools reheader aadr.vcf.gz -s aadr_clean.samples | bcftools view -O z -o aadr.reheaded.vcf.gz
```


### Test Step № 8

Run preprocessing

```bash
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \
    grape:latest launcher.py preprocess \
        --ref-directory /media/ref --cores 8 --directory /media/data/aadr \
        --vcf-file /media/data/aadr/aadr.reheaded.vcf.gz --assembly hg37 \
        --het-samples 0.0 --real-run
```


### Test Step № 9

Run relationship inference

```bash
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \
    grape:latest launcher.py find --ref-directory /media/ref --cores 8 \
        --directory /media/runs/aadr --flow ibis --real-run
```

#### Test result

`result.csv` file with no relatives between ancient samples and modern samples. Also, the running time of `preprocess` and `find` command.

#### Status

Success

#### Notes


## Post-conditions
