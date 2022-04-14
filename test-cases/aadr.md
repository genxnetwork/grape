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

| Step number | Test step | Expected result | Actual result | Status |  Notes|
|:--|:--|:--|:--|:--|:--|
| 1  | go to aadr: `cd /media/data/aadr` |   |   | success  |   |
| 2  | download dataset with the command <br/>```wget https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V50/V50.0/SHARE/public.dir/v50.0_1240K_public.tar``` |   |   | success  |   |
| 3  | unpack dataset with the command <br/>```tar -xvf v50.0_1240K_public.tar``` |   |   | success  |   |
| 4  | install converter from packed ancestry map format to plink ped with the command <br/>```conda install -c bioconda eigensoft``` |   |   | success  |   |
| 5  | convert to plink ped format with the command <br/>```convertf -p par.PACKEDANCESTRYMAP.PED```  |   |   | success  |   |
| 6  | convert to vcf format with the command <br/>```plink --ped v50.0_1240k_public.ped --map v50.0_1240k_public.pedsnp --alleleACGT --recode vcf-iid bgz --out aadr``` |   |   | success  |   |
| 7  | remove underscores from sample ids<br/>```bcftools query --list-samples aadr.vcf.gz > aadr.samples```<br/>```bcftools reheader aadr.vcf.gz -s aadr.samples \ bcftools view -O z -o aadr.reheaded.vcf.gz```  |   |   | success  |   |
| 8  | run preprocessing <br>```docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py preprocess --ref-directory /media/ref --cores 8 --directory /media/data/aadr --vcf-file /media/data/aadr/aadr.reheaded.vcf.gz --assembly hg37 --real-run```  |   |   | success |   |  
| 9  | run relationship inference <br/> ```docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py find --ref-directory /media/ref --cores 8 --directory /media/runs/aadr --flow ibis --real-run``` | `result.csv` file with no relatives between ancient samples and modern samples. Also, the running time of `preprocess` and `find` command.  |   |  success |   |

## Post-conditions
