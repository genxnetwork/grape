# Khazar Dataset

| Info | Description |
|:--|:--|
| Test case ID  |   |
| Test priority  |   |
| Module name  | GRAPE   |
| Test executed by  |   |
| Test execution date  |   |
| Description  | We are checking how GRAPE performs on the dataset with 314 samples with 350K SNPs in hg37 assembly. We perform preprocessing and search for SNPs.  |

## Pre-conditions

- /media is a parent directory to all data we use in the pipeline;
- /media/ref is a reference directory;
- /media/data is a working pipeline directory.

## Dependencies

## Steps

| Step number | Test step | Expected result | Actual result | Status |  Notes|
|:--|:--|:--|:--|:--|:--|
| 1  | download data with the command <br/>```wget https://evolbio.ut.ee/khazar/new_data_in_paper.{bed,bim,fam}``` |  |   | success  |   |
| 2  | replace _ in ids with -, because otherwise plink fails in our preprocessing <br/>```sed -i 's/_/-/g' new_data_in_paper.fam```  |   |   | success  |   |
| 3  | convert data to vcf.gz <br/>```plink --bfile new_data_in_paper --recode vcf-iid bgz --out /media/data/khazar/khazar314```  |   |   | success  |   |
| 4  | run preprocessing with the command <br/>```docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py preprocess --ref-directory /media/ref --cores 8 --directory /media/runs/khazar --vcf-file /media/data/khazar/khazar314.vcf.gz --assembly hg37 --real-run```  |   |   | success  |   |
| 5  | run relationship inference with the command <br/>```docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py find --ref-directory /media/ref --cores 8 --directory /media/runs/khazar --flow ibis --real-run```  | `result.csv` file with approximately 56 relatives.  |   | success  |   |

## Post-conditions
