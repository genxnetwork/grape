# Accuracy check via simulation: IBIS workflow

| Info | Description |
|:--|:--|
| Test case ID  |   |
| Test priority  |   |
| Module name  | IBIS  |
| Test executed by  |   |
| Test execution date  |   |
| Description  | We are checking accuracy of IBIS workflow with `simulate` command. It takes roughly an hour.  |

## Pre-conditions

- /media is a parent directory to all data we use in the pipeline;
- /media/ref is a reference directory;
- /media/data is a working pipeline directory.

## Dependencies

## Steps

| Step number | Test step | Test result | Status |  Notes|
|:--|:--|:--|:--|:--|:--|
| 1  | enter command in console <br/>```docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py simulate --ref-directory /media/ref --cores 8 --directory /media/data --flow ibis --assembly hg37 --real-run``` | Recall of almost 100% for the first 3 degrees, around 95% for 4-6 degrees, and better than 0% for 7-9 degrees.  |  success |   | 
