# Accuracy check via simulation: KING workflow

| Info | Description |
|:--|:--|
| Test case ID  |   |
| Test priority  |   |
| Module name  | KING  |
| Test executed by  |   |
| Test execution date  |   |
| Description  | We are checking accuracy of IBIS+KING workflow with `simulate` command and option `--flow ibis-king`. It takes roughly an hour.  |

## Pre-conditions

- /media is a parent directory to all data we use in the pipeline;
- /media/ref is a reference directory;
- /media/data is a working pipeline directory.

## Dependencies

## Steps


### Test Step № 1

Enter command in console

```bash
docker run --rm --privileged -it -v /media:/media -v /etc/localtime:/etc/localtime:ro genx_relatives:latest launcher.py simulate --ref-directory /media/ref --cores 8 --directory /media/data --flow ibis-king --assembly hg37 --real-run
```

#### Test result

After executing the command, the file  `metrics.tsv` containing metrics will appear:

| True Degree | Precision | Recall |
|:--|:--|:--|
| 1  | 1.0  | 1.0  |
| 2  | 1.0  | 1.0  |
| 3  | 0.9944444444444445  | 0.9944444444444445  |
| 4  | 0.9720670391061452 | 0.9666666666666667  |
| 5  | 0.9935064935064936  | 0.9444444444444444  |
| 6  | 0.952  | 0.8263888888888888  |
| 7  | 1.0  | 0.5793650793650794  |
| 8  | 0.9090909090909091  | 0.17094017094017094  |
| 9  | 1.0  | 0.027777777777777776  |
| 10  | 0.0  | 0.0  |
| 11  | 0.0  | 0.0  |
| 12  | 0.0  | 0.0  |
| 13  | 0.0  | 0.0  |

Recall of almost 100% for the first 3 degrees, around 95% for 4-6 degrees, and better than 0% for 7-9 degrees.

#### Status

Success

#### Notes




## Post-conditions
