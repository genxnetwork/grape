=======================================
Installation
=======================================

Requirements
----------------------------------

**GRAPE** requires a 64-bit version of Linux.

Total amount of required disk space must be at least **50GB**.

Installation
----------------------------------

Download the `latest release <https://github.com/genxnetwork/releases>`_ and extract it to your selected
directory.

Alternatively, you can clone the most recent (unstable) version from the
GitHub repository:

::

    git clone https://github.com/genxnetwork/grape.git

Build docker container from the GRAPE directory:
::

    docker build -t grape:latest -f containers/snakemake/Dockerfile -m 8GB .


Download Reference Datasets
----------------------------------
Reference datasets are required for correct pipeline operataion.

To download reference dataset one should use `reference` command of the pipeline launcher.

The command below downloads needed references to the directory specified by the flag `--ref-directory`.
After that `--ref-directory` argument should be used in all subsequent commands.

::

    docker run --rm -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \
    grape:latest launcher.py reference --ref-directory /media/ref --real-run

If phasing and imputation are required (mainly for the GERMLINE workflow) one should specify additional `--phase` and `--impute` flags to previous command to download additional reference datasets.

::

    docker run --rm -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \
    grape:latest launcher.py reference \
        --ref-directory /media/ref --phase --impute --real-run


There are another options to download all required reference data as a single file.
This file is prepared by the GRAPE developers and preloaded on our side in the cloud.

It can be done by specifying additional flag `--use-bundle` to the `reference` command.
This way is faster, since all the post-processing procedures have been already performed.

::

    docker run --rm -it -v /media:/media -v /etc/localtime:/etc/localtime:ro \
    grape:latest launcher.py reference --use-bundle \
        --ref-directory /media/ref --phase --impute --real-run

.. note::
    This step only needs to be done once.
