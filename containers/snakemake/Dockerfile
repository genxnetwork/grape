FROM debian:buster

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get clean && apt-get update && apt-get install -y ca-certificates libseccomp-dev squashfs-tools \
    && rm -rf /tmp/*

ADD envs/snakemake.yaml envs/snakemake.yaml

ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash

RUN apt-get install -y wget bzip2 gnupg2 git libgomp1 && \
    wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    conda install -c conda-forge mamba && \
    mamba env create -f envs/snakemake.yaml && \
    conda clean --all -y

# Intall Minimac3
RUN wget "https://bioinformatics.file.core.windows.net/bundles/Minimac3Executable.tar.gz?sv=2020-08-04&ss=f&srt=sco&sp=r&se=2022-08-08T14:35:53Z&st=2021-08-27T06:35:53Z&spr=https&sig=SjxrSn2KBuQYjYgT2ZZTHQ6IOhA%2BRUSvLIgog%2FH2Tnk%3D" -O /opt/Minimac3Executable.tar.gz
RUN tar -xzvf /opt/Minimac3Executable.tar.gz -C /opt/
RUN rm /opt/Minimac3Executable.tar.gz
ENV PATH "$PATH:/opt/Minimac3Executable/bin"

# Install Minimac4
RUN apt-get install -y minimac4

# Install Eagle
RUN wget "https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4.1" -O /usr/bin/eagle
RUN chmod a+x /usr/bin/eagle

# Install Germline
# conda version of germline has a problem with non-zero error code
# https://github.com/gusevlab/germline/issues/8
RUN apt-get install -y make g++ git && \
    git clone https://github.com/gusevlab/germline && \
    cd /germline && \
    mkdir -p bin && \
    make all && \
    mv germline bin/germline && \
    apt-get --purge -y remove g++ make git && \
    apt-get -y autoremove

# Workaround of NonWritableError when conda tries to create environments for the first time
# funnel launches docker containers with --read-only and snakemake cannot create conda envs
# because it has to do something with urls.txt
# https://github.com/conda/conda/issues/7267#issuecomment-458661530
RUN mkdir -p /root/.conda/pkgs && touch /root/.conda/pkgs/urls.txt

ADD . /src/repo
WORKDIR /src/repo

RUN echo "source activate snakemake" > ~/.bashrc
ENV PATH /opt/conda/envs/snakemake/bin:${PATH}

ENTRYPOINT ["python"]
