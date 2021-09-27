FROM golang:1.14.0-buster as builder

ARG SINGULARITY_TAG="v3.7.1"

ENV DEBIAN_FRONTEND noninteractive

WORKDIR $GOPATH/src/github.com/sylabs
#RUN apt-get update && apt-get install -y gawk gcc git libc-dev linux-headers libressl-dev libuuid libseccomp-dev make util-linux-dev
RUN apt-get update && apt-get install -y gawk gcc git libseccomp-dev make
RUN git clone https://github.com/sylabs/singularity.git \
    && cd singularity \
    && git checkout "tags/$SINGULARITY_TAG" \
    && ./mconfig -p /usr/local/singularity \
    && cd builddir \
    && make \
    && make install
RUN apt-get clean

FROM debian:buster

ENV DEBIAN_FRONTEND noninteractive

COPY --from=builder /usr/local/singularity /usr/local/singularity
ENV PATH="/usr/local/singularity/bin:$PATH" \
    SINGULARITY_TMPDIR="/tmp-singularity"
RUN apt-get clean && apt-get update && apt-get install -y ca-certificates libseccomp-dev squashfs-tools \
    && rm -rf /tmp/* \
    && mkdir -p $SINGULARITY_TMPDIR

#WORKDIR /work

#ENTRYPOINT ["/usr/local/singularity/bin/singularity"]

ADD envs/snakemake.yaml envs/snakemake.yaml

ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash

#ENV CONDA_PKGS_DIR /tmp/conda/pgks

RUN apt-get install -y wget bzip2 gnupg2 git openjdk-11-jre openjdk-11-jdk libgomp1 && \
    wget -nv https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    conda install -c conda-forge mamba && \
    mamba env create -f envs/snakemake.yaml && \
    conda clean --all -y

RUN wget https://bioinformatics.file.core.windows.net/bundles/Minimac3Executable.tar.gz?sv=2020-08-04&ss=f&srt=sco&sp=r&se=2022-08-08T14:35:53Z&st=2021-08-27T06:35:53Z&spr=https&sig=SjxrSn2KBuQYjYgT2ZZTHQ6IOhA%2BRUSvLIgog%2FH2Tnk%3D -O Minimac3Executable.tar.gz

RUN tar -xzvf Minimac3Executable.tar.gz

ENV PATH "$PATH:/Minimac3Executable/bin"

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