FROM debian:sid-slim

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update

RUN apt-get install -y wget bcftools libgomp1

RUN wget ftp://share.sph.umich.edu/minimac3/Minimac3Executable.tar.gz -O Minimac3Executable.tar.gz

RUN tar -xzvf Minimac3Executable.tar.gz

ENV PATH "$PATH:/Minimac3Executable/bin"