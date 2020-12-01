FROM debian:buster-slim


ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update

RUN apt-get install -y wget unzip && \
    wget -c http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip && \
    mkdir /plink && \
    mv plink2_linux_avx2.zip /plink/plink2_linux_avx2.zip && \
    cd plink && unzip plink2_linux_avx2.zip && mv plink2 /usr/bin/plink2 && \
    cd / && \
    chmod 755 /usr/bin/plink2

