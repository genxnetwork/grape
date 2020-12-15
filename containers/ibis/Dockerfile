FROM debian:buster-slim

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y wget zlib1g-dev libgomp1 && apt-get clean autoclean && apt-get autoremove -y

RUN wget https://github.com/williamslab/ibis/releases/download/1.20.9/ibis_release.tar.gz && \
    tar -xzvf ibis_release.tar.gz

ENV PATH "$PATH:/ibis_release"