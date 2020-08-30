FROM debian:buster-slim

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update

RUN apt-get install -y make g++ git libboost-all-dev libz-dev && \
    git clone https://github.com/williamslab/ped-sim && \
    cd /ped-sim && \
    make && \
    mkdir -p bin && \
    mv ped-sim bin/pedsim && \
    apt-get --purge -y remove g++ make git libboost-all-dev && \
    apt-get -y autoremove

ENV PATH "$PATH:/ped-sim/bin"
