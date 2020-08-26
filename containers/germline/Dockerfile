FROM debian:buster-slim

RUN apt-get update

RUN apt-get install -y make g++ git && \
    git clone https://github.com/gusevlab/germline && \
    cd /germline && \
    mkdir -p bin && \
    make all && \
    mv germline bin/germline && \
    apt-get --purge -y remove g++ make git && \
    apt-get -y autoremove

ENV PATH "$PATH:/germline/bin"
