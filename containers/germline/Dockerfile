FROM debian:buster

RUN apt-get update && apt-get install -y wget

RUN apt-get install -y make g++

RUN wget -q http://gusevlab.org/projects/germline/release/germline-1-5-3.tar.gz \
        -O /germline.tar.gz && \
        tar -xzf /germline.tar.gz&& \
        mv germline-1-5-3 /germline && \
        cd /germline && \
        make all

ENV PATH "$PATH:/germline/bin"
