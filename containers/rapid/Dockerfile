FROM ubuntu:16.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update

RUN apt-get install -y wget && \
    wget https://github.com/ZhiGroup/RaPID/raw/master/RaPID_v.1.7 -O rapid && \
    chmod u+x rapid && \
    mkdir app && \
    mv rapid /app/ && \
    apt-get --purge -y remove wget && \
    apt-get -y autoremove

ENV PATH "$PATH:/app"
