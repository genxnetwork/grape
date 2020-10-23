FROM openjdk:8
MAINTAINER Aleksandr Medvedev
# http://faculty.washington.edu/browning/refined-ibd/refined-ibd.17Jan20.102.jar

RUN apt-get update && \
    apt-get --no-install-recommends install -y --force-yes \
        wget && \
    apt-get clean autoclean && \
    apt-get autoremove -y

RUN mkdir /app && \
    wget http://faculty.washington.edu/browning/refined-ibd/refined-ibd.17Jan20.102.jar -O app/refined-ibd.jar

WORKDIR /app