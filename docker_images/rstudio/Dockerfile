FROM --platform=linux/amd64 rocker/rstudio
USER root

RUN apt update && \
    apt install -y --no-install-recommends \
        wget \
        unzip

WORKDIR /root
RUN wget -nv https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && \
    mv awscli-exe-linux-x86_64.zip awscliv2.zip && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

RUN mkdir /root/javier
WORKDIR /root/javier

ENV PASSWORD=passw0rd
