FROM jupyter/datascience-notebook:ubuntu-22.04

USER root

WORKDIR /home/jovyan
RUN wget -nv https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && \
    mv awscli-exe-linux-x86_64.zip awscliv2.zip && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws
