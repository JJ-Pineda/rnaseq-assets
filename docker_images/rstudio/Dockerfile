FROM rocker/rstudio

RUN apt update && apt install -y --no-install-recommends wget unzip

RUN if [ $(dpkg --print-architecture) = "amd64" ]; then \
         wget -nv https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && \
         mv awscli-exe-linux-x86_64.zip awscliv2.zip; \
    else \
         wget https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip && \
         mv awscli-exe-linux-aarch64.zip awscliv2.zip; \
    fi

RUN unzip awscliv2.zip && ./aws/install && rm -rf awscliv2.zip aws

ENV PASSWORD=passw0rd
ENV ROOT=true