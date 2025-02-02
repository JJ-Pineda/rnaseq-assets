FROM --platform=linux/amd64 ubuntu:latest

RUN apt update && \
    apt install -y --no-install-recommends \
        fastqc \
        multiqc \
        sambamba \
        bbmap \
        rna-star \
        dialog \
        hisat2 \
        gffread \
        salmon \
        rsem \
        subread \
        ca-certificates wget tar pigz unzip

RUN wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/root/miniconda3/bin:$PATH"

WORKDIR /
RUN wget -nv https://gitlab.com/genomeinformatics/xengsort/-/archive/2.0.5/xengsort-2.0.5.tar.gz && \
    tar xf xengsort-2.0.5.tar.gz && \
    rm xengsort-2.0.5.tar.gz

WORKDIR xengsort-2.0.5
RUN sed -i 's/numpy>=1.24/numpy>=1.24,<2.0/g' environment.yml

RUN conda init bash && \
    . ~/.bashrc && \
    conda env create && \
    conda activate xengsort && \
    pip install -e .

WORKDIR /
RUN wget -nv https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz && \
    tar -xf arriba_v2.4.0.tar.gz --exclude=arriba*/.git && \
    rm arriba_v2.4.0.tar.gz
ENV PATH="/arriba_v2.4.0:$PATH"

RUN wget -nv https://github.com/gpertea/stringtie/releases/download/v2.2.3/stringtie-2.2.3.Linux_x86_64.tar.gz && \
    tar -xf stringtie-2.2.3.Linux_x86_64.tar.gz && \
    mv stringtie-2.2.3.Linux_x86_64 stringtie-2.2.3 && \
    rm stringtie-2.2.3.Linux_x86_64.tar.gz
ENV PATH="/stringtie-2.2.3:$PATH"

WORKDIR /usr/share/java
RUN unzip -d bbmap bbmap.jar
ENV PATH="/usr/share/bbmap:$PATH"

WORKDIR /root
RUN mkdir javier

RUN wget -nv https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && \
    mv awscli-exe-linux-x86_64.zip awscliv2.zip && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

WORKDIR /root/javier
