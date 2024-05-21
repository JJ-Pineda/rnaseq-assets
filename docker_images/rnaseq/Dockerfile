FROM public.ecr.aws/ubuntu/ubuntu:latest as installer

RUN apt update && \
    apt install -y --no-install-recommends \
        fastqc \
        multiqc \
        sambamba \
        bbmap \
        rna-star \
        stringtie \
        salmon \
        rsem \
        subread \
        ca-certificates wget tar gcc g++ make

RUN wget -nv 'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' && \
    tar -xzf 'arriba_v2.4.0.tar.gz' --exclude='arriba*/.git' && \
    rm 'arriba_v2.4.0.tar.gz'

RUN cd arriba_v2.4.0 && make

RUN echo 'PATH="$PATH:/arriba_v2.4.0/"' > ~/.bashrc

RUN cd && mkdir data && mkdir scripts
