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

RUN echo 'export PATH="$PATH:/arriba_v2.4.0/"' > ~/.bashrc
RUN echo 'export PATH="$PATH:/usr/share/bbmap/"' >> ~/.bashrc

RUN cd /usr/share/java && unzip -d bbmap bbmap.jar

RUN cd && mkdir data scripts gencode_references
RUN cd gencode_references && mkdir grch38 grcm39
