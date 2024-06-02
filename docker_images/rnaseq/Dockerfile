FROM --platform=linux/amd64 public.ecr.aws/ubuntu/ubuntu:latest as installer

RUN apt update && \
    apt install -y --no-install-recommends \
        fastqc \
        multiqc \
        sambamba \
        bbmap \
        rna-star \
        dialog \
        hisat2 \
        salmon \
        rsem \
        subread \
        ca-certificates wget tar gcc g++ make

RUN wget -nv https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz && \
    tar -xzf arriba_v2.4.0.tar.gz --exclude=arriba*/.git && \
    rm arriba_v2.4.0.tar.gz

RUN cd arriba_v2.4.0 && make

RUN wget -nv https://github.com/gpertea/stringtie/releases/download/v2.2.3/stringtie-2.2.3.Linux_x86_64.tar.gz && \
    tar -xzf stringtie-2.2.3.Linux_x86_64.tar.gz && \
    mv stringtie-2.2.3.Linux_x86_64 stringtie-2.2.3 && \
    rm stringtie-2.2.3.Linux_x86_64.tar.gz

RUN echo 'export PATH="$PATH:/arriba_v2.4.0/"' > ~/.bashrc
RUN echo 'export PATH="$PATH:/usr/share/bbmap/"' >> ~/.bashrc
RUN echo 'export PATH="$PATH:/stringtie-2.2.3/"' >> ~/.bashrc

RUN cd /usr/share/java && unzip -d bbmap bbmap.jar

RUN cd && mkdir data scripts indexes
RUN cd /root/indexes && mkdir star bbsplit salmon hisat2
RUN cd /root/indexes/star && mkdir grch38
RUN cd /root/indexes/bbsplit && mkdir grch38_grcm39
RUN cd /root/indexes/salmon && mkdir grch38
