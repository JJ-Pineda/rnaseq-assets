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
        ca-certificates wget tar

RUN wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

RUN /root/miniconda3/bin/conda create -n xengsort -c bioconda xengsort -y

RUN wget -nv https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz && \
    tar -xzf arriba_v2.4.0.tar.gz --exclude=arriba*/.git && \
    rm arriba_v2.4.0.tar.gz

RUN wget -nv https://github.com/gpertea/stringtie/releases/download/v2.2.3/stringtie-2.2.3.Linux_x86_64.tar.gz && \
    tar -xzf stringtie-2.2.3.Linux_x86_64.tar.gz && \
    mv stringtie-2.2.3.Linux_x86_64 stringtie-2.2.3 && \
    rm stringtie-2.2.3.Linux_x86_64.tar.gz

RUN echo 'PATH="$PATH:/root/miniconda3/bin:/arriba_v2.4.0:/usr/share/bbmap:/stringtie-2.2.3"' > ~/.bashrc

RUN cd /usr/share/java && unzip -d bbmap bbmap.jar

RUN cd && mkdir data scripts indexes
RUN cd /root/indexes && mkdir star xengsort salmon hisat2
