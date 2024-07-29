FROM rocker/rstudio

# Extra packages necessary for proper Bioconductor functionality
RUN apt update && apt install -y --no-install-recommends \
    zlib1g-dev liblzma-dev libbz2-dev libxml2-dev libpng-dev

# Install Bioconductor package manager
RUN Rscript -e 'install.packages("BiocManager", lib="/usr/local/lib/R/site-library", repos="https://cloud.r-project.org")'

RUN Rscript -e 'BiocManager::install(c("Biobase", "AnnotationHub", "ensembldb", "DESeq2", "DEXSeq"), lib="/usr/local/lib/R/site-library")'

# Install AWS CLI
RUN apt install -y --no-install-recommends wget unzip

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