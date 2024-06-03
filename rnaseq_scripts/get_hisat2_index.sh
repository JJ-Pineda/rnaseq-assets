#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# Takes ~40 minutes
cd /root/indexes/hisat2
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar xvf grch38_genome.tar.gz
rm grch38_genome.tar.gz
