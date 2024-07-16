#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

GENOME=$1

cd /root/javier
mkdir -p gencode_references
cd gencode_references

if [ $GENOME = "grch38" ]
then
  mkdir -p grch38 && cd grch38
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
elif [ $GENOME = "grcm39" ]
then
  mkdir -p grcm39 && cd grcm39
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.primary_assembly.genome.fa.gz
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.transcripts.fa.gz
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.primary_assembly.annotation.gtf.gz
fi
