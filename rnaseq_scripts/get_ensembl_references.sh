#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

GENOME=$1
ENSEMBL_RELEASE=$2

if [ -z "$ENSEMBL_RELEASE" ]
then
  ENSEMBL_RELEASE=112
fi

cd /root
mkdir -p ensembl_references
cd ensembl_references

if [ $GENOME = "grch38" ]
then
  mkdir -p grch38 && cd grch38
  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz
elif [ $GENOME = "grcm39" ]
then
  mkdir -p grcm39 && cd grcm39
  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz
fi
