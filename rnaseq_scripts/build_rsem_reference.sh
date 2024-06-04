#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# "grch38" or "grcm39"
GENOME=$1
REF_DIR="/root/ensembl_references/$GENOME"

SECONDS=0

ASSEMBLY=$(ls ${REF_DIR}/*primary_assembly.fa.gz)
ANNOTATION=$(ls ${REF_DIR}/*gtf.gz)

gunzip -k "$ASSEMBLY" "$ANNOTATION"
ASSEMBLY=$(echo "$ASSEMBLY" | sed -r "s/.gz//g")
ANNOTATION=$(echo "$ANNOTATION" | sed -r "s/.gz//g")

mkdir /root/rsem_references

RSEM_REF_DIR="/root/rsem_references/$GENOME"
mkdir "$RSEM_REF_DIR"

# Fast: ~1 minute
rsem-prepare-reference --gtf "$ANNOTATION" "$ASSEMBLY" "${RSEM_REF_DIR}/${GENOME}"

# Gzip RSEM fasta files
gzip ${RSEM_REF_DIR}/*.fa

# Remove uncompressed ensembl files
rm "$ASSEMBLY" "$ANNOTATION"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
