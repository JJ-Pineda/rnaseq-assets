#!/bin/bash

# "grch38" or "grcm39"
GENOME=$1

STAR_INDEX="/root/indexes/star/$GENOME"
REF_DIR="/root/ensembl_references/$GENOME"

mkdir "$STAR_INDEX"

SECONDS=0

# Tell bash to abort on error
set -o pipefail
set -e -u

ASSEMBLY=${REF_DIR}/$(cd "$REF_DIR" && ls *primary_assembly.fa.gz)
ANNOTATION=${REF_DIR}/$(cd "$REF_DIR" && ls *gtf.gz)

STAR --runThreadN 8 --runMode genomeGenerate --readFilesCommand "gunzip -c" --genomeDir "$STAR_INDEX" --genomeFastaFiles "$ASSEMBLY" --sjdbGTFfile "$ANNOTATION" --sjdbOverhang 100

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
