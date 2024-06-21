#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# "grch38" or "grcm39"
GENOME_ASSEMBLY=$1
ANNOTATION=$2
STAR_INDEX=$3

mkdir -p "$STAR_INDEX"

SECONDS=0

gunzip -k "$GENOME_ASSEMBLY" "$ANNOTATION"

GENOME_ASSEMBLY=$(echo "$GENOME_ASSEMBLY" | sed -r "s/.gz//g")
ANNOTATION=$(echo "$ANNOTATION" | sed -r "s/.gz//g")

# Takes ~30 minutes
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX" \
  --genomeFastaFiles "$GENOME_ASSEMBLY" \
  --sjdbGTFfile "$ANNOTATION" \
  --sjdbOverhang 100

rm "$GENOME_ASSEMBLY" "$ANNOTATION"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
