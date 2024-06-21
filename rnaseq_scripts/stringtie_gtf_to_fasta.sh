#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# "grch38" or "grcm39"
STRINGTIE_GTF=$1
GENOME_ASSEMBLY=$2
OUTPUT_DIR=$3

mkdir -p "$OUTPUT_DIR"

SECONDS=0

gunzip -k "$GENOME_ASSEMBLY" "$STRINGTIE_GTF"

GENOME_ASSEMBLY=$(echo "$GENOME_ASSEMBLY" | sed -r "s/.gz//g")
STRINGTIE_GTF=$(echo "$STRINGTIE_GTF" | sed -r "s/.gz//g")

gffread "$GENOME_ASSEMBLY" -g "$GENOME_ASSEMBLY" -w "${OUTPUT_DIR}/stringtie_transcriptome.fa"

rm "$GENOME_ASSEMBLY" "$STRINGTIE_GTF"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
