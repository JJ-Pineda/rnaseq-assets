#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

ASSEMBLY=$1    # Path to compressed genome assembly
ANNOTATION=$2  # Path to compressed genome annotation
OUTPUT_DIR=$3  # Path to output directory

SECONDS=0

gunzip -k "$ASSEMBLY" "$ANNOTATION"
ASSEMBLY=$(echo "$ASSEMBLY" | sed -r "s/.gz//g")
ANNOTATION=$(echo "$ANNOTATION" | sed -r "s/.gz//g")

mkdir -p "$OUTPUT_DIR"

# Fast: ~1-3 minutes
rsem-prepare-reference --gtf "$ANNOTATION" "$ASSEMBLY" "$OUTPUT_DIR/rsem"

# Gzip RSEM fasta files
gzip ${OUTPUT_DIR}/*.fa

# Remove uncompressed ensembl files
rm "$ASSEMBLY" "$ANNOTATION"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
