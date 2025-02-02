#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

GENOME_ASSEMBLY=$1
TRANSCRIPTS=$2
INDEX_DIR=$3

SECONDS=0

mkdir -p "$INDEX_DIR"

grep "^>" <(gunzip -c "$GENOME_ASSEMBLY") | cut -d " " -f 1 > /tmp/decoys.txt
sed -i.bak -e 's/>//g' /tmp/decoys.txt

# Concatenate transcripts and genome targets (transcripts MUST come first)
cat "$TRANSCRIPTS" "$GENOME_ASSEMBLY" > /tmp/gentrome.fa.gz

# Build Salmon index (takes ~21 minutes and creates ~15gb worth of stuff)
# If using gencode references, "--gencode" flag is required for removing extra metadata in the target header
# "--threads" parameter specifies number of threads to use for index creation
# "-k" parameter specifies k-mer size (defaulted to 31 which is best for reads >=75bp)
salmon index \
  -t /tmp/gentrome.fa.gz \
  --decoys /tmp/decoys.txt \
  --threads 12 \
  --index "$INDEX_DIR"

rm /tmp/decoys.txt /tmp/decoys.txt.bak /tmp/gentrome.fa.gz

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
