#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# "grch38" or "grcm39"
GENOME=$1

INDEX_DIR="/root/indexes/salmon"
REF_DIR="/root/ensembl_references/$GENOME"

SECONDS=0

cd "$REF_DIR"

GENOME_ASSEMBLY=$(ls *primary_assembly.fa.gz)
TRANSCRIPTS=$(ls *cdna.all.fa.gz)

grep "^>" <(gunzip -c "$GENOME_ASSEMBLY") | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# Concatenate transcripts and genome targets (transcripts MUST come first)
cat "$TRANSCRIPTS" "$GENOME_ASSEMBLY" > gentrome.fa.gz

# Build Salmon index (takes ~21 minutes and creates ~15gb worth of stuff)
# If using gencode references, "--gencode" flag is required for removing extra metadata in the target header
# "--threads" parameter specifies number of threads to use for index creation
# "-k" parameter specifies k-mer size (defaulted to 31 which is best for reads >=75bp)
cd "$INDEX_DIR"

salmon index \
  -t "${REF_DIR}/gentrome.fa.gz" \
  --decoys "${REF_DIR}/decoys.txt" \
  --threads 12 \
  --index "$GENOME"

echo "Created Salmon index from $GENOME transcriptome using genome for decoys"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
