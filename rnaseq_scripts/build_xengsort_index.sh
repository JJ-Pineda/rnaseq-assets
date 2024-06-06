#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

SECONDS=0

HUMAN_GENOME=$(ls /root/ensembl_references/grch38/*primary_assembly.fa.gz)
MOUSE_GENOME=$(ls /root/ensembl_references/grcm39/*primary_assembly.fa.gz)

INDEX_PATH=/root/indexes/xengsort/grch38_grcm39

# Build the directory if it doesn't exist
mkdir -p "$INDEX_PATH"

# Activate xengsort conda environment
source /root/miniconda3/etc/profile.d/conda.sh
conda activate xengsort

# Suggested number of 25-k-mers for human/mouse is 4.5 billion
# Likely possible to also create a transcriptome-based index (e.g. for downstream Salmon)
echo "Building xengsort index..."
xengsort index --index $INDEX_PATH/idx --host "$MOUSE_GENOME" --graft "$HUMAN_GENOME" -n 4_500_000_000 -k 25

conda deactivate

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
