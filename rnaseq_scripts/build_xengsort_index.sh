#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# "g" for genome-based index (i.e. default) or "t"
# In general, use "g"
# Only use "t" if you absolutely don't care about finding gene fusions or novel transcripts (i.e. using Salmon)
METHOD=$1

SECONDS=0

HUMAN_PATH=/root/javier/ensembl_references/grch38
MOUSE_PATH=/root/javier/ensembl_references/grcm39
INDEX_PATH=/root/javier/indexes/xengsort/grch38_grcm39

# Build the index directory if it doesn't exist
mkdir -p "$INDEX_PATH"

# Suggested number of 25-k-mers for human/mouse genome is 4.5 billion
# Using 210 million 25-k-mers for transcriptome index
if [ -z "$METHOD" ] || [ $METHOD = "g" ]
then
  HUMAN_REF=$(ls $HUMAN_PATH/*primary_assembly.fa.gz)
  MOUSE_REF=$(ls $MOUSE_PATH/*primary_assembly.fa.gz)
  INDEX_PATH=$INDEX_PATH/genome
  N_KMERS=4_500_000_000
else
  HUMAN_REF=$(ls $HUMAN_PATH/*cdna.all.fa.gz)
  MOUSE_REF=$(ls $MOUSE_PATH/*cdna.all.fa.gz)
  INDEX_PATH=$INDEX_PATH/transcriptome
  N_KMERS=210_000_000
fi

# Activate xengsort conda environment
source /root/miniconda3/etc/profile.d/conda.sh
conda activate xengsort

# Benchmark: ~15 minutes to build genome index
# Benchmark: ~1 minutes to build transcriptome index
xengsort index --index "$INDEX_PATH" --host "$MOUSE_REF" --graft "$HUMAN_REF" -n $N_KMERS -k 25

conda deactivate

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
