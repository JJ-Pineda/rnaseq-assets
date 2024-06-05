#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# Note: genome-based --> "build 1"; transcriptome-based --> "build 2"
METHOD=$1 # "g" for genome-based index (i.e. default) or "t"

SECONDS=0

# For reference: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbsplit.sh
# Path to ENSEMBL files
HUMAN_PATH=/root/ensembl_references/grch38
MOUSE_PATH=/root/ensembl_references/grcm39
INDEX_PATH=/root/indexes/bbsplit/grch38_grcm39

# Build the directory if it doesn't exist
mkdir -p "$INDEX_PATH"

# Benchmark: <5 minutes
# Note: uses ~22gb of storage for a genome-based human index (much smaller for transcriptome-based)
if [ -z "$METHOD" ] || [ $METHOD = "g" ]
then
  HUMAN_GENOME=$(ls $HUMAN_PATH/*primary_assembly.fa.gz)
  MOUSE_GENOME=$(ls $MOUSE_PATH/*primary_assembly.fa.gz)

  bbsplit.sh -Xmx50g threads=12 build=1 path="$INDEX_PATH" ref_x="$HUMAN_GENOME" ref_y="$MOUSE_GENOME"
else
  HUMAN_TRANSCRIPTS=$(ls $HUMAN_PATH/*cdna.all.fa.gz)
  MOUSE_TRANSCRIPTS=$(ls $MOUSE_PATH/*cdna.all.fa.gz)
  bbsplit.sh -Xmx10g threads=12 build=2 path="$INDEX_PATH" ref_x="$HUMAN_TRANSCRIPTS" ref_y="$MOUSE_TRANSCRIPTS"
fi

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
