#!/bin/bash

FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3
METHOD=$4

# Note: method should be either "t" for transcriptome (i.e. default) or "g" for genome
# In most RNAseq cases, "t" should be used, whereas "g" would be used when you need to look for novel transcripts
# But in those cases, you wouldn't be looking at data from a xenograft so BBSplit is probably not necessary at that point...
# For WES, "g" must be used

SECONDS=0

# For reference: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbsplit.sh
# Path to GENCODE files
HUMAN_PATH=/Users/javier/CompBioAssets/gencode_references/grch38
MOUSE_PATH=/Users/javier/CompBioAssets/gencode_references/grcm39
INDEX_PATH=/Users/javier/CompBioAssets/bbsplit_indices/grch38_grcm39

# Tell bash to abort on error
set -o pipefail
set -e -u

if [ -z "$METHOD" ] || [ $METHOD = "t" ]
then
  REF_X="${HUMAN_PATH}/gencode.v45.transcripts.fa.gz"
  REF_Y="${MOUSE_PATH}/gencode.vM34.transcripts.fa.gz"
else
  REF_X="${HUMAN_PATH}/GRCh38.primary_assembly.genome.fa.gz"
  REF_Y="${MOUSE_PATH}/GRCm39.primary_assembly.genome.fa.gz"
fi

# Alters the open-file limit
ulimit -n 10000

# Build BBSplit index (<5 minutes)
# Note: uses ~22gb of storage for a genome-based human index (much smaller for transcriptome-based)
bbsplit.sh -Xmx50g threads=12 build=1 path="$INDEX_PATH" ref_x="$REF_X" ref_y="$REF_Y"

# Perform the actual read splitting
cd "$FASTQ_DIR"
READ1_FILES=$(ls *$READ1_SUFFIX)

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"
  echo "Splitting reads for $BASE_NAME"

  # Check if variable is empty
  if [ -z "$READ2_SUFFIX" ]
  then
    IN_NAME=$BASE_NAME$READ1_SUFFIX

    bbsplit.sh -Xmx10g threads=12 build=1 path="$INDEX_PATH" in="$IN_NAME" basename="${BASE_NAME}_%.fastq.gz" refstats="${BASE_NAME}_stat.txt"

  else
    IN_NAME1=$BASE_NAME$READ1_SUFFIX
    IN_NAME2=$BASE_NAME$READ2_SUFFIX

    # Transcriptome benchmarks: 10gb seems to be sufficient
    # Benchmark: 12 threads --> ~48 minutes for one set of paired reads
    bbsplit.sh -Xmx10g threads=12 build=1 path="$INDEX_PATH" in1="$IN_NAME1" in2="$IN_NAME2" basename="${BASE_NAME}_%_#.fastq.gz" refstats="${BASE_NAME}_stat.txt"
  fi

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

cd "$INDEX_PATH"
rm -rf ref
