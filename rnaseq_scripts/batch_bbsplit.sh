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
# Path to ENSEMBL files
HUMAN_PATH=/root/ensembl_references/grch38
MOUSE_PATH=/root/ensembl_references/grcm39
INDEX_PATH=/root/indexes/bbsplit/grch38_grcm39

# Tell bash to abort on error
set -o pipefail
set -e -u

# Build BBSplit index (<5 minutes)
# Note: uses ~22gb of storage for a genome-based human index (much smaller for transcriptome-based)
if [ -z "$METHOD" ] || [ $METHOD = "t" ]
then
  HUMAN_TRANSCRIPTS=$(cd "$HUMAN_PATH" && ls *cdna.all.fa.gz)
  MOUSE_TRANSCRIPTS=$(cd "$MOUSE_PATH" && ls *cdna.all.fa.gz)
  bbsplit.sh -Xmx10g threads=12 build=1 path="$INDEX_PATH" ref_x="$HUMAN_TRANSCRIPTS" ref_y="$MOUSE_TRANSCRIPTS"
else
  HUMAN_GENOME=$(cd "$HUMAN_PATH" && ls *primary_assembly.fa.gz)
  MOUSE_GENOME=$(cd "$MOUSE_PATH" && ls *primary_assembly.fa.gz)

  # Note: "maxindel" parameter set to "200k" when performing for genome sequencing
  bbsplit.sh -Xmx50g threads=12 build=1 maxindel=200k path="$INDEX_PATH" ref_x="$HUMAN_GENOME" ref_y="$MOUSE_GENOME"
fi

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
    # Benchmark: 12 threads over 11 CPUs --> ~21 minutes for one set of paired reads (including index creation)
    bbsplit.sh -Xmx10g threads=12 build=1 path="$INDEX_PATH" in1="$IN_NAME1" in2="$IN_NAME2" basename="${BASE_NAME}_%_#.fastq.gz" refstats="${BASE_NAME}_stat.txt"
  fi

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

cd "$INDEX_PATH"
rm -rf ref
