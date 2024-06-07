#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3
METHOD=$4

# Note: method should be either "g" for genome (i.e. default), or "t" for transcriptome
# If using Salmon downstream, "t" should be used whereas "g" should be used otherwise
# The genome based method should especially be used if you care about novel transcripts or gene fusions
# For WES, "g" must be used

SECONDS=0

# For reference: https://github.com/BioInfoTools/BBMap/blob/master/sh/bbsplit.sh
INDEX_PATH=/root/indexes/bbsplit/grch38_grcm39

if [ -z "$METHOD" ] || [ $METHOD = "g" ]
then
  BUILD=1
  MAXINDEL=200000
  AMBIG2="all"
  ./build_bbsplit_index.sh g
else
  BUILD=2
  MAXINDEL=20   # i.e. BBSplit's default
  AMBIG2="best" # i.e. BBSplit's default
  ./build_bbsplit_index.sh t
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

    bbsplit.sh \
      -Xmx10g \
      threads=12 \
      build=$BUILD \
      path="$INDEX_PATH" \
      maxindel=$MAXINDEL \
      ambiguous2=$AMBIG2 \
      in="$IN_NAME" \
      basename="${BASE_NAME}_%.fastq.gz" \
      refstats="${BASE_NAME}_stat.txt"
  else
    IN_NAME1=$BASE_NAME$READ1_SUFFIX
    IN_NAME2=$BASE_NAME$READ2_SUFFIX

    # Transcriptome benchmarks: 10gb seems to be sufficient
    # Benchmark: 12 threads over 11 CPUs --> ~21 minutes for one set of paired reads (including index creation)
    bbsplit.sh \
      -Xmx10g \
      threads=12 \
      build=$BUILD \
      path="$INDEX_PATH" \
      maxindel=$MAXINDEL \
      ambiguous2=$AMBIG2 \
      in1="$IN_NAME1" \
      in2="$IN_NAME2" \
      basename="${BASE_NAME}_%_#.fastq.gz" \
      refstats="${BASE_NAME}_stat.txt"
  fi

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

# Remove BBSplit index (very quick to build and takes up a LOT of disk space)
rm -rf ref
