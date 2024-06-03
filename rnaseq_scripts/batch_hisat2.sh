#!/bin/bash

STRANDNESS=$1
FASTQ_DIR=$2
READ1_SUFFIX=$3
READ2_SUFFIX=$4

SECONDS=0

HISAT2_INDEX="/root/indexes/hisat2/grch38/genome"

# Tell bash to abort on error
set -o pipefail
set -e -u

if [ -z "$STRANDNESS" ]
then
  STRANDNESS_ARG=""
else
  STRANDNESS_ARG="--rna-strandness $STRANDNESS"
fi


cd "$FASTQ_DIR"
mkdir hisat2

READ1_FILES=$(ls *$READ1_SUFFIX)
for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"
  echo "Processing reads for $BASE_NAME"

  # Check for paired-end files and decompress
  if [ -z "$READ2_SUFFIX" ]
  then
    READ_ARGUMENT="-U $BASE_NAME$READ1_SUFFIX"
  else
    READ_ARGUMENT="-1 $BASE_NAME$READ1_SUFFIX -2 $BASE_NAME$READ2_SUFFIX"
  fi

  # "--dta" is short for "--downstream-transcriptome-assembly" (i.e. StringTie)
  hisat2 -q -p 8 $STRANDNESS_ARG -x "$HISAT2_INDEX" $READ_ARGUMENT --dta |
  sambamba view --sam-input --format bam --compression-level=0 /dev/stdin |
  sambamba sort -t 8 -o "${FASTQ_DIR}/hisat2/${BASE_NAME}_Sorted.bam" /dev/stdin

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done
