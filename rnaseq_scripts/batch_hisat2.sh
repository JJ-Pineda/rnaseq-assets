#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

HISAT2_INDEX=$1
STRANDNESS=$2   # Set to empty string (i.e. "") if unstranded
FASTQ_DIR=$3
READ1_SUFFIX=$4
READ2_SUFFIX=$5 # Leave blank for single-end reads

# Ensure that HISAT2 index exists
if [ -z "$(ls $HISAT2_INDEX)" ]
then
  echo "No HISAT2 index detected. Exiting..."
  exit
fi

SECONDS=0

if [ -z "$STRANDNESS" ]
then
  STRANDNESS_ARG=""
else
  STRANDNESS_ARG="--rna-strandness $STRANDNESS"
fi


cd "$FASTQ_DIR"
mkdir -p hisat2

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
