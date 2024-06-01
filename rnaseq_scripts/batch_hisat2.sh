#!/bin/bash

STRANDEDNESS=$1
FASTQ_DIR=$2
READ1_SUFFIX_GZ=$3
READ2_SUFFIX_GZ=$4

SECONDS=0

HISAT2_INDEX="/root/indexes/hisat2/grch38/genome"
READ1_SUFFIX=$(echo "$READ1_SUFFIX_GZ" | sed -r "s/.gz//g")
READ2_SUFFIX=$(echo "$READ2_SUFFIX_GZ" | sed -r "s/.gz//g")

# Tell bash to abort on error
set -o pipefail
set -e -u

cd "$FASTQ_DIR"
mkdir hisat2

READ1_FILES=$(ls *$READ1_SUFFIX_GZ)
for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX_GZ/}"
  echo "Processing reads for $BASE_NAME"

  # Check for paired-end files and decompress
  if [ -z "$READ2_SUFFIX_GZ" ]
  then
    gunzip -k "$BASE_NAME$READ1_SUFFIX_GZ"
    READ_ARGUMENT="-U $BASE_NAME$READ1_SUFFIX"
  else
    gunzip -k "$BASE_NAME$READ1_SUFFIX_GZ" "$BASE_NAME$READ2_SUFFIX_GZ"
  	READ_ARGUMENT="-1 $BASE_NAME$READ1_SUFFIX -2 $BASE_NAME$READ2_SUFFIX"
  fi

  # "--dta" is short for "--downstream-transcriptome-assembly" (i.e. StringTie)
  hisat2 -q -p 8 --rna-strandness "$STRANDEDNESS" -x "$HISAT2_INDEX" $READ_ARGUMENT --dta |
  sambamba view --sam-input --format bam --compression-level=0 /dev/stdin |
  sambamba sort -t 8 -o "${FASTQ_DIR}/hisat2/${BASE_NAME}_Sorted.bam" --compression-level=6 /dev/stdin

  # Remove decompressed read files
  rm "$BASE_NAME$READ1_SUFFIX" "$BASE_NAME$READ2_SUFFIX"

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done
