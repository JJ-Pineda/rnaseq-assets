#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

FASTQ_DIR=$1
FASTQ_SUFFIX=$2
OUTPUT_FOLDER=$3

if [ -z "$FASTQ_DIR" ] || [ -z "$FASTQ_SUFFIX" ] || [ -z "$OUTPUT_FOLDER" ]; then
  echo "Either fastq data path, fastq suffix, or output folder is missing"
  exit 1
fi

SECONDS=0

cd "$FASTQ_DIR"
mkdir -p "$OUTPUT_FOLDER"
FASTQ_FILES=$(ls *$FASTQ_SUFFIX | tr '\n' ' ')
FASTQ_FILES=$(echo $FASTQ_FILES)

# Benchmark: 8 threads --> ~7 minutes to process 8 PDMR fastq files (i.e. 4 paired-end sets)
fastqc $FASTQ_FILES -t 8 -o "${OUTPUT_FOLDER}/" > batch_fastqc.log

multiqc "${OUTPUT_FOLDER}/" -o "${OUTPUT_FOLDER}/"
 
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
