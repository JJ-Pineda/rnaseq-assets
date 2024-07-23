#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# Here we are assuming that no downstream transcriptome assembly will be performed
# This script has also been written for human samples
STAR_INDEX=$1
FASTQ_DIR=$2
READ1_SUFFIX=$3
READ2_SUFFIX=$4

# Ensure that STAR index exists
if [ -z "$(ls $STAR_INDEX)" ]
then
  echo "No STAR index detected. Exiting..."
  exit
fi

SECONDS=0

cd "$FASTQ_DIR"

# Create directories for STAR and Arriba
STAR_OUT_DIR="star"
mkdir $STAR_OUT_DIR

# Adapted from "run_arriba.sh" script provided with Arriba installation
echo "Will use \"--twopassMode None\" for STAR alignment. This bash script will only retain transcriptomic coordinate BAM files."

READ1_FILES=$(ls *$READ1_SUFFIX)

echo -n > batch_star.log

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"

  echo "Processing reads for $BASE_NAME" >> batch_star.log

  # Check for paired-end file
  if [ -z "$READ2_SUFFIX" ]
  then
    READ_FILES="$BASE_NAME$READ1_SUFFIX"
  else
  	READ_FILES="$BASE_NAME$READ1_SUFFIX $BASE_NAME$READ2_SUFFIX"
  fi

  STAR \
    --runThreadN 8 \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad NoSharedMemory \
    --readFilesIn $READ_FILES \
    --readFilesCommand gunzip -c \
    --twopassMode None \
    --quantMode TranscriptomeSAM \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outBAMcompression 0 \
    --outFilterMultimapNmax 50 \
    --outFileNamePrefix "${STAR_OUT_DIR}/${BASE_NAME}_" \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 > "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.out.bam"

  # Remove unnecessary genomic coordinate BAM file
  rm "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.out.bam"

  # Filter and sort BAM file
  echo "Filtering and sorting transcriptomics coordinate BAM file" >> batch_star.log

  if [ -z "$READ2_SUFFIX" ]
  then
    sambamba view -F "not chimeric" -f bam --compression-level=0 "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.toTranscriptome.out.bam" |
    sambamba sort --sort-by-name -o "${STAR_OUT_DIR}/${BASE_NAME}_Filtered_Sorted.toTranscriptome.out.bam" /dev/stdin
  else
    # sambamba and samtools sometimes mess up multimapping paired-end reads when sorting --> use unix sort instead
    (
      sambamba view -F "not chimeric" --with-header "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.toTranscriptome.out.bam" |
      sort -S 24G --parallel=8 -t $'\t' -k 1,1 -k 13,13
    ) | sambamba view --sam-input -f bam -o "${STAR_OUT_DIR}/${BASE_NAME}_Filtered_Sorted.toTranscriptome.out.bam" /dev/stdin
  fi
  rm "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.toTranscriptome.out.bam"

  echo "Finished processing reads for $BASE_NAME" >> batch_star.log

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done
