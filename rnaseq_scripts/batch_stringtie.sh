#!/bin/bash

BAM_DIR=$1
BAM_SUFFIX=$2
STRANDNESS=$3
ENSEMBL_DIR=$4

SECONDS=0

if [ -z "$STRANDNESS" ]
then
  STRANDNESS=""
fi

if [ -z "$ENSEMBL_DIR" ]
then
  ENSEMBL_DIR="/root/ensembl_references/grch38"
fi

if [ -z "$ENSEMBL_DIR" ]
then
  ENSEMBL_DIR="/root/ensembl_references/grch38"
fi

ANNOTATION=${ENSEMBL_DIR}/$(cd "$ENSEMBL_DIR" && ls *gtf.gz)

# Tell bash to abort on error
set -o pipefail
set -e -u

# Unzip GTF
gunzip -k "$ANNOTATION"
ANNOTATION=$(echo "$ANNOTATION" | sed -r "s/.gz//g")

# Perform actual transcriptome assembly
cd "$BAM_DIR"
BAM_FILES=$(ls *$BAM_SUFFIX)

for f in $BAM_FILES
do
  BASE_NAME="${f//$BAM_SUFFIX/}"
  echo "Assembling transcriptome for $BASE_NAME"
  stringtie -p 8 -G "$ANNOTATION" -o "${BASE_NAME}_stringtie.gtf" $STRANDNESS "$f"
done

# Merge StringTie outputs
STRINGTIE_FILES=$(ls *stringtie.gtf | tr '\n' ' ' | xargs)
stringtie --merge -G "$ANNOTATION" -o "stringtie_transcriptome.gtf" "$STRINGTIE_FILES"

# Remove decompressed annotation file and intermediate stringtie GTF files
rm "$ANNOTATION" *stringtie.gtf

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
