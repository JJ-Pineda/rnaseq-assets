#!/bin/bash

BAM_DIR=$1
BAM_SUFFIX=$2

SECONDS=0

GENCODE_DIR="/root/gencode_references/grch38"
GENCODE_FILE="gencode.v46.primary_assembly.annotation.gtf.gz"

# Tell bash to abort on error
set -o pipefail
set -e -u

# Gunzip annotation
cd "$GENCODE_DIR"
gunzip -k "$GENCODE_FILE"
GENCODE_FILE=$(echo "$GENCODE_FILE" | sed -r "s/.gz//g")

# Remove header info from reference GTF
echo "test 1"
FIRST_LINE=$(grep -n "chr1" "$GENCODE_FILE" | head -n 1 | cut -d: -f1)
echo "test 2"
tail -n "+$FIRST_LINE" "$GENCODE_FILE" > "truncated_gencode_annotation.gtf"
echo "test 3"
mv "truncated_gencode_annotation.gtf" "$GENCODE_FILE"
echo "test 4"
ANNOTATION="${GENCODE_DIR}/${GENCODE_FILE}"

cd "$BAM_DIR"
BAM_FILES=$(ls *$BAM_SUFFIX)

for f in $BAM_FILES
do
  BASE_NAME="${f//$BAM_SUFFIX/}"
  echo "Assembling transcriptome for $BASE_NAME"
  stringtie -p 8 -G "$ANNOTATION" -o "${BASE_NAME}_stringtie.gtf" "$f"
done

# Merge StringTie outputs
STRINGTIE_FILES=$(ls *stringtie.gtf | tr '\n' ' ' | xargs)
stringtie -G "$ANNOTATION" -o "stringtie_transcriptome.gtf" "$STRINGTIE_FILES"

# Remove decompressed annotation file and intermediate stringtie GTF files
rm "$ANNOTATION" *stringtie.gtf

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
