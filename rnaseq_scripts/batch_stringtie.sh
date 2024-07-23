#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

BAM_DIR=$1
BAM_SUFFIX=$2
STRANDNESS=$3   # e.g. "--rf" or "--fr"
ENSEMBL_DIR=$4

SECONDS=0

if [ -z "$STRANDNESS" ]
then
  STRANDNESS=""
fi

if [ -z "$ENSEMBL_DIR" ]
then
  ENSEMBL_DIR="/root/javier/ensembl_references/grch38"
fi

ANNOTATION=$(ls ${ENSEMBL_DIR}/*gtf.gz)

# Unzip GTF
gunzip -k "$ANNOTATION"
ANNOTATION=$(echo "$ANNOTATION" | sed -r "s/.gz//g")

# Perform actual transcriptome assembly
cd "$BAM_DIR"
BAM_FILES=$(ls *$BAM_SUFFIX)

echo -n > batch_stringtie.log

for f in $BAM_FILES
do
  BASE_NAME="${f//$BAM_SUFFIX/}"
  echo "Assembling transcriptome for $BASE_NAME" >> batch_stringtie.log

  stringtie -p 8 -G "$ANNOTATION" -o "${BASE_NAME}_stringtie.gtf" $STRANDNESS "$f"

  echo "Finished assembling transcriptome for $BASE_NAME" >> batch_stringtie.log
done

# Merge StringTie outputs
echo "Merging StringTie transcriptomes" >> batch_stringtie.log

STRINGTIE_FILES=$(ls *stringtie.gtf | tr '\n' ' ' | xargs)
stringtie --merge -G "$ANNOTATION" -o "stringtie_transcriptome.gtf" "$STRINGTIE_FILES"

# Check number of lines in the merged GTF that are missing strand information
UNSTRANDED=$(awk -F'\t' '$7 == "."' stringtie_transcriptome.gtf | wc -l)
STRANDED=$(awk -F'\t' '$7 == "+" || $7 == "-"' stringtie_transcriptome.gtf | wc -l)
LINE_COUNT=$(wc -l stringtie_transcriptome.gtf)
echo "$UNSTRANDED lines are missing strand information (out of $LINE_COUNT)." >> batch_stringtie.log
echo "$STRANDED lines contain strand information (out of $LINE_COUNT)." >> batch_stringtie.log
echo "Will remove the $UNSTRANDED offending lines that are missing strand information." >> batch_stringtie.log

awk '$7 == "+" || $7 == "-"' stringtie_transcriptome.gtf > stringtie_transcriptome_filtered.gtf

# Remove decompressed annotation file and intermediate stringtie GTF files
rm "$ANNOTATION" *stringtie.gtf stringtie_transcriptome.gtf

# Finally, compress final annotation file
gzip stringtie_transcriptome_filtered.gtf

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
