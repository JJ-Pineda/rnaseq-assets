#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

# Note: method should be either "g" for genome (i.e. default), or "t" for transcriptome
# If using Salmon downstream, "t" should be used whereas "g" should be used otherwise
# The genome based method should especially be used if you care about novel transcripts or gene fusions
# For WES, "g" must be used
FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3
METHOD=$4

SECONDS=0

INDEX_PATH=/root/indexes/xengsort/grch38_grcm39
if [ -z "$METHOD" ] || [ $METHOD = "g" ]
then
  METHOD="g"
  INDEX_PATH=$INDEX_PATH/genome
else
  INDEX_PATH=$INDEX_PATH/transcriptome
fi

mkdir -p "$INDEX_PATH"
if [ -z "$(ls $INDEX_PATH)" ]
then
  ./build_xengsort_index.sh "$METHOD"
fi

cd "$FASTQ_DIR"
READ1_FILES=$(ls *$READ1_SUFFIX)

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"
  echo "Splitting reads for $BASE_NAME"

  gunzip -k "$BASE_NAME"*
  READ_ARG="--fastq $f"
  if [ -n "$READ2_SUFFIX" ]
  then
    READ_ARG="$READ_ARG --pairs $BASE_NAME$READ2_SUFFIX"
  fi

  source /root/miniconda3/etc/profile.d/conda.sh
  conda activate xengsort

  # Note: can only specify threads for the "classify" method not "index"
  xengsort classify --index "$INDEX_PATH" $READ_ARG --prefix "$BASE_NAME" --threads 8

  conda deactivate

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

#rm -rf "$INDEX_PATH"
