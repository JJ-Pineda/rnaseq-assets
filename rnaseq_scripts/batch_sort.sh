#!/bin/bash

BAM_DIR=$1
BAM_SUFFIX=$2

SECONDS=0

# Tell bash to abort on error
set -o pipefail
set -e -u

cd "$BAM_DIR"
BAM_FILES=$(ls *$BAM_SUFFIX)

for f in $BAM_FILES
do
  BASE_NAME="${f//$BAM_SUFFIX/}"
  echo "Sorting records for $BASE_NAME"
  samtools sort -@8 -m4g "$f" -o "${BASE_NAME}_Sorted${BAM_SUFFIX}"
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
