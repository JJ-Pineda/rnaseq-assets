#!/bin/bash

FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3  # Leave blank for single-end reads

SECONDS=0

# Tell bash to abort on error
set -o pipefail
set -e -u

# Check for Salmon index
INDEX_PATH="/root/indexes/salmon"
cd "$INDEX_PATH"

# Build Salmon index if needed
if [ -z "$(ls -A $INDEX_PATH)" ]
then
  echo "No salmon index detected --> building salmon index"

  cd /root/scripts
  build_salmon_index.sh "grch38"
fi

# Grab files
cd "$FASTQ_DIR"
READ1_FILES=$(ls *$READ1_SUFFIX)

# Automatically detect library type for first file
# Then use that library type for the remaining files
LIB_TYPE="A"

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"
  IN_NAME1=$BASE_NAME$READ1_SUFFIX

  echo "Quantifying reads for $BASE_NAME"

  if [ -z "$READ2_SUFFIX" ]
  then
    READ_ARGUMENT="-r $IN_NAME1"
    echo "Using default '--fldMean 250' and '--fldStd 25' values for single-end fragment size. Correct these as needed!!!"
  else
    READ_ARGUMENT="-1 $IN_NAME1 -2 $BASE_NAME$READ2_SUFFIX"
  fi

  # Note: using incompatPrior=0 to better match RSEM quantification
  # I.e. throw out incompatible alignments even if they are the only ones for a given fragment
  # Important: if FastQC has determined that the samples have GC bias, add the "--gcBias" and "--seqBias" flags here
  salmon quant \
    --index="${INDEX_PATH}/grch38" \
    --threads 12 \
    --libType $LIB_TYPE \
    --incompatPrior 0.0 \
    --output "${BASE_NAME}.salmon" \
    $READ_ARGUMENT

  if [ $LIB_TYPE = "A" ]
  then
  	LIB_TYPE=$(egrep "Automatically detected most likely library type as [A-Z]+" "${BASE_NAME}.salmon/logs/salmon_quant.log" | rev | cut -d " " -f 1 | rev)
    echo "Will use $LIB_TYPE as the library type for all remaining files"
  fi

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

conda deactivate
