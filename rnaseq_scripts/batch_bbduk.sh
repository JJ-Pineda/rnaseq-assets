#!/bin/bash

FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3

SECONDS=0
REF_PATH="/Applications/bbmap/resources/adapters.fa"
TRIM_TAG="_trimmed_"

# Tell bash to abort on error
set -o pipefail
set -e -u

# Don't trim reads that have already been trimmed
cd "$FASTQ_DIR"
READ1_FILES=$(ls *$READ1_SUFFIX | grep -v $TRIM_TAG)

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"

  echo "Trimming reads for $BASE_NAME"

  # Check if variable is empty
  if [ -z "$READ2_SUFFIX" ]
  then
    IN_NAME=$BASE_NAME$READ1_SUFFIX
    OUT_NAME=$BASE_NAME$TRIM_TAG$READ1_SUFFIX

    # Adapter trimming arguments: ref=$REF_PATH ktrim=r k=23 mink=11 hdist=1
    # Quality trimming arguments should be avoided unless there's significantly low quality
    # In which case only trim on right end of read with weak trimming --> qtrim=r trimq=8
    bbduk.sh -Xmx4g threads=12 in=$IN_NAME out=$OUT_NAME ref=$REF_PATH ktrim=rl k=23 mink=11 hdist=1
  else
    IN_NAME1=$BASE_NAME$READ1_SUFFIX
    IN_NAME2=$BASE_NAME$READ2_SUFFIX
    OUT_NAME1=$BASE_NAME$TRIM_TAG$READ1_SUFFIX
    OUT_NAME2=$BASE_NAME$TRIM_TAG$READ2_SUFFIX

    # Benchmark: 12 threads --> ~5 minutes to process 1 paired-end file set
    # For paired-end reads, add "tpe tbo" arguments to trim the mates consistently
    bbduk.sh -Xmx4g threads=12 in1=$IN_NAME1 in2=$IN_NAME2 out1=$OUT_NAME1 out2=$OUT_NAME2 ref=$REF_PATH ktrim=rl k=23 mink=11 hdist=1 tpe tbo
  fi

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done
