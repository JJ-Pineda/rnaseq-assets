#!/bin/bash

# Tell bash to abort on error
set -o pipefail
set -e -u

## STRANDNESS
# "RF" for reverse and paired-end
# "R" for reverse and single-end
# "FR" for forward and paired-end
# "F" for forward and single-end

RSEM_REF_DIR=$1
BAM_DIR=$2
BAM_SUFFIX=$3
STRANDNESS=$4

if [ $STRANDNESS = "RF" ]
then
  STRAND_ARG="--strandedness reverse"
  PAIRED_ARG="--paired-end"
elif [ $STRANDNESS = "FR" ]
then
  STRAND_ARG="--strandedness forward"
  PAIRED_ARG="--paired-end"
elif [ $STRANDNESS = "R" ]
then
  STRAND_ARG="--strandedness reverse"
  PAIRED_ARG=""
elif [ $STRANDNESS = "F" ]
then
  STRAND_ARG="--strandedness forward"
  PAIRED_ARG=""
else
  STRAND_ARG="--strandedness none"
  PAIRED_ARG=""
fi

SECONDS=0

# Unzip RSEM fasta files
gunzip -k $RSEM_REF_DIR/*fa.gz

# Perform actual transcriptome assembly
cd "$BAM_DIR"
BAM_FILES=$(ls *$BAM_SUFFIX)

for f in $BAM_FILES
do
  BASE_NAME="${f//$BAM_SUFFIX/}"
  echo "Running quantification for $BASE_NAME"

  rsem-calculate-expression \
    -p 8 $PAIRED_ARG $STRAND_ARG \
    --bam \
    --no-bam-output \
	--estimate-rspd \
	--append-names \
	"$f" \
	"${RSEM_REF_DIR}/rsem" \
	"$BASE_NAME"

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

# Remove decompressed RSEM fasta files
rm "$RSEM_REF_DIR"/*fa
