#!/bin/bash

FASTQ_DIR=$1
READ1_FILE=$2
READ2_FILE=$3

SECONDS=0

# Tell bash to abort on error
set -o pipefail
set -e -u

cd "$FASTQ_DIR"

# Check for paired-end file
if [ -z "$READ2_FILE" ]
then
  READ_ARGUMENT="-r $READ1_FILE"
else
  READ_ARGUMENT="-1 $READ1_FILE -2 $READ2_FILE"
fi

salmon quant \
  --index="/root/indexes/grch38" \
  --threads 12 \
  --libType A \
  --incompatPrior 0.0 \
  --output "${BASE_NAME}.salmon" \
  "$READ_ARGUMENT"

LIB_TYPE=$(egrep "Automatically detected most likely library type as [A-Z]+" "${BASE_NAME}.salmon/logs/salmon_quant.log" | rev | cut -d " " -f 1 | rev)

if [[ $LIB_TYPE == *"R"* ]]
then
  if [ ! -z "$READ2_FILE" ]
  then
    STRANDEDNESS="RF"
  else
  	STRANDEDNESS="R"
  fi
elif [[ $LIB_TYPE == *"F"* ]]
then
  if [ ! -z "$READ2_FILE" ]
  then
    STRANDEDNESS="FR"
  else
  	STRANDEDNESS="F"
  fi
else
  STRANDEDNESS="Unstranded"
fi

echo "Salmon detected the library type as '$LIBTYPE' --> Strandedness = '$STRANDEDNESS'"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
