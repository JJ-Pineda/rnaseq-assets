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
FILE_EXTENSION=$4  # E.g. "FASTQ.gz", "fa.gz", "fq.gz", etc
METHOD=$5

SECONDS=0

INDEX_PATH=/root/indexes/xengsort/grch38_grcm39
if [ -z "$METHOD" ] || [ $METHOD = "g" ]
then
  METHOD="g"
  INDEX_PATH=$INDEX_PATH/genome
else
  INDEX_PATH=$INDEX_PATH/transcriptome
fi

if [ -z "$(ls $INDEX_PATH*)" ]
then
  ./build_xengsort_index.sh "$METHOD"
fi

cd "$FASTQ_DIR"
READ1_FILES=$(ls *$READ1_SUFFIX)

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"
  echo "Splitting reads for $BASE_NAME"

  if [ $FILE_EXTENSION != "fq.gz" ]
  then
    READ1_INPUT=$(echo "$BASE_NAME$READ1_SUFFIX" | sed -r "s/$FILE_EXTENSION/fq.gz/g")
    cp "$BASE_NAME$READ1_SUFFIX" "$READ1_INPUT"
  else
    READ1_INPUT="$BASE_NAME$READ1_SUFFIX"
  fi

  READ_ARG="--fastq $READ1_INPUT"
  if [ -n "$READ2_SUFFIX" ]
  then
    if [ $FILE_EXTENSION != "fq.gz" ]
    then
      READ2_INPUT=$(echo "$BASE_NAME$READ2_SUFFIX" | sed -r "s/$FILE_EXTENSION/fq.gz/g")
      cp "$BASE_NAME$READ2_SUFFIX" "$READ2_INPUT"
    else
      READ2_INPUT="$BASE_NAME$READ2_SUFFIX"
    fi

    READ_ARG="$READ_ARG --pairs $READ2_INPUT"
  fi

  source /root/miniconda3/etc/profile.d/conda.sh
  conda activate xengsort

  # Note: can only specify threads for the "classify" method not "index"
  # Benchmark: ~10 minutes to split reads using genome-based index
  # Benchmark: ~7 minutes to split reads using transcriptome-based index
  xengsort classify --index "$INDEX_PATH" $READ_ARG --prefix "$BASE_NAME" --threads 8 > xengsort.log

  conda deactivate

  # Remove unnecessary file copies
  if [ $FILE_EXTENSION != "fq.gz" ]
  then
    rm "$READ1_INPUT"
    if [ -n "$READ2_SUFFIX" ]
    then
      rm "$READ2_INPUT"
    fi
  fi

  # Concatenate desired classifications
  # Note: for a transcriptome-based index, novel transcripts or reads spanning fusion junctions could get classified as "neither"
  # First remove reads we don't care to keep
  rm "${BASE_NAME}-host."* "${BASE_NAME}-ambiguous."*

  if [ -n "$READ2_SUFFIX" ]
  then
    cat "${BASE_NAME}-graft.1.fq.gz" "${BASE_NAME}-both.1.fq.gz" "${BASE_NAME}-neither.1.fq.gz" > "${BASE_NAME}-human${READ1_SUFFIX}"
    cat "${BASE_NAME}-graft.2.fq.gz" "${BASE_NAME}-both.2.fq.gz" "${BASE_NAME}-neither.2.fq.gz" > "${BASE_NAME}-human${READ2_SUFFIX}"
  else
    cat "${BASE_NAME}-"*".fq.gz" > "${BASE_NAME}-human${READ1_SUFFIX}"
  fi

  # At this point, remove all remaining intermediate files from xengsort
  rm -f "${BASE_NAME}-graft."* "${BASE_NAME}-both."* "${BASE_NAME}-neither."*

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

#rm -rf "$INDEX_PATH"
