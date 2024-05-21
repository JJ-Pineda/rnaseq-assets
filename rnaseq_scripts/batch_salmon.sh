#!/bin/bash

FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3  # Leave blank for single-end reads

SECONDS=0

# Tell bash to abort on error
set -o pipefail
set -e -u

# Check for Salmon index
ASSETS_DIR="/Users/javier/CompBioAssets"
cd "${ASSETS_DIR}/salmon_indices/grch38"

# Create Salmon index as needed
if [ ! -d "salmon_index" ]; then
  echo "No salmon index detected --> building salmon index"

  # Extract names of genome targets for decoys construction
  cd "${ASSETS_DIR}/gencode_references/grch38/"
  
  grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
  sed -i.bak -e 's/>//g' decoys.txt
  
  # Concatenate transcripts and genome targets (transcripts MUST come first)
  cat gencode.v45.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
  
  # Build Salmon index (takes ~30 minutes and creates ~15gb worth of stuff --> maybe re-do with every batch)
  # --gencode flag is for removing extra metadata in the target header separated by | from the gencode reference
  # --threads parameter specifies number of threads to use for index creation
  # -k parameter specifies k-mer size (defaulted to 31 which is best for reads >=75bp)
  cd "${ASSETS_DIR}/salmon_indices/grch38/"

  salmon index \
    -t ../../gencode_references/grch38/gentrome.fa.gz \
    --decoys ../../gencode_references/grch38/decoys.txt \
    --threads 12 \
    --index salmon_index \
    --gencode

  echo "Created Salmon index from GRCh38 transcriptome using genome for decoys"

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
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
    --index="${ASSETS_DIR}/salmon_indices/grch38/salmon_index" \
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
