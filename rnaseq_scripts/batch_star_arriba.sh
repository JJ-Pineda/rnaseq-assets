#!/bin/bash

# If only STAR index creation desired, run script with no specified arguments
# If performing downstream transcriptome assembly, specify "Basic" for "TWO_PASS_MODE"
FASTQ_DIR=$1
READ1_SUFFIX=$2
READ2_SUFFIX=$3
TWO_PASS_MODE=$4   # "None" or "Basic"

SECONDS=0

STAR_INDEX="/root/indexes/star/grch38"
REFERENCE_DIR="/root/gencode_references/grch38"
ASSEMBLY="${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa.gz"
ANNOTATION="${REFERENCE_DIR}/gencode.v46.primary_assembly.annotation.gtf.gz"
ARRIBA_PATH=$(echo "$(type -p arriba)" | sed -r "s/\/arriba$//g")
BLACKLIST_TSV="${ARRIBA_PATH}/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
KNOWN_FUSIONS_TSV="${ARRIBA_PATH}/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
PROTEIN_DOMAINS_GFF3="${ARRIBA_PATH}/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"

# Tell bash to abort on error
set -o pipefail
set -e -u

# Gunzip assembly, annotation, blacklisted fusions, and known fusions files
gunzip -k "$ASSEMBLY" "$ANNOTATION" "$BLACKLIST_TSV" "$KNOWN_FUSIONS_TSV"
ASSEMBLY=$(echo "$ASSEMBLY" | sed -r "s/.gz//g")
ANNOTATION=$(echo "$ANNOTATION" | sed -r "s/.gz//g")
BLACKLIST_TSV=$(echo "$BLACKLIST_TSV" | sed -r "s/.gz//g")
KNOWN_FUSIONS_TSV=$(echo "$KNOWN_FUSIONS_TSV" | sed -r "s/.gz//g")

# Build STAR index if needed
if [ -z "$(ls -A $STAR_INDEX)" ]
then
  echo "No STAR index detected --> building STAR index"
  STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "$STAR_INDEX" --genomeFastaFiles "$ASSEMBLY" --sjdbGTFfile "$ANNOTATION" --sjdbOverhang 100
  echo "Finished building STAR index"
fi

cd "$FASTQ_DIR"

# Create directories for STAR and Arriba
STAR_OUT_DIR="star_alignment"
ARRIBA_OUT_DIR="arriba"
mkdir $STAR_OUT_DIR $ARRIBA_OUT_DIR

# Adapted from "run_arriba.sh" script provided with Arriba installation
# Add "--twopassMode Basic \" if you plan to do transcriptome assembly with StringTie
echo "Will use \"--twopassMode $TWO_PASS_MODE\" for STAR alignment"

READ1_FILES=$(ls *$READ1_SUFFIX)

for f in $READ1_FILES
do
  BASE_NAME="${f//$READ1_SUFFIX/}"
  echo "Processing reads for $BASE_NAME"

  # Check for paired-end file
  if [ -z "$READ2_SUFFIX" ]
  then
    READ_FILES="$BASE_NAME$READ1_SUFFIX"
  else
  	READ_FILES="$BASE_NAME$READ1_SUFFIX $BASE_NAME$READ2_SUFFIX"
  fi

  STAR \
    --runThreadN 8 \
    --genomeDir "$STAR_INDEX" \
    --genomeLoad NoSharedMemory \
    --readFilesIn $READ_FILES \
    --readFilesCommand gunzip -c \
    --twopassMode "$TWO_PASS_MODE" \
    --quantMode TranscriptomeSAM \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outSAMstrandField intronMotif \
    --outBAMcompression 0 \
    --outFilterMultimapNmax 50 \
    --outFileNamePrefix "${STAR_OUT_DIR}/${BASE_NAME}_" \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM HardClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreSeparation 1 \
    --chimSegmentReadGapMax 3 \
    --chimMultimapNmax 50 |

  tee "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.out.bam" |

  arriba \
	-x /dev/stdin \
	-o "${ARRIBA_OUT_DIR}/${BASE_NAME}_fusions.tsv" \
	-O "${ARRIBA_OUT_DIR}/${BASE_NAME}_fusions.discarded.tsv" \
	-a "$ASSEMBLY" \
	-g "$ANNOTATION" \
	-b "$BLACKLIST_TSV" \
	-k "$KNOWN_FUSIONS_TSV" \
	-t "$KNOWN_FUSIONS_TSV" \
	-p "$PROTEIN_DOMAINS_GFF3"

  # Filter and sort BAM files
  if [[ "$TWO_PASS_MODE" == "Basic"  ]]
  then
    echo "Filtering and sorting genomics coordinate BAM file"
    sambamba view -F "not chimeric" -f bam --compression-level=0 "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.out.bam" |
    sambamba sort --compression-level=6 -o "${STAR_OUT_DIR}/${BASE_NAME}_Filtered_Sorted.out.bam" /dev/stdin
  fi

  echo "Filtering and sorting transcriptomics coordinate BAM file"
  sambamba view -F "not chimeric" -f bam --compression-level=0 "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.toTranscriptome.out.bam" |
  sambamba sort --compression-level=6 -o "${STAR_OUT_DIR}/${BASE_NAME}_Filtered_Sorted.toTranscriptome.out.bam" /dev/stdin

  # Free up disk space
  rm "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.out.bam" "${STAR_OUT_DIR}/${BASE_NAME}_Aligned.toTranscriptome.out.bam"

  duration=$SECONDS
  echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds have elapsed."
done

# Remove decompressed assembly/annotation files
rm "$ASSEMBLY" "$ANNOTATION" "$BLACKLIST_TSV" "$KNOWN_FUSIONS_TSV"
