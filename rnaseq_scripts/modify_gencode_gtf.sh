#!/bin/bash

GENCODE_DIR=$1
GTF_FILE=$2

cd "$GENCODE_DIR"
echo 'Removing "chr" prefix from all chromosomes'
sed -i -E 's/chr([0-9]+)/\1/g' "$GTF_FILE"
