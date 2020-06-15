#!/bin/bash

##########
#name:          call_plotFingerprint_D8.sh
#description:   plotFingerprint for ATAC-seq bam files
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 27, 2019
##########


OUTPUT_DIR=$1
QUALITY_DIR=$2


echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "#####################"


FILES=($OUTPUT_DIR/ATAC/D*/mapping/*_nodup.bam)
echo "${FILES[@]}"

plotFingerprint -b "${FILES[@]}" --labels "D0 rep1" "D0 rep2" "D10 rep1" \
    "D10 rep2" "D2 rep1" "D2 rep2" "D5 rep1" "D5 rep2" "D7 rep1" "D7 rep2" \
    "D8 rep1" "D8 rep2" -T "ATAC-seq bam files" \
    --plotFile $QUALITY_DIR/fingerprint_ATAC.pdf


exit

