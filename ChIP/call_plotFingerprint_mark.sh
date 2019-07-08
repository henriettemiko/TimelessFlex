#!/bin/bash

##########
#name:          call_plotFingerprint_mark.sh
#description:   call plotFingerprint from deeptools
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


MARK=$1
OUTPUT_DIR=$2
QUALITY_DIR=$3

echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "current mark: $MARK"
echo "#####################"


FILES=($OUTPUT_DIR/ChIP/$MARK/D*/mapping/*nodup.bam)
echo "${FILES[@]}"

plotFingerprint -b "${FILES[@]}" --labels "D0 rep1" "D0 rep2" "D10 rep1" "D10 rep2" "D2 rep1" "D2 rep2" "D5 rep1" "D5 rep2" "D7 rep1" "D7 rep2" -T "ChIP bam files $MARK" --plotFile $QUALITY_DIR/fingerprint_${MARK}.pdf

exit

