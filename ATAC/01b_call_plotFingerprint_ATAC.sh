#!/bin/bash

##########
#name:          01b_call_deeptools_ATAC.sh
#description:   call plotFingerprint from deeptools
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh

QUALITY_DIR=$OUTPUT_DIR/quality_ATAC
mkdir -p $QUALITY_DIR

qsub -V -j y -o $QUALITY_DIR/call_plotFingerprint.txt -cwd -pe smp 1 \
    -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ATAC/call_plotFingerprint.sh \
    $OUTPUT_DIR $QUALITY_DIR


exit

