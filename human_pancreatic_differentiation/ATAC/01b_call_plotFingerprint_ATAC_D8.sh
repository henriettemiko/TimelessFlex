#!/bin/bash

##########
#name:          01b_call_deeptools_ATAC_D8.sh
#description:   call plotFingerprint from deeptools
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 27, 2019
##########


source ../set_variables_hg19.sh

QUALITY_DIR=$OUTPUT_DIR/quality_ATAC
mkdir -p $QUALITY_DIR

qsub -hold_jid "process_ATAC_PE_fastq_*" -N call_plotFingerprint \
    -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_D8.txt -cwd -pe smp 1 \
    -l mem_free=50G $SCRIPT_DIR/ATAC/call_plotFingerprint_D8.sh \
    $OUTPUT_DIR $QUALITY_DIR


exit

