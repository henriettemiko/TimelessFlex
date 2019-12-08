#!/bin/bash

##########
#name:          01b_call_plotFingerprint_ChIP_D8.sh
#description:   call plotFingerprint from deeptools for each histone mark
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 27, 2019
##########


source ../set_variables_hg19.sh


QUALITY_DIR=$OUTPUT_DIR/quality_ChIP
mkdir -p $QUALITY_DIR

for f in $OUTPUT_DIR/ChIP/*
do

    MARK=$(echo $f | rev | cut -d "/" -f 1 | rev)
    echo $MARK


qsub -hold_jid "process_ChIP_fastq_*" -N call_plotFingerprint_${MARK} \
    -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_${MARK}_D8.txt -cwd -pe smp 1 \
    -l mem_free=50G $SCRIPT_DIR/ChIP/call_plotFingerprint_mark_D8.sh \
    $MARK $OUTPUT_DIR $QUALITY_DIR

done


exit


