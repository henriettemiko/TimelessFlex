#!/bin/bash

##########
#name:          01b_call_plotFingerprint_ChIP.sh
#description:   call plotFingerprint from deeptools for each histone mark
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh


QUALITY_DIR=$OUTPUT_DIR/quality_ChIP
mkdir -p $QUALITY_DIR

cd $QUALITY_DIR

qsub -hold_jid "process_ChIP_fastq_*" -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_H3K27ac.txt -cwd -pe smp 1 \
    -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ChIP/call_plotFingerprint_mark.sh \
    H3K27ac $OUTPUT_DIR $QUALITY_DIR

qsub -hold_jid "process_ChIP_fastq_*" -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_H3K27me3.txt -cwd -pe smp 1 \
    -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ChIP/call_plotFingerprint_mark.sh \
    H3K27me3 $OUTPUT_DIR $QUALITY_DIR

qsub -hold_jid "process_ChIP_fastq_*" -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_H3K4me1.txt -cwd -pe smp 1 \
    -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ChIP/call_plotFingerprint_mark.sh \
    H3K4me1 $OUTPUT_DIR $QUALITY_DIR

qsub -hold_jid "process_ChIP_fastq_*" -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_H3K4me3.txt -cwd -pe smp 1 \
    -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ChIP/call_plotFingerprint_mark.sh \
    H3K4me3 $OUTPUT_DIR $QUALITY_DIR

qsub -hold_jid "process_ChIP_fastq_*" -V -j y \
    -o $QUALITY_DIR/call_plotFingerprint_Input.txt -cwd -pe smp 1 \
    -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ChIP/call_plotFingerprint_mark.sh \
    Input $OUTPUT_DIR $QUALITY_DIR


exit

