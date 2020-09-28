#!/bin/bash

##########
#name:          02a_compute_ATAC_MACS2_peaks.sh
#description:   call peak calling with MACS2 for each time point
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          April 2, 2020
##########


source ../set_variables_mm10.sh


for f in $OUTPUT_DIR/ATAC/*
do
    echo $f

    TIME=$(echo $f | rev | cut -d "/" -f 1 | rev)
    echo $TIME

    MACS2_DIR=$OUTPUT_DIR/ATAC/$TIME/MACS2_peaks
    mkdir -p $MACS2_DIR

    qsub -hold_jid "process_ATAC_PE_fastq_*" -N call_MACS2_ATAC_${TIME} \
        -V -j y -o $MACS2_DIR/call_MACS2_ATAC_${TIME}.txt -cwd -pe smp 1 \
        -l mem_free=20G,h_vmem=20G $SCRIPT_DIR/ATAC/call_MACS2_ATAC.sh \
        $MACS2_DIR $TIME $OUTPUT_DIR $MACS2_ORGANISM $CHR_SIZES

done


exit

