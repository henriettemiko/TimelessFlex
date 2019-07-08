#!/bin/bash

##########
#name:          02a_compute_ATAC_MACS2_peaks.sh
#description:   call peak calling with MACS2 for each time point
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh


for f in $OUTPUT_DIR/ATAC/D*
do
    echo $f

    TIME=$(echo $f | rev | cut -d "/" -f 1 | rev)
    echo $TIME

    MACS2_DIR=$OUTPUT_DIR/ATAC/$TIME/MACS2_peaks
    mkdir -p $MACS2_DIR

    qsub -V -j y -o $MACS2_DIR/call_MACS2_ATAC.txt -cwd -pe smp 1 \
        -l mem_free=20G,h_vmem=20G $SCRIPT_DIR/ATAC/call_MACS2_ATAC.sh \
        $MACS2_DIR $TIME $OUTPUT_DIR $MACS2_ORGANISM $CHR_SIZES

done


exit

