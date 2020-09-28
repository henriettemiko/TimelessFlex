#!/bin/bash

##########
#name:          02a_compute_ChIP_JAMM_peaks.sh
#description:   call peak calling for each time point with JAMM D8
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 8, 2019
##########


source ../set_variables_hg19.sh


for f in $OUTPUT_DIR/ChIP/H3K*/D8
do

    echo $f

    MARK=$(echo $f | rev | cut -d "/" -f 2 | rev)
    TIME=$(echo $f | rev| cut -d "/" -f 1 | rev)

    echo $MARK
    echo $TIME

    JAMM_DIR=$OUTPUT_DIR/ChIP/$MARK/$TIME/JAMM_peaks
    mkdir -p $JAMM_DIR

    qsub -N call_JAMM_${MARK}_${TIME} \
        -hold_jid "process_ChIP_fastq_*" -V -j y \
        -o $JAMM_DIR/call_JAMM.txt -cwd -pe smp 1 \
        -l mem_free=30G,h_vmem=30G -l h_rt=24:00:00 \
        $SCRIPT_DIR/ChIP/call_JAMM.sh $JAMM_DIR $MARK $TIME $OUTPUT_DIR \
        $CHR_SIZES 

done



exit

